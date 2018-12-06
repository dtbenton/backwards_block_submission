########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 3 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
# load all relevant libraries:
library(lme4)
library(nlme)
library(boot)
library(car) 
library(reshape2)
library(ggplot2)
library(ez)
library(plyr)
library(ggsignif)
library(lsr)
library(sjmisc)
library(sjstats)
options(scipen=9999)


# DATA CLEAN UP AND RESTRUCTURING #
D = read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
D = D[c(1:24),]

D_tall = reshape(D, varying = 4:39, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

D_tall$measure = as.numeric(D_tall$measure)
D_tall$sex = as.factor(D_tall$sex)
D_tall$condition_order = as.factor(D_tall$condition_order)
D_tall$condition = as.factor(D_tall$condition)

D_tall = D_tall[order(D_tall$ID),]

# ADD A CONDITION NAME COLUMN
D_tall$condition_names = as.factor(rep(1:4, each = 9, times = 24))
D_tall$condition_names = revalue(x = as.factor(D_tall$condition_names), 
                                 c("1" = "BB", "2"="IS", "3" = "1C", 
                                   "4" = "2C"))

# ADD A 'PHASE' COLUMN
D_tall$phase = as.factor(rep(1:3, each = 3, times = 96))
D_tall$phase = revalue(x = as.factor(D_tall$phase), 
                       c("1" = "Pre", "2"="Mid", "3" = "Post"))

# RENAME SEX COLUMN
D_tall$sex = revalue(x = as.factor(D_tall$sex), 
                     c("1" = "M", "2"="F"))

# OBJECT COLUMN
D_tall$objects = as.factor(rep(1:3, times = 288))
D_tall$objects = revalue(x = as.factor(D_tall$objects), 
                         c("1" = "A", "2"="B", "3"="C"))
# REORDER COLUMNS'
D_tall$condition = NULL
D_tall$row.names = NULL
D_tall = D_tall[,c(1,2,4,5,6,3)]


# 1,2,3,

########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################
# NORMALITY CHECK
D_tall$norm_col = rep(1:24, times=24)
hist(D_tall$measure[D_tall$norm_col==ii], breaks=5)
shapiro.ps = rep(0,24)
for(i in 1:16) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$norm_col==i])
  shapiro.ps[i] = shap.calc$p.value
}


# EQUAL VARIANCE CHECK
#box plots
boxplot(D_tall$measure~D_tall$norm_col) 

# formal test of equal variance
leveneTest(D_tall$measure, as.factor(D_tall$norm_col), center=median) # used 'median' because it's a better measure of central tendency given the non-normality



# ASSUMPTION CHECK SUMMARY
# Based on the analyses above, there is a clear violation of the multi-variate normality and 
# the homoskedasticity assumptions. 
# Violations were indicated by a p-value of less than .005 for 22 of the 24 tests.
# Conventional parametric tests, therefore, are not appropriate, and so subsequent
# confidence intervals will be estimated using boostrapping and p-values will be
# obtained using permutation testing. Planned comparisons were also conducted using
# Wilcoxon paired sign-ranked tests.

########################################################
########################################################
########################################################
#############                              #############
#############            Models            #############
#############                              #############
########################################################
########################################################
########################################################

########################
### GLOBAL FUNCTIONS ###
########################
#  PERMUTATION FUNCTION
perm_func = function(s1,s2,p1,p2,o1,o2){ # s=condition name; p=phase; o=objects
  set.seed(2018)                         # do NOT forget to put the arguments in quote.
  b = rep(0,4000) 
  for(i in 1:4000){
    x = sample(D_tall$measure)
    dif = x[D_tall$condition_names==s1 & D_tall$phase==p1 & D_tall$objects==o1] - 
      x[D_tall$condition_names==s2 & D_tall$phase==p2 & D_tall$objects==o2]
    b[i] = mean(dif)
  }
  
  # compute the actual difference beween BB pre and BB post
  bb_diff = mean(D_tall$measure[D_tall$condition_names==s1 & D_tall$phase==p1 & D_tall$objects==o1]-
                   D_tall$measure[D_tall$condition_names==s2 & D_tall$phase==p2 & D_tall$objects==o2])
  
  # 1- and 2-tailed p-values
  c(bb_diff, sum(abs(b) > bb_diff)/4000, sum(abs(b) < bb_diff)/4000,
    sum(b > bb_diff)/4000, sum(b < bb_diff)/4000)
}



# BOOTSTRAP FUNCTION 
# Single-factor bootstrap function
global_boot = function(s,p,o){
  set.seed(2018)
  boot_fit = function(data,b,formula){ 
    d= data[b,] 
    dif.1 =  mean(d$measure[d$condition_names==s & d$phase==p & d$objects==o], 
                  data=D_tall) 
    return(dif.1)
  }
  
  boot_obj = boot(D_tall, boot_fit, R=4000) 
  c(boot_obj$t0, boot_obj$t0  + 1.96*-sd(boot_obj$t), 
    boot_obj$t0  + 1.96*sd(boot_obj$t))
}




# BOOTSTRAP FUNCTION 
# Condition-difference bootstrap function
global_boot_2 = function(s1,s2,p1,p2,o1,o2){
  set.seed(2018)
  boot_fit = function(data,b,formula){ 
    d= data[b,] 
    dif.1 =  mean(d$measure[d$condition_names==s1 & d$phase==p1 & d$objects==o1], 
                  data=D_tall) - mean(d$measure[d$condition_names==s2 & d$phase==p2 & d$objects==o2], 
                                      data=D_tall)
    return(dif.1)
  }
  
  boot_obj = boot(D_tall, boot_fit, R=4000) 
  c(boot_obj$t0, boot_obj$t0  + 1.96*-sd(boot_obj$t), 
    boot_obj$t0  + 1.96*sd(boot_obj$t))
}




############################
### PRELIMINARY ANALYSES ###
############################

# DETERMINING WHETHER 'SEX' OR 'TEST TRIAL ORDER' INTERACTED WITH ANY OF THE 
# REMAINING FACTORS.
prelim_analysis = lme(measure~(sex+condition_order+condition_names+objects+phase)^5,
                      random=~1|ID,
                      data=D_tall)
anova.lme(prelim_analysis)


########################################################
####          CONTROL CONDITITION ANALYSES          ####
########################################################

#####################################################################################
# CONDITION (1C vs 2C) x OBJECT (A vs B) x PHASE (Pre vs Mid vs Post) OMNIBUS ANOVA #
#####################################################################################
# create a data frame in which the 1C condition is subsetted
one__and_two_cause_subset = subset(D_tall, ! condition_names %in% c("BB","IS")) # creating a smaller
# data set by removing the 
# BB and IS conditions.

# 1C condition
lme_one__and_two_cause_subset = lme(measure~(condition_names+phase+objects)^3, 
                                    random=~1|ID, 
                                    data=one__and_two_cause_subset)

# omnibus ANOVA
anova.lme(lme_one__and_two_cause_subset)


#######################
# ONE-CAUSE CONDITION #
#######################

#### A RATINGS AND MEASURES ####
# Apre:
# Mean: 51.46; 95%CI[45.76,57.16]
global_boot("1C","Pre","A")

# Amid:
# Mean: 98.13; 95%CI[95.84,100.41]
global_boot("1C","Mid","A")

# Apost:
# Mean: 99.38; 95%CI[98.49,100.26]
global_boot("1C","Post","A")


# Apre vs Amid
perm_func("1C","1C","Pre","Mid","A","A")
# -46.66667   1.00000   0.00000   1.00000   0.00000
global_boot_2("1C","1C","Pre","Mid","A","A") 
# -46.66667 -52.77777 -40.55556

# Apre vs Apost
perm_func("1C","1C","Pre","Post","A","A")
# -47.91667   1.00000   0.00000   1.00000   0.00000
global_boot_2("1C","1C","Pre","Post","A","A")
# -47.91667 -53.68569 -42.14764

# Amid vs Apost
perm_func("1C","1C","Mid","Post","A","A")
# -1.25000  1.00000  0.00000  0.54775  0.44600
global_boot_2("1C","1C","Mid","Post","A","A")
# -1.250000 -3.708972  1.208972



#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 50.20833; 95%CI[44.48784,55.92883]
global_boot("1C","Pre","B")

# Bmid:
# Mean: 2.083333; 95%CI[-1.257449,5.424116]
global_boot("1C","Mid","B")

# Bpost:
# Mean: 0.8333333; 95%CI[-0.2966464,1.9633130]
global_boot("1C","Post","B")




# Bpre vs Bmid
perm_func("1C","1C","Pre","Mid","B","B")
# 48.125  0.000  1.000  0.000  1.000
global_boot_2("1C","1C","Pre","Mid","B","B") 
# 48.12500 41.56021 54.68979

# Bpre vs Bpost
perm_func("1C","1C","Pre","Post","B","B")
# 49.375  0.000  1.000  0.000  1.000
global_boot_2("1C","1C","Pre","Post","B","B")
# 49.37500 43.53044 55.21956

# Bmid vs Bpost
perm_func("1C","1C","Mid","Post","B","B")
# 1.250 0.883 0.104 0.455 0.537
global_boot_2("1C","1C","Mid","Post","B","B")
# 1.25000 -2.30575  4.80575



#######################
# TWO-CAUSE CONDITION #
#######################
#### A RATINGS AND MEASURES ####
# Apre:
# Mean: 46.45833; 95%CI[40.94497,51.97169]
global_boot("2C","Pre","A")

# Amid:
# Mean: 96.46; 95%CI[93.03,99.88]
global_boot("2C","Mid","A")

# Apost:
# Mean: 98.33333; 95%CI[96.41997,100.24670]
global_boot("2C","Post","A")


# Apre vs Amid
perm_func("2C","2C","Pre","Mid","A","A")
# -50   1   0   1   0
global_boot_2("2C","2C","Pre","Mid","A","A") 
# -50.00000 -56.49915 -43.50085

# Apre vs Apost
perm_func("2C","2C","Pre","Post","A","A")
# -51.875   1.000   0.000   1.000   0.000
global_boot_2("2C","2C","Pre","Post","A","A")
# -51.87500 -57.71247 -46.03753

# Amid vs Apost
perm_func("2C","2C","Mid","Post","A","A")
# -1.87500  1.00000  0.00000  0.57700  0.41525
global_boot_2("2C","2C","Mid","Post","A","A")
# -1.875000 -5.809238  2.059238






#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 48.33333; 95%CI[43.11416,53.55251]
global_boot("2C","Pre","B")

# Bmid:
# Mean: 44.16667; 95%CI[38.73088,49.60245]
global_boot("2C","Mid","B")

# Bpost:
# Mean: 83.58333; 95%CI[70.77127,96.39539]
global_boot("2C","Post","B")




# Bpre vs Bmid
perm_func("2C","2C","Pre","Mid","B","B")
# 4.166667 0.670500 0.320500 0.331750 0.665000
global_boot_2("2C","2C","Pre","Mid","B","B") 
# 4.166667 -3.337112 11.670445

# Bpre vs Bpost
perm_func("2C","2C","Pre","Post","B","B")
# -35.25   1.00   0.00   1.00   0.00
global_boot_2("2C","2C","Pre","Post","B","B")
# -35.25000 -49.04484 -21.45516


# Bmid vs Bpost
perm_func("2C","2C","Mid","Post","B","B")
# -39.41667   1.00000   0.00000   0.99975   0.00025
global_boot_2("2C","2C","Mid","Post","B","B")
# -39.41667 -53.30272 -25.53061






#####################################################
####          MAIN CONDITITION ANALYSES          ####
#####################################################
################
# IS CONDITION #
################

# create a data frame in which the IS condition is subsetted
IS_subset = subset(D_tall, ! condition_names %in% c("BB","1C", "2C"))

# 1C condition
IS_subset_lme = lme(measure~(phase+objects)^2, 
                    random=~1|ID, 
                    data=IS_subset)

# omnibus ANOVA
anova.lme(IS_subset_lme)



#######################
# PLANNED COMPARISONS #
#######################
#### A RATINGS AND MEASURES ####
# Apre:
# Mean: 47.58333; 95%CI[42.09714,53.06952]
global_boot("IS","Pre","A")

# Amid:
# Mean: 56.12500; 95%CI[49.41984,62.83016]
global_boot("IS","Mid","A")

# Apost:
# Mean: 1.8750000; 95%CI[-0.3801578,4.1301578]
global_boot("IS","Post","A")


# Apre vs Amid
perm_func("IS","IS","Pre","Mid","A","A")
# -8.541667  1.000000  0.000000  0.815000  0.179750
global_boot_2("IS","IS","Pre","Mid","A","A") 
# -8.541667 -17.264271   0.180938

# Apre vs Apost
perm_func("IS","IS","Pre","Post","A","A")
# 45.70833  0.00000  1.00000  0.00000  1.00000
global_boot_2("IS","IS","Pre","Post","A","A")
# 45.70833 39.75218 51.66449

# Amid vs Apost
perm_func("IS","IS","Mid","Post","A","A")
# 54.25  0.00  1.00  0.00  1.00
global_boot_2("IS","IS","Mid","Post","A","A")
# 54.25000 47.10903 61.39097






#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 51.33333; 95%CI[44.04646,58.62021]
global_boot("IS","Pre","B")

# Bmid:
# Mean: 62.79167; 95%CI[55.48057,70.10277]
global_boot("IS","Mid","B")

# Bpost:
# Mean: 97.08333; 95%CI[94.14307,100.02360]
global_boot("IS","Post","B")




# Bpre vs Bmid
perm_func("IS","IS","Pre","Mid","B","B")
# -11.45833   1.00000   0.00000   0.88600   0.11200
global_boot_2("IS","IS","Pre","Mid","B","B") 
# -11.458333 -21.882283  -1.034384

# Bpre vs Bpost
perm_func("IS","IS","Pre","Post","B","B")
# -45.75   1.00   0.00   1.00   0.00
global_boot_2("IS","IS","Pre","Post","B","B")
# -45.75000 -53.57834 -37.92166


# Bmid vs Bpost
perm_func("IS","IS","Mid","Post","B","B")
# -34.29167   1.00000   0.00000   1.00000   0.00000
global_boot_2("IS","IS","Mid","Post","B","B")
# -34.29167 -42.11646 -26.46687


#### BAYES FACTOR TO COMPARE PRE- AND MID RATINGS OF OBJECT B IN THE IS CONDITION ####
IS_subset_2 = subset(IS_subset, ! phase %in% c("Post"))
IS_subset_3 = subset(IS_subset_2, ! objects %in% c("A","C"))
# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=IS_subset_3)
lm.alt = lme(measure~phase, random=~1|ID, data=IS_subset_3)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01



################
# BB CONDITION #
################
# create a data frame in which the IS condition is subsetted
BB_subset = subset(D_tall, ! condition_names %in% c("IS","1C", "2C"))

# 1C condition
BB_subset_lme = lme(measure~(phase+objects)^2, 
                    random=~1|ID, 
                    data=BB_subset)

# omnibus ANOVA
anova.lme(BB_subset_lme)



#######################
# PLANNED COMPARISONS #
#######################
#### A RATINGS AND MEASURES ####
# Apre:
# Mean: 49.3750; 95%CI[42.3325,56.4175]
global_boot("BB","Pre","A")

# Amid:
# Mean: 57.50000; 95%CI[47.82757,67.17243]
global_boot("BB","Mid","A")

# Apost:
# Mean: 98.54167; 95%CI[96.35058,100.73276]
global_boot("BB","Post","A")


# Apre vs Amid
perm_func("BB","BB","Pre","Mid","A","A")
# -8.12500  1.00000  0.00000  0.80975  0.18800
global_boot_2("BB","BB","Pre","Mid","A","A") 
# -8.12500 -20.02575   3.77575

# Apre vs Apost
perm_func("BB","BB","Pre","Post","A","A")
# -49.16667   1.00000   0.00000   1.00000   0.00000
global_boot_2("BB","BB","Pre","Post","A","A")
# -49.16667 -56.52526 -41.80807

# Amid vs Apost
perm_func("BB","BB","Mid","Post","A","A")
# -41.04167   1.00000   0.00000   1.00000   0.00000
global_boot_2("BB","BB","Mid","Post","A","A")
# -41.04167 -51.01843 -31.06491


#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 53.95833; 95%CI[48.22159,59.69507]
global_boot("BB","Pre","B")

# Bmid:
# Mean: 57.70833; 95%CI[50.14964,65.26703]
global_boot("BB","Mid","B")

# Bpost:
# Mean: 49.16667; 95%CI[37.61650,60.71683]
global_boot("BB","Post","B")




# Bpre vs Bmid
perm_func("BB","BB","Pre","Mid","B","B")
# -3.7500  1.0000  0.0000  0.6440  0.3515
global_boot_2("BB","BB","Pre","Mid","B","B") 
# -3.750000 -13.165556   5.665556

# Bpre vs Bpost
perm_func("BB","BB","Pre","Post","B","B")
# 4.791667 0.617250 0.375000 0.304250 0.692250
global_boot_2("BB","BB","Pre","Post","B","B")
# 4.791667 -8.062515 17.645849


# Bmid vs Bpost
perm_func("BB","BB","Mid","Post","B","B")
# 8.541667 0.371750 0.619000 0.184250 0.812250
global_boot_2("BB","BB","Mid","Post","B","B")
# 8.541667 -5.214364 22.297697





#############################################################
# COMPARE POST-RATING OF B BETWEEN THE BB AND IS CONDITIONS #
#############################################################

perm_func("IS","IS","Mid","Post","B","B")
# -34.29167   1.00000   0.00000   1.00000   0.00000
global_boot_2("IS","IS","Mid","Post","B","B")
# -34.29167 -42.11646 -26.46687


perm_func("IS","BB","Post","Post","B","B")
# 47.91667  0.00000  1.00000  0.00000  1.00000

global_boot_2("IS","BB","Post","Post","B","B")
# 47.91667 35.91580 59.91753




################################################################
################################################################
################################################################
#############                                      #############
#############            OMNIBUS FIGURE            #############
#############                                      #############
################################################################
################################################################
################################################################
condition_barplot = ggplot(D_tall, aes(objects, measure, fill = phase)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  facet_wrap(~condition_names, scales = 'free') + # scales='free' ensures that each blot has x labels
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 110)) +
  theme_classic() +
  scale_fill_manual(values = c("white","gray68", "black")) +
  theme(strip.text = element_text(colour = 'black', size = 12)) + # this changes the size and potentially weight of the facet labels
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + # this adds a vertical & horizontal line to each plot
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + # ditto
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")







##############################################################################
##############################################################################
##############################################################################
#############                                                    #############
#############            INDIVIDUAL DIFFERENCE FIGURE            #############
#############                                                    #############
##############################################################################
##############################################################################
##############################################################################
condition_barplot = ggplot(D_tall, aes(objects, measure, fill = phase)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  facet_wrap(ID~condition_names, labeller = label_wrap_gen(multi_line=FALSE)) + # scales='free' ensures that each blot has x labels
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0)) + # ensure that bars hit the x-axis
  coord_cartesian(ylim=c(0, 110)) +
  theme_classic() +
  scale_fill_manual(values = c("white","gray68", "black")) +
  theme(strip.text = element_text(colour = 'black', size = 12)) + # this changes the size and potentially weight of the facet labels
  theme(axis.title=element_text(size="12"),axis.text=element_text(size=12)) + 
  theme(legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6)) +
  theme(legend.text = element_text(size = 12)) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + # this adds a vertical & horizontal line to each plot
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + # ditto
  theme(legend.title=element_blank()) +
  labs(x = "Test trials")


