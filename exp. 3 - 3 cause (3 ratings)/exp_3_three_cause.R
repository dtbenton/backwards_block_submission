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
D_tall$condition = as.factor(D_tall$condition)
D_tall$condition_order = as.factor(D_tall$condition_order)

D_tall = D_tall[order(D_tall$ID),]

# ADD A CONDITION NAME COLUMN
D_tall$condition_names = as.factor(rep(1:4, each = 9, times = 24))
D_tall$condition_names = revalue(x = as.factor(D_tall$condition_names), 
                                 c("1" = "BB", "2"="IS", "3" = "1C", 
                                   "4" = "2C"))


# ADD A CONDITION ORDER COLUMN
D_tall$condition_order = revalue(x = as.factor(D_tall$condition_order), 
                     c("1" = "1234", "2"="2413", "3"="3142", "4"="4321"))

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
D_tall = D_tall[,c(1,2,3,5,6,7,4)]

########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################
# NORMALITY CHECK
# plot norm plots for each condition
par(mfrow=c(2,2)) 
for (ii in c("BB","IS","1C","2C"))hist(D_tall$measure[D_tall$condition_names==ii], breaks=5)
par(mfrow=c(1,1))

# get p-values for multi-variate norm test
shapiro.ps = rep(0,4)
for(i in c("BB","IS","1C","2C")) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition_names==i])
  shapiro.ps[i] = shap.calc$p.value
}


# EQUAL VARIANCE CHECK
#box plots
boxplot(D_tall$measure~D_tall$condition_names) 

# formal test of equal variance
leveneTest(D_tall$measure, as.factor(D_tall$condition_names), center=median) # used 'median' because it's a better measure of central tendency given the non-normality



# ASSUMPTION CHECK SUMMARY
# Based on the analyses above, there is a clear violation of the multi-variate normality and 
# the homoskedasticity assumptions. 
# Violations were indicated by a p-value of less than .005 for 22 of the 24 tests.
# Conventional parametric tests, therefore, are not appropriate, and so subsequent
# confidence intervals will be estimated using boostrapping and p-values will be
# obtained using permutation testing. Planned comparisons were also conducted using
# permutation tests.

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
prelim_analysis = lme(measure~(sex+condition_names+condition_order+objects+phase)^5,
                      random=~1|ID,
                      data=D_tall)
anova.lme(prelim_analysis)


########################################################
####          CONTROL CONDITITION ANALYSES          ####
########################################################

# NOTE THAT FORMAL ANALYSIS WERE NOT INCLUDED IN THE MANUSCRIPT FOR EXPERIMENT 3.
# HOWEVER, THE CODE WILL BE KEPT HERE IN CASE I'M REQUIRED TO REPORT THEM IN THE 
# REVIEW OF THE MS.
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

global_boot("1C","Pre","A")

# Amid:

global_boot("1C","Mid","A")

# Apost:

global_boot("1C","Post","A")


# Apre vs Amid
perm_func("1C","1C","Pre","Mid","A","A")

global_boot_2("1C","1C","Pre","Mid","A","A") 


# Apre vs Apost
perm_func("1C","1C","Pre","Post","A","A")

global_boot_2("1C","1C","Pre","Post","A","A")


# Amid vs Apost
perm_func("1C","1C","Mid","Post","A","A")

global_boot_2("1C","1C","Mid","Post","A","A")




#### B RATINGS AND MEASURES ####
# Bpre:

global_boot("1C","Pre","B")

# Bmid:

global_boot("1C","Mid","B")

# Bpost:

global_boot("1C","Post","B")




# Bpre vs Bmid
perm_func("1C","1C","Pre","Mid","B","B")

global_boot_2("1C","1C","Pre","Mid","B","B") 


# Bpre vs Bpost
perm_func("1C","1C","Pre","Post","B","B")

global_boot_2("1C","1C","Pre","Post","B","B")


# Bmid vs Bpost
perm_func("1C","1C","Mid","Post","B","B")

global_boot_2("1C","1C","Mid","Post","B","B")


#### C RATINGS AND MEASURES ####
# Cpre:

global_boot("1C","Pre","C")

# Cmid:

global_boot("1C","Mid","C")

# Cpost:

global_boot("1C","Post","C")




# Cpre vs Cmid
perm_func("1C","1C","Pre","Mid","C","C")

global_boot_2("1C","1C","Pre","Mid","C","C") 


# Cpre vs Cpost
perm_func("1C","1C","Pre","Post","C","C")

global_boot_2("1C","1C","Pre","Post","C","C")


# Cmid vs Cpost
perm_func("1C","1C","Mid","Post","C","C")

global_boot_2("1C","1C","Mid","Post","C","C")

#######################
# TWO-CAUSE CONDITION #
#######################
#### A RATINGS AND MEASURES ####
# Apre:

global_boot("2C","Pre","A")

# Amid:

global_boot("2C","Mid","A")

# Apost:

global_boot("2C","Post","A")


# Apre vs Amid
perm_func("2C","2C","Pre","Mid","A","A")

global_boot_2("2C","2C","Pre","Mid","A","A") 


# Apre vs Apost
perm_func("2C","2C","Pre","Post","A","A")

global_boot_2("2C","2C","Pre","Post","A","A")


# Amid vs Apost
perm_func("2C","2C","Mid","Post","A","A")

global_boot_2("2C","2C","Mid","Post","A","A")







#### B RATINGS AND MEASURES ####
# Bpre:

global_boot("2C","Pre","B")

# Bmid:

global_boot("2C","Mid","B")

# Bpost:

global_boot("2C","Post","B")




# Bpre vs Bmid
perm_func("2C","2C","Pre","Mid","B","B")

global_boot_2("2C","2C","Pre","Mid","B","B") 


# Bpre vs Bpost
perm_func("2C","2C","Pre","Post","B","B")

global_boot_2("2C","2C","Pre","Post","B","B")



# Bmid vs Bpost
perm_func("2C","2C","Mid","Post","B","B")

global_boot_2("2C","2C","Mid","Post","B","B")




#### C RATINGS AND MEASURES ####
# Cpre:

global_boot("2C","Pre","C")

# Cmid:

global_boot("2C","Mid","C")

# Cpost:

global_boot("2C","Post","C")




# Cpre vs Cmid
perm_func("2C","2C","Pre","Mid","C","C")

global_boot_2("2C","2C","Pre","Mid","C","C") 


# Cpre vs Cpost
perm_func("2C","2C","Pre","Post","C","C")

global_boot_2("2C","2C","Pre","Post","C","C")



# Cmid vs Cpost
perm_func("2C","2C","Mid","Post","C","C")

global_boot_2("2C","2C","Mid","Post","C","C")


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
# Mean: 54.58; 95%CI[48.62,60.54]
global_boot("IS","Pre","A")

# Amid:
# Mean: 65.83; 95%CI[57.94,73.73]
global_boot("IS","Mid","A")

# Apost:
# Mean: 10.21; 95%CI[2.23,18.19]
global_boot("IS","Post","A")


# Apre vs Amid
perm_func("IS","IS","Pre","Mid","A","A")
# -11.25000   1.00000   0.00000   0.89700   0.09925
global_boot_2("IS","IS","Pre","Mid","A","A") 
# -11.250000 -21.112594  -1.387406

# Apre vs Apost
perm_func("IS","IS","Pre","Post","A","A")
# 44.375  0.000  1.000  0.000  1.000
global_boot_2("IS","IS","Pre","Post","A","A")
# 44.37500 34.41221 54.33779

# Amid vs Apost
perm_func("IS","IS","Mid","Post","A","A")
# 55.625  0.000  1.000  0.000  1.000
global_boot_2("IS","IS","Mid","Post","A","A")
# 55.6250 44.4737 66.7763






#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 52.08333; 95%CI[44.25578,59.91089]
global_boot("IS","Pre","B")
# 52.08333 44.25578 59.91089

# Bmid:
# Mean: 67.91667; 95%CI[59.40407,76.42926]
global_boot("IS","Mid","B")
# 67.91667 59.40407 76.42926

# Bpost:
# Mean: 98.75000; 95%CI[97.55206,99.94794]
global_boot("IS","Post","B")
# 98.75000 97.55206 99.94794



# Bpre vs Bmid
perm_func("IS","IS","Pre","Mid","B","B")
# -15.83333   1.00000   0.00000   0.96675   0.03200
global_boot_2("IS","IS","Pre","Mid","B","B") 
# -15.83333 -27.53501  -4.13166

# Bpre vs Bpost
perm_func("IS","IS","Pre","Post","B","B")
# -46.66667   1.00000   0.00000   1.00000   0.00000
global_boot_2("IS","IS","Pre","Post","B","B")
# -46.66667 -54.62106 -38.71227


# Bmid vs Bpost
perm_func("IS","IS","Mid","Post","B","B")
# -30.83333   1.00000   0.00000   1.00000   0.00000
global_boot_2("IS","IS","Mid","Post","B","B")
# -30.83333 -39.43753 -22.22914



#### C RATINGS AND MEASURES ####
# Cpre:
# Mean: 52.08; 95%CI[44.25578,59.91089]
global_boot("IS","Pre","B")
# 52.08333 44.25578 59.91089

# Cmid:
# Mean: 67.91667; 95%CI[59.40407,76.42926]
global_boot("IS","Mid","B")
# 67.91667 59.40407 76.42926

# Cpost:
# Mean: 48.33333; 95%CI[39.38302,57.28365]
global_boot("IS","Post","C")
# 48.33333 39.38302 57.28365



# Cpre vs Cmid
perm_func("IS","IS","Pre","Mid","C","C")
# 0.4166667 0.9492500 0.0370000 0.4725000 0.5207500
global_boot_2("IS","IS","Pre","Mid","C","C") 
# 0.4166667 -10.8830480  11.7163813

# Cpre vs Cpost
perm_func("IS","IS","Pre","Post","C","C")
# 4.37500 0.60675 0.38125 0.29950 0.69575
global_boot_2("IS","IS","Pre","Post","C","C")
# 4.375000 -7.359509 16.109509


# Cmid vs Cpost
perm_func("IS","IS","Mid","Post","C","C")
# 3.958333 0.658000 0.327500 0.328000 0.664500
global_boot_2("IS","IS","Mid","Post","C","C")
# 3.958333 -8.177740 16.094407




#### BAYES FACTOR TO COMPARE PRE- AND MID RATINGS OF OBJECT B IN THE IS CONDITION ####
IS_subset_2 = subset(IS_subset, ! phase %in% c("Post"))
IS_subset_3 = subset(IS_subset_2, ! objects %in% c("B","C"))
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
# Mean: 47.95833; 95%CI[41.68105,54.23562]
global_boot("BB","Pre","A")
# 47.95833 41.68105 54.23562

# Amid:
# Mean: 60.66667; 95%CI[51.60960,69.72373]
global_boot("BB","Mid","A")
# 60.66667 51.60960 69.72373

# Apost:
# Mean: 92.95833; 95%CI[84.18414,101.73252]
global_boot("BB","Post","A")
# 92.95833  84.18414 101.73252

# Apre vs Amid
perm_func("BB","BB","Pre","Mid","A","A")
# -12.70833   1.00000   0.00000   0.93325   0.06350
global_boot_2("BB","BB","Pre","Mid","A","A") 
# -12.708333 -23.643729  -1.772937

# Apre vs Apost
perm_func("BB","BB","Pre","Post","A","A")
# -45   1   0   1   0
global_boot_2("BB","BB","Pre","Post","A","A")
# -45.00000 -55.80508 -34.19492

# Amid vs Apost
perm_func("BB","BB","Mid","Post","A","A")
# -32.29167   1.00000   0.00000   1.00000   0.00000
global_boot_2("BB","BB","Mid","Post","A","A")
# -32.29167 -44.86835 -19.71499


#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: ; 95%CI[,]
global_boot("BB","Pre","B")
#45.41667 40.49397 50.33937

# Bmid:
# Mean: 57.70833; 95%CI[47.99414,67.42252]
global_boot("BB","Mid","B")
# 57.70833 47.99414 67.42252

# Bpost:
# Mean: 45.20833; 95%CI[34.07616,56.34051]
global_boot("BB","Post","B")
#  45.20833 34.07616 56.34051



# Bpre vs Bmid
perm_func("BB","BB","Pre","Mid","B","B")
# -12.29167   1.00000   0.00000   0.92550   0.07125
global_boot_2("BB","BB","Pre","Mid","B","B") 
# -12.291667 -23.180229  -1.403104

# Bpre vs Bpost
perm_func("BB","BB","Pre","Post","B","B")
# 0.2083333 0.9760000 0.0120000 0.4815000 0.5120000
global_boot_2("BB","BB","Pre","Post","B","B")
# 0.2083333 -11.9132980  12.3299647


# Bmid vs Bpost
perm_func("BB","BB","Mid","Post","B","B")
# 12.50000  0.14775  0.85075  0.06675  0.93250
global_boot_2("BB","BB","Mid","Post","B","B")
# 12.500000 -2.068473 27.068473


#### C RATINGS AND MEASURES ####
# Cpre:
# Mean: 46.66667; 95%CI[38.38119,54.95214]
global_boot("BB","Pre","C")
# 46.66667 38.38119 54.95214

# Cmid:
# Mean: 44.16667; 95%CI[35.19274,53.14059]
global_boot("BB","Mid","C")
# 44.16667 35.19274 53.14059

# Cpost:
# Mean: 45.41667; 95%CI[35.35893,55.47440]
global_boot("BB","Post","C")
# 45.41667 35.35893 55.47440




# Cpre vs Cmid
perm_func("BB","BB","Pre","Mid","C","C")
# 2.50000 0.75575 0.23000 0.37875 0.61425
global_boot_2("BB","BB","Pre","Mid","C","C") 
# 2.500000 -9.690093 14.690093

# Cpre vs Cpost
perm_func("BB","BB","Pre","Post","C","C")
# 1.25000 0.86675 0.11950 0.43250 0.56175
global_boot_2("BB","BB","Pre","Post","C","C")
# 1.25000 -11.66125  14.16125


# Cmid vs Cpost
perm_func("BB","BB","Mid","Post","C","C")
# -1.25000  1.00000  0.00000  0.55300  0.43825
global_boot_2("BB","BB","Mid","Post","C","C")
# -1.25000 -14.58547  12.08547



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
  facet_wrap(condition_names, labeller = label_wrap_gen(multi_line=FALSE)) + # scales='free' ensures that each blot has x labels
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


# FOR THE BB CONDITION ONLY
condition_barplot = ggplot(BB_subset, aes(objects, measure, fill = phase)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge", colour = "black") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  ylab("ratings (scale: 0-100)") + # change the label of the y-axis
  facet_wrap(~ID, scales = 'free') + # scales='free' ensures that each blot has x labels
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

