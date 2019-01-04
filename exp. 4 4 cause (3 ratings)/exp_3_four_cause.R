########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 4 SCRIPT     #############
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
D = D[c(1:20),]

D_tall = reshape(D, varying = 4:51, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

D_tall$measure = as.numeric(D_tall$measure)
D_tall$sex = as.factor(D_tall$sex)
D_tall$condition = as.factor(D_tall$condition)
D_tall$condition_order = as.factor(D_tall$condition_order)

D_tall = D_tall[order(D_tall$ID),]

# ADD A CONDITION NAME COLUMN
D_tall$condition_names = as.factor(rep(1:4, each = 12, times = 20))
D_tall$condition_names = revalue(x = as.factor(D_tall$condition_names), 
                                 c("1" = "BB", "2"="IS", "3" = "1C", 
                                   "4" = "2C"))


# ADD A CONDITION ORDER COLUMN
D_tall$condition_order = revalue(x = as.factor(D_tall$condition), 
                                 c("1" = "1234", "2"="2413", "3"="3142", "4"="4321"))

# ADD A 'PHASE' COLUMN
D_tall$phase = as.factor(rep(1:3, each = 4, times = 80))
D_tall$phase = revalue(x = as.factor(D_tall$phase), 
                       c("1" = "Pre", "2"="Mid", "3" = "Post"))

# RENAME SEX COLUMN
D_tall$sex = revalue(x = as.factor(D_tall$sex), 
                     c("1" = "M", "2"="F"))

# OBJECT COLUMN
D_tall$objects = as.factor(rep(1:4, times = 240))
D_tall$objects = revalue(x = as.factor(D_tall$objects), 
                         c("1" = "A", "2"="B", "3"="C", "4"="D"))
# REORDER COLUMNS'
D_tall$condition = NULL
D_tall$row.names = NULL
D_tall = D_tall[,c(1,2,4,5,6,3)]

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
prelim_analysis = lme(measure~(sex+condition_names+objects+phase)^4,
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
# Mean: 53.5; 95%CI[48.8,62.19]
global_boot("IS","Pre","A")

# Amid:
# Mean: 78.15; 95%CI[71.46,84.84]
global_boot("IS","Mid","A")

# Apost:
# Mean: 8.25; 95%CI[-1.86,18.36]
global_boot("IS","Post","A")


# Apre vs Amid
perm_func("IS","IS","Pre","Mid","A","A")
# -24.650   1.000   0.000   0.996   0.004
global_boot_2("IS","IS","Pre","Mid","A","A") 
# -24.65000 -35.59173 -13.70827

# Apre vs Apost
perm_func("IS","IS","Pre","Post","A","A")
# 45.25  0.00  1.00  0.00  1.00
global_boot_2("IS","IS","Pre","Post","A","A")
# 45.25000 31.86622 58.63378

# Amid vs Apost
perm_func("IS","IS","Mid","Post","A","A")
# 69.9  0.0  1.0  0.0  1.0
global_boot_2("IS","IS","Mid","Post","A","A")
# 69.90000 57.58002 82.21998






#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 51.7500; 95%CI[42.8935,60.6065]
global_boot("IS","Pre","B")
# 51.7500 42.8935 60.6065

# Bmid:
# Mean: 64.40000; 95%CI[51.92206,76.87794]
global_boot("IS","Mid","B")
# 64.40000 51.92206 76.87794

# Bpost:
# Mean: 99.50000; 95%CI[98.52844,100.47156]
global_boot("IS","Post","B")
# 99.50000  98.52844 100.47156



# Bpre vs Bmid
perm_func("IS","IS","Pre","Mid","B","B")
# -12.65000   1.00000   0.00000   0.93725   0.06175
global_boot_2("IS","IS","Pre","Mid","B","B") 
#-12.650000 -27.765826   2.465826

# Bpre vs Bpost
perm_func("IS","IS","Pre","Post","B","B")
# -47.75   1.00   0.00   1.00   0.00
global_boot_2("IS","IS","Pre","Post","B","B")
# -47.75000 -56.67078 -38.829227


# Bmid vs Bpost
perm_func("IS","IS","Mid","Post","B","B")
# -35.1   1.0   0.0   1.0   0.0
global_boot_2("IS","IS","Mid","Post","B","B")
# -35.10000 -47.61257 -22.58743



#### C RATINGS AND MEASURES ####
# Cpre:
# Mean: 44.75000; 95%CI[38.16464,51.33536]
global_boot("IS","Pre","C")
# 44.75000 38.16464 51.33536

# Cmid:
# Mean: 47.25000; 95%CI[38.87933,55.62067]
global_boot("IS","Mid","C")
# 47.25000 38.87933 55.62067

# Cpost:
# Mean: 49.75000; 95%CI[46.53611,52.96389]
global_boot("IS","Post","C")
# 49.75000 46.53611 52.96389



# Cpre vs Cmid
perm_func("IS","IS","Pre","Mid","C","C")
# -2.50000  1.00000  0.00000  0.61650  0.37775
global_boot_2("IS","IS","Pre","Mid","C","C") 
# -2.500000 -13.248049   8.248049

# Cpre vs Cpost
perm_func("IS","IS","Pre","Post","C","C")
# -5.00000  1.00000  0.00000  0.71075  0.28425
global_boot_2("IS","IS","Pre","Post","C","C")
# -5.000000 -12.269478   2.269478


# Cmid vs Cpost
perm_func("IS","IS","Mid","Post","C","C")
# 3-2.5000  1.0000  0.0000  0.6015  0.3920
global_boot_2("IS","IS","Mid","Post","C","C")
# -2.500000 -11.491956   6.491956



#### D RATINGS AND MEASURES ####
# Dpre:
# Mean: 58.25000; 95%CI[51.34072,65.15928]
global_boot("IS","Pre","D")
# 58.25000 51.34072 65.15928

# Dmid:
# Mean: 52.25000; 95%CI[45.73075,58.76925]
global_boot("IS","Mid","D")
# 52.25000 45.73075 58.76925

# Dpost:
# Mean: 51.0000; 95%CI[47.9345,54.0655]
global_boot("IS","Post","D")
# 51.0000 47.9345 54.0655



# Dpre vs Dmid
perm_func("IS","IS","Pre","Mid","D","D")
# 6.00000 0.45250 0.53650 0.22950 0.76525
global_boot_2("IS","IS","Pre","Mid","C","C") 
# -2.500000 -13.248049   8.248049

# Dpre vs Dpost
perm_func("IS","IS","Pre","Post","D","D")
# 7.25000 0.38400 0.60850 0.19425 0.80175
global_boot_2("IS","IS","Pre","Post","D","D")
# 7.2500000 -0.3149281 14.8149281


# Dmid vs Dpost
perm_func("IS","IS","Mid","Post","D","D")
# 1.25000 0.86950 0.11175 0.43800 0.55225
global_boot_2("IS","IS","Mid","Post","D","D")
# 1.250000 -5.948724  8.448724



#### BAYES FACTOR TO COMPARE PRE- AND MID RATINGS OF OBJECT B IN THE IS CONDITION ####
IS_subset_2 = subset(IS_subset, ! phase %in% c("Post"))
IS_subset_3 = subset(IS_subset_2, ! objects %in% c("A","C","D"))
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
# Mean: 50.2500; 95%CI[45.2578,55.2422]
global_boot("BB","Pre","A")
# 50.2500 45.2578 55.2422

# Amid:
# Mean: 71.00000; 95%CI[64.64633,77.35367]
global_boot("BB","Mid","A")
# 71.00000 64.64633 77.35367

# Apost:
# Mean: 99.75000; 95%CI[99.25304,100.24696]
global_boot("BB","Post","A")
# 92.95833  99.25304 100.24696

# Apre vs Amid
perm_func("BB","BB","Pre","Mid","A","A")
# -20.75000   1.00000   0.00000   0.99325   0.00675
global_boot_2("BB","BB","Pre","Mid","A","A") 
# -20.7500 -28.7261 -12.7739

# Apre vs Apost
perm_func("BB","BB","Pre","Post","A","A")
# -49.5   1.0   0.0   1.0   0.0
global_boot_2("BB","BB","Pre","Post","A","A")
# -49.50000 -54.51708 -44.48292

# Amid vs Apost
perm_func("BB","BB","Mid","Post","A","A")
# -28.7500   1.0000   0.0000   0.9995   0.0005
global_boot_2("BB","BB","Mid","Post","A","A")
# -28.75000 -35.12404 -22.37596


#### B RATINGS AND MEASURES ####
# Bpre:
# Mean: 50.90000; 95%CI[44.56488,57.23512]
global_boot("BB","Pre","B")
# 50.90000 44.56488 57.23512

# Bmid:
# Mean: 63.40000; 95%CI[54.73142,72.06858]
global_boot("BB","Mid","B")
# 63.40000 54.73142 72.06858

# Bpost:
# Mean: 50.25000; 95%CI[37.71325,62.78675]
global_boot("BB","Post","B")
#  50.25000 37.71325 62.78675



# Bpre vs Bmid
perm_func("BB","BB","Pre","Mid","B","B")
# -12.5000   1.0000   0.0000   0.9340   0.0635
global_boot_2("BB","BB","Pre","Mid","B","B") 
# -12.500000 -23.271454  -1.728546

# Bpre vs Bpost
perm_func("BB","BB","Pre","Post","B","B")
# 0.6500 0.9370 0.0600 0.4655 0.5330
global_boot_2("BB","BB","Pre","Post","B","B")
# 0.65000 -13.36115  14.66115


# Bmid vs Bpost
perm_func("BB","BB","Mid","Post","B","B")
# 13.15000  0.10725  0.89225  0.05350  0.94600
global_boot_2("BB","BB","Mid","Post","B","B")
# 13.150000 -2.159301 28.459301


#### BAYES FACTOR TO COMPARE POST- AND MID RATINGS OF OBJECT B IN THE IS CONDITION ####
BB_subset_2 = subset(BB_subset, ! phase %in% c("Pre"))
BB_subset_3 = subset(BB_subset_2, ! objects %in% c("A","C"))
# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=BB_subset_3)
lm.alt = lme(measure~phase, random=~1|ID, data=BB_subset_3)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01

#### C RATINGS AND MEASURES ####
# Cpre:
# Mean: 50.0000; 95%CI[44.6435,55.3565]
global_boot("BB","Pre","C")
# 50.0000 44.6435 55.3565

# Cmid:
# Mean: 51.25000; 95%CI[44.41516,58.08484]
global_boot("BB","Mid","C")
# 51.25000 44.41516 58.08484

# Cpost:
# Mean: 51.15000; 95%CI[45.86878,56.43122]
global_boot("BB","Post","C")
# 51.15000 45.86878 56.43122




# Cpre vs Cmid
perm_func("BB","BB","Pre","Mid","C","C")
# -1.25000  1.00000  0.00000  0.55450  0.44025
global_boot_2("BB","BB","Pre","Mid","C","C") 
# -1.250000 -9.957799  7.457799

# Cpre vs Cpost
perm_func("BB","BB","Pre","Post","C","C")
# -1.1500  1.0000  0.0000  0.5595  0.4390
global_boot_2("BB","BB","Pre","Post","C","C")
# -1.1500 -8.6424  6.3424


# Cmid vs Cpost
perm_func("BB","BB","Mid","Post","C","C")
# 0.10000 0.99000 0.00875 0.50175 0.49750
global_boot_2("BB","BB","Mid","Post","C","C")
# 0.100000 -8.604116  8.804116



#### D RATINGS AND MEASURES ####
# Dpre:
# Mean: 53.30000; 95%CI[46.58798,60.01202]
global_boot("BB","Pre","D")
# 53.30000 46.58798 60.01202

# Dmid:
# Mean: 50.50000; 95%CI[46.76551,54.23449]
global_boot("BB","Mid","D")
# 50.50000 46.76551 54.23449

# Dpost:
# Mean: 47.55000; 95%CI[41.99166,53.10834]
global_boot("BB","Post","D")
# 47.55000 41.99166 53.10834




# Dpre vs Dmid
perm_func("BB","BB","Pre","Mid","D","D")
# 2.80000 0.72400 0.27325 0.37125 0.62775
global_boot_2("BB","BB","Pre","Mid","D","D") 
# 2.800000 -4.939667 10.539667

# Dpre vs Dpost
perm_func("BB","BB","Pre","Post","D","D")
# 5.7500 0.4815 0.5090 0.2510 0.7440
global_boot_2("BB","BB","Pre","Post","D","D")
# 5.750000 -2.980983 14.480983


# Dmid vs Dpost
perm_func("BB","BB","Mid","Post","D","D")
# 2.95000 0.74150 0.25675 0.37450 0.62450
global_boot_2("BB","BB","Mid","Post","D","D")
# 2.950000 -3.732969  9.632969

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

