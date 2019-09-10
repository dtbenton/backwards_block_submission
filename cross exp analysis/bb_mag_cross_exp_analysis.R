#############################################################
#############################################################
#############################################################
#############                                   #############
#############      CROSS EXPERIMENT ANLYSIS     #############
#############                                   #############
#############################################################
#############################################################
#############################################################
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
library(BayesFactor)
options(scipen=9999)

# load data
D = read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
D = D[c(1:68),c(1:3)]
D$ID = c(1:68)

# reorder 'D' columns so 'ID' is first
D = D[,c(4,1,2,3)]

# data restructring
D_tall = reshape(D, varying = 3:4, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

# add 'exp' name column
D$Exp = revalue(x = as.factor(D$Exp), 
                c("1" = "Exp1", "2"="Exp2",
                  "3" = "Exp3"))

# rename 'condition' column to 'ratingType'
D_tall$condition = revalue(x = as.factor(D_tall$condition), 
                c("1" = "modifiedRW", "2"="Bayesian"))

# reorder 'Exp' column
D_tall = D_tall[order(D_tall$Exp),]


# remove 'row.names' column
D_tall$row.names = NULL



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
perm_func = function(s1,s2,p1,p2){ # s=condition; p=Experiment
  set.seed(2018)                         # do NOT forget to put the arguments in quote.
  b = rep(0,4000) 
  for(i in 1:4000){
    x = sample(D_tall$measure)
    dif = x[D_tall$condition==s1 & D_tall$Exp==p1] - 
      x[D_tall$condition==s2 & D_tall$Exp==p2]
    b[i] = mean(dif)
  }
  
  # compute the actual difference beween BB pre and BB post
  bb_diff = mean(D_tall$measure[D_tall$condition==s1 & D_tall$Exp==p1]-
                   D_tall$measure[D_tall$condition==s2 & D_tall$Exp==p2])
  
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
    dif.1 =  mean(d$measure[d$condition==s & d$Exp==p], 
                  data=D_tall) 
    return(dif.1)
  }
  
  boot_obj = boot(D_tall, boot_fit, R=4000) 
  c(boot_obj$t0, boot_obj$t0  + 1.96*-sd(boot_obj$t), 
    boot_obj$t0  + 1.96*sd(boot_obj$t))
}


# BOOTSTRAP FUNCTION 
# Condition-difference bootstrap function
global_boot_2 = function(s1,s2,p1,p2){
  set.seed(2018)
  boot_fit = function(data,b,formula){ 
    d= data[b,] 
    dif.1 =  mean(d$measure[d$condition==s1 & d$Exp==p1], 
                  data=D_tall) - mean(d$measure[d$condition==s2 & d$Exp==p2], 
                                      data=D_tall)
    return(dif.1)
  }
  
  boot_obj = boot(D_tall, boot_fit, R=4000) 
  c(boot_obj$t0, boot_obj$t0  + 1.96*-sd(boot_obj$t), 
    boot_obj$t0  + 1.96*sd(boot_obj$t))
}

############################
# MAIN ANALYSIS DIFFERENCE #
############################

# Experiment 1 vs Experiment 2 BB comparison #
onevstwosub = subset(D_tall, ! Exp %in% c("Exp3")) 
onevstwosub = subset(onevstwosub, ! condition %in% c("modifiedRW")) 

# creating a smaller
# data set by removing the 
# BB and IS conditions.

# 1C condition
onevstwosub_lme = lme(measure~Exp, 
                      random=~1|ID,
                      data=onevstwosub)

# omnibus ANOVA
summary(onevstwosub_lme)
anova.lme(onevstwosub_lme)

perm_func("Bayesian","Bayesian","Exp1","Exp2")
global_boot_2("Bayesian","Bayesian","Exp1","Exp2")

# Experiment 1 vs Experiment 2 Bayes' factor
x = D_tall$measure[D_tall$condition=="Bayesian" & D_tall$Exp=="Exp1"]
y = D_tall$measure[D_tall$condition=="Bayesian" & D_tall$Exp=="Exp2"]

BF1 = ttestBF(x=x,y=y,paired=TRUE)
BF1


# Experiment 1 vs Experiment 3 BB comparison #
onevsthreesub = subset(D_tall, ! Exp %in% c("Exp2")) 
onevsthreesub = subset(onevsthreesub, ! condition %in% c("modifiedRW")) 

# creating a smaller
# data set by removing the 
# BB and IS conditions.

# 1C condition
onevsthreesub_lme = lme(measure~Exp, 
                      random=~1|ID,
                      data=onevsthreesub)

# omnibus ANOVA
summary(onevsthreesub_lme)
anova.lme(onevsthreesub_lme)


perm_func("Bayesian","Bayesian","Exp1","Exp3")
global_boot_2("Bayesian","Bayesian","Exp1","Exp3")

# Experiment 1 vs Experiment 2 Bayes' factor
x2 = D_tall$measure[D_tall$condition=="Bayesian" & D_tall$Exp=="Exp1"]
y2 = D_tall$measure[D_tall$condition=="Bayesian" & D_tall$Exp=="Exp4"]

BF2 = ttestBF(x=x2,y=y2,paired=TRUE)
BF2
