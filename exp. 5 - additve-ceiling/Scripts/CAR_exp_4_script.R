########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 4 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
## INITIAL SET UP ##
# import "4_cause_CSV.csv"
D = read.csv(file.choose(), header = TRUE)

D = D[-c(44:1013),]

# convert data from "wide" to "tall" format
D_tall = reshape(D, varying = 3:58, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:1176, direction = "long")

# order by ID
D_tall = D_tall[order(D_tall$ID),]

# remove scientific notation
options(scipen=999)

# install relevant packages
install.packages('lazyeval')
install.packages('ggplot2', dep=TRUE)

# libraries:
library(boot)
library(car)
library(lme4)
library(nlme)
library(pgirmess)
library(ez)
library(ggplot2)


########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

##==> Normality check <==##
# formal test of normality
shapiro.ps = rep(0,56)
for(i in 1:56) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition==i])
  shapiro.ps[i] = shap.calc$p.value
}

shapiro.mat = matrix(NA, nrow=6, ncol=4, byrow = TRUE)
for(ii in 1:24){
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition==ii])
  shapiro.mat[i] = shap.calc$p.value
}


##==> summary of normality check <==##
# Based on the analyses above, there is no evidence of normality for each of the 
# univariate histograms. Violations was indicated by a p-value of less than .005.
# Conventional parametric tests, therefore, are not appropriate, and so subsequent
# confidence intervals will be estimated using boostrapping and p-values will be
# obtained using permutation testing. Planned comparisons were also conducted using
# Wilcoxon paired sign-ranked tests.
#####################################




##==> Equal variance check <==##

# plot the boxplots
boxplot(D_tall$measure~D_tall$condition)


# formal test of equal variance
leveneTest(D_tall$measure, as.factor(D_tall$condition), center=median) # used 'median' because it's a better measure of central tendency given the non-normality


##==> summary of equal variance check <==##
# Based on the analysis above, there is evidence of unequal variance.
# Violations were indicated by p-values that were less than .005.
# Corroborating the 'normality" analyses above, then, subsequent analyses
# will use bootstrapping and parametric tests.
###########################################



##################
## IS CONDITION ##
##################
# A+ events across pre, mid, and post

subD1 = subset(D_tall, ! condition %in% c(1:20,22:24,26:28,30:56))
subD1$condition = factor(subD1$condition)

# run model
lme.fit.aplusesIS = lme(measure~condition, random=~1|ID, data=subD1)

# summary of model 
summary(lme.fit.aplusesIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aplusesIS)

### One-factor Repeated Measures ANOVA ###
install.packages("ez")
library(ez)

# subset data to include pre- and post-ratings of A and B in the 1C condition 
# run repeated-measures ANOVA
repeat_AOV.aplusesIS = ezANOVA(subD1, dv = measure, wid = ID, within = as.factor(condition))
print(repeat_AOV.aplusesIS) 

####################################
## Permutation test for contrasts ##
####################################

    ############
# A+ pre vs. A+ mid
    ############
# create a function that computes the difference between BB pre and BB post
subApluspremid = subset(D_tall, ! condition %in% c(1:20,22:24,26:56))
subApluspremid$condition = factor(subApluspremid$condition)

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subApluspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subApluspremid)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed

# bootstrapped for A+ pre in IS
set.seed(2017)
dif.lm.fit.apluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==21,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.apluspre = boot(D_tall, dif.lm.fit.apluspre, R=5000) 
dif.lm.Bootobj.apluspre
dif.lm.Bootobj.apluspre.confint = 50.2381  + 1.96*c(-3.676017, 3.676017)

# bootstrapped for A+ mid in IS
set.seed(2017)
dif.lm.fit.aplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==25,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aplusmid = boot(D_tall, dif.lm.fit.aplusmid, R=5000) 
dif.lm.Bootobj.aplusmid
dif.lm.Bootobj.aplusmid.confint = 60.47619  + 1.96*c(-4.389916, 4.389916)





    ############
# A+ mid vs. A+ post
    ############

# create a function that computes the difference between BB pre and BB post
subAplusmidpost = subset(D_tall, ! condition %in% c(1:24,26:28,30:56))
subAplusmidpost$condition = factor(subAplusmidpost$condition)

set.seed(2017)
b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusmidpost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b < beta_actual)/5000 #p-value, one-tailed

# bootstrapped for A+ post in IS
set.seed(2017)
dif.lm.fit.apluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==29,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aplusmidpost = boot(D_tall, dif.lm.fit.apluspost, R=5000) 
dif.lm.Bootobj.aplusmidpost
dif.lm.Bootobj.apluspost.confint = 26.57143  + 1.96*c(-5.141666, 5.141666)

# bootstrapped for A+ mid in IS
set.seed(2017)
dif.lm.fit.aplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==25,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aplusmid = boot(D_tall, dif.lm.fit.aplusmid, R=5000) 
dif.lm.Bootobj.aplusmid
dif.lm.Bootobj.aplusmid.confint = 60.47619  + 1.96*c(-4.389916, 4.389916)






  ############
# A+ pre vs. A+ post
  ############
subAplusprepost = subset(D_tall, ! condition %in% c(1:20,22:28,30:56))
subAplusprepost$condition = factor(subAplusprepost$condition)

set.seed(2017)
b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusprepost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b < beta_actual)/5000 #p-value, one-tailed





  ############
# A- pre vs. A- mid
  ############
# A- events across pre, mid, and post

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)

# create a function that computes the difference between BB pre and BB post
subAminuspremid = subset(D_tall, ! condition %in% c(1:21,23:25,27:56))
subAminuspremid$condition = factor(subAminuspremid$condition)

# permutation test for a- pre vs mid ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminuspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAminuspremid)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed

# bootstrapped for A- pre in IS
set.seed(2017)
dif.lm.fit.aminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==22,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminuspre = boot(D_tall, dif.lm.fit.aminuspre, R=5000) 
dif.lm.Bootobj.aminuspre
dif.lm.Bootobj.aminuspre.confint = 46.90476  + 1.96*c(-6.021357, 6.021357)

# bootstrapped for A- mid in IS
set.seed(2017)
dif.lm.fit.aminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==26,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminusmid = boot(D_tall, dif.lm.fit.aminusmid, R=5000) 
dif.lm.Bootobj.aminusmid
dif.lm.Bootobj.aminusmid.confint = 46.19048  + 1.96*c(-4.571999, 4.571999)


}

# bootstrapped for A- post in IS
set.seed(2017)
dif.lm.fit.aminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==30,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminuspost = boot(D_tall, dif.lm.fit.aminuspost, R=5000) 
dif.lm.Bootobj.aminuspost
dif.lm.Bootobj.aminuspost.confint = 71.66667  + 1.96*c(-6.957788, 6.957788)


# permutation test for a- pre vs post ratings
subAminusprepost = subset(D_tall, ! condition %in% c(1:21,23:29,31:56))
subAminusprepost$condition = factor(subAminusprepost$condition)

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusprepost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed


# permutation test for a- mid vs post ratings
subAminusmidpost = subset(D_tall, ! condition %in% c(1:25,27:29,31:56))
subAminusmidpost$condition = factor(subAminusmidpost$condition)

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusmidpost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed






  ############
# B+ pre vs. B+ mid
  ############
# B+ events across pre, mid, and post

# run model
subD3 = subset(D_tall, ! condition %in% c(1:22,24:26,28:30,32:56))
lme.fit.bpluspremidIS = lme(measure~condition, random=~1|ID, data=subD3)

# summary of model 
summary(lme.fit.bpluspremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.bpluspremidIS)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)

# create a function that computes the difference between B+ pre and mid
subBpluspremid = subset(D_tall, ! condition %in% c(1:22,24:26,28:56))
subBpluspremid$condition = factor(subBpluspremid$condition)

# permutation test for B+ pre vs mid ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBpluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBpluspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBpluspremid)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed

# create a function that computes the difference between B+ mid and post
subBplusmidpost = subset(D_tall, ! condition %in% c(1:26,28:30,32:56))
subBplusmidpost$condition = factor(subBplusmidpost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusmidpost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed


# create a function that computes the difference between B+ mid and post
subBplusprepost = subset(D_tall, ! condition %in% c(1:22,24:30,32:56))
subBplusprepost$condition = factor(subBplusprepost$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusprepost)
beta_actual = fixed.effects(lm.fit)[2]


# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed







# bootstrapped for B+ pre in IS
set.seed(2017)
dif.lm.fit.bpluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==23,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bpluspre = boot(D_tall, dif.lm.fit.bpluspre, R=5000) 
dif.lm.Bootobj.bpluspre
dif.lm.Bootobj.bpluspre.confint = 51.42857 + 1.96*c(-3.940875, 3.940875)

# bootstrapped for B+ mid in IS
set.seed(2017)
dif.lm.fit.bplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==27,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bplusmid = boot(D_tall, dif.lm.fit.bplusmid, R=5000) 
dif.lm.Bootobj.bplusmid
dif.lm.Bootobj.bplusmid.confint = 60.47619  + 1.96*c(-4.448825, 4.448825)


# bootstrapped for B+ post in IS
set.seed(2017)
dif.lm.fit.bpluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==31,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bpluspost = boot(D_tall, dif.lm.fit.bpluspost, R=5000) 
dif.lm.Bootobj.bpluspost
dif.lm.Bootobj.bpluspost.confint = 76.19048  + 1.96*c(-5.238196, 5.238196)




############
# B- data
############
# B- events across pre, mid, and post

# run model
subD4 = subset(D_tall, ! condition %in% c(1:23,25:27,29:31,33:56))
lme.fit.bminus = lme(measure~condition, random=~1|ID, data=subD4)

# summary of model 
summary(lme.fit.bminus)

# get ANOVA results for lme model
anova.lme(lme.fit.bminus)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)


## BOOTSTRAPPED ANALYSES FOR THE B- EVENTS ##

# B- pre 
set.seed(2017)
dif.lm.fit.bminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==24,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminuspre = boot(D_tall, dif.lm.fit.bminuspre, R=5000) 
dif.lm.Bootobj.bminuspre
dif.lm.Bootobj.bminuspre.confint = 47.85714 + 1.96*c(-6.125307, 6.125307)


# B- mid 
set.seed(2017)
dif.lm.fit.bminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==28,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminusmid = boot(D_tall, dif.lm.fit.bminusmid, R=5000) 
dif.lm.Bootobj.bminusmid
dif.lm.Bootobj.bminusmid.confint = 46.19048 + 1.96*c(-4.649999, 4.649999)


# B- post 
set.seed(2017)
dif.lm.fit.bminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==32,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminuspost = boot(D_tall, dif.lm.fit.bminuspost, R=5000) 
dif.lm.Bootobj.bminuspost
dif.lm.Bootobj.bminuspost.confint = 32.14286 + 1.96*c(-6.65209, 6.65209)

## PERMUTATION ANALYSES FOR THE B- EVENTS ##

# B- pre vs mid
# create a function that computes the difference between B- pre vs mid
subBminuspremid = subset(D_tall, ! condition %in% c(1:23,25:27,29:56))
subBminuspremid$condition = factor(subBminuspremid$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminuspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminuspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBminuspremid)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed





# B- mid vs post
# create a function that computes the difference between B- pre vs mid
subBminusmidpost = subset(D_tall, ! condition %in% c(1:27,29:31,33:56))
subBminusmidpost$condition = factor(subBminusmidpost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusmidpost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b < beta_actual)/5000 #p-value, one-tailed



# B- pre vs post
# create a function that computes the difference between B- pre vs mid
subBminusprepost = subset(D_tall, ! condition %in% c(1:23,25:31,33:56))
subBminusprepost$condition = factor(subBminusprepost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusprepost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed


##################
## BB CONDITION ##
##################
############
# A+ data
############
# A+ events across pre, mid, and post

# run model
subD5 = subset(D_tall, ! condition %in% c(1:8,10:12,14:16,18:56))
subD5$condition = factor(subD5$condition)
lme.fit.aplusBB = lme(measure~condition, random=~1|ID, data=subD5)

# summary of model 
summary(lme.fit.aplusBB)

# get ANOVA results for lme model
anova.lme(lme.fit.aplusBB)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)


## BOOTSTRAPPED ANALYSES FOR THE A+ EVENTS ##

# A+ pre 
set.seed(2017)
dif.lm.fit.apluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==9,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.apluspre = boot(D_tall, dif.lm.fit.apluspre, R=5000) 
dif.lm.Bootobj.apluspre
dif.lm.Bootobj.apluspre.confint = 76.42857 + 1.96*c(-4.851917, 4.851917)


# A+ mid 
set.seed(2017)
dif.lm.fit.apluspmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==13,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.apluspmid  = boot(D_tall, dif.lm.fit.apluspmid, R=5000) 
dif.lm.Bootobj.apluspmid 
dif.lm.Bootobj.apluspmid.confint = 62.14286 + 1.96*c(-3.639573, 3.639573)

# A+ post 
set.seed(2017)
dif.lm.fit.aplusppost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==17,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aplusppost  = boot(D_tall, dif.lm.fit.aplusppost, R=5000) 
dif.lm.Bootobj.aplusppost 
dif.lm.Bootobj.aplusppost.confint = 84.7619 + 1.96*c(-4.128162, 4.128162)


## PERMUTATION ANALYSIS FOR THE A+ EVENTS ##

# A+ pre vs mid
# create a function that computes the difference between B- pre vs mid
subApluspremid = subset(D_tall, ! condition %in% c(1:8,10:12,14:56))
subApluspremid$condition = factor(subApluspremid$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subApluspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subApluspremid)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed




# A+ pre vs post
# create a function that computes the difference between A+ pre vs post
subAplusprepost = subset(D_tall, ! condition %in% c(1:8,10:16,18:56))
subAplusprepost$condition = factor(subAplusprepost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

hist(b)

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusprepost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed




# A+ mid vs post
# create a function that computes the difference between A+ pre vs post
subAplusmidpost = subset(D_tall, ! condition %in% c(1:12,14:16,18:56))
subAplusmidpost$condition = factor(subAplusmidpost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusmidpost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed



############
# A- data
############
# A+ events across pre, mid, and post

# run model
subD6 = subset(D_tall, ! condition %in% c(1:9,11:13,15:17,19:56))
subD6$condition = factor(subD6$condition)
lme.fit.aminusBB = lme(measure~condition, random=~1|ID, data=subD6)

# summary of model 
summary(lme.fit.aminusBB)

# get ANOVA results for lme model
anova.lme(lme.fit.aminusBB)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)


## BOOTSTRAPPED ANALYSES FOR THE A- EVENTS ##

# A- pre 
set.seed(2017)
dif.lm.fit.aminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==10,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminuspre = boot(D_tall, dif.lm.fit.aminuspre, R=5000) 
dif.lm.Bootobj.aminuspre
dif.lm.Bootobj.aminuspre.confint = 48.09524 + 1.96*c(-5.306536, 5.306536)


# A- mid 
set.seed(2017)
dif.lm.fit.aminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==14,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminusmid  = boot(D_tall, dif.lm.fit.aminusmid, R=5000) 
dif.lm.Bootobj.aminusmid 
dif.lm.Bootobj.aminusmid.confint = 48.09524 + 1.96*c(-4.203841, 4.203841)

# A- post 
set.seed(2017)
dif.lm.fit.aminusppost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==18,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.aminusppost  = boot(D_tall, dif.lm.fit.aminusppost, R=5000) 
dif.lm.Bootobj.aminusppost
dif.lm.Bootobj.aminusppost.confint = 29.09524 + 1.96*c(-6.674694, 6.674694)



## PERMUTATION ANALYSIS FOR THE A+ EVENTS ##
# A- pre vs post
# create a function that computes the difference between A- pre vs post
subAminusprepost = subset(D_tall, ! condition %in% c(1:9,11:17,19:56))
subAminusprepost$condition = factor(subAminusprepost$condition)

# permutation test for A- mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusprepost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.008
sum(b > beta_actual)/5000 #p-value, one-tailed


# A- mid vs post
# create a function that computes the difference between B- pre vs mid
subAminusmidpost = subset(D_tall, ! condition %in% c(1:13,15:17,19:56))
subAminusmidpost$condition = factor(subAminusmidpost$condition)

# permutation test for B+ mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusmidpost)
beta_actual = fixed.effects(lm.fit)[2]

sum(b > beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #1
sum(b > beta_actual)/5000 #p-value, one-tailed





############
# B+ data
############
# A+ events across pre, mid, and post

# run model
subD7 = subset(D_tall, ! condition %in% c(1:10,12:14,16:18,20:56))
lme.fit.bplusBB = lme(measure~condition, random=~1|ID, data=subD7)

# summary of model 
summary(lme.fit.bplusBB)

# get ANOVA results for lme model
anova.lme(lme.fit.bplusBB)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)

mean(D_tall$measure[D_tall$condition==11])
mean(D_tall$measure[D_tall$condition==15])
mean(D_tall$measure[D_tall$condition==19])

## BOOTSTRAPPED ANALYSES FOR THE B+ EVENTS ##

# B+ pre 
set.seed(2017)
dif.lm.fit.bpluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==11,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bpluspre = boot(D_tall, dif.lm.fit.bpluspre, R=5000) 
dif.lm.Bootobj.bpluspre
dif.lm.Bootobj.bpluspre.confint = 73.09524 + 1.96*c(-4.797063, 4.797063)
# 63.69300 82.49748

# B+ mid 
set.seed(2017)
dif.lm.fit.bplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==15,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bplusmid = boot(D_tall, dif.lm.fit.bplusmid, R=5000) 
dif.lm.Bootobj.bplusmid
dif.lm.Bootobj.bplusmid.confint = 59.28571 + 1.96*c(-4.301547, 4.301547)
# 50.85468 67.71674

# B+ post 
set.seed(2017)
dif.lm.fit.bpluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==19,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bpluspost = boot(D_tall, dif.lm.fit.bpluspost, R=5000) 
dif.lm.Bootobj.bpluspost
dif.lm.Bootobj.bpluspost.confint = 45.38095 + 1.96*c(-5.848017, 5.848017)
# 33.91884 56.84306




## PERMUTATION ANALYSIS FOR THE B+ EVENTS ##
# B+ pre vs mid
# create a function that computes the difference between A- pre vs post
subBpluspremid = subset(D_tall, ! condition %in% c(1:10,12:14,16:56))
subBpluspremid$condition = factor(subBpluspremid$condition)

# permutation test for A- mid vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBpluspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBpluspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBpluspremid)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001


# B+ pre vs post
# create a function that computes the difference between A- pre vs post
subBplusprepost = subset(D_tall, ! condition %in% c(1:10,12:18,20:56))
subBplusprepost$condition = factor(subBplusprepost$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusprepost)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001

# B+ pre vs post
# create a function that computes the difference between A- pre vs post
subBplusmidpost = subset(D_tall, ! condition %in% c(1:14,16:18,20:56))
subBplusmidpost$condition = factor(subBplusmidpost$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusmidpost)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001







############
# B- data
############
# B- events across pre, mid, and post

# run model
subD8 = subset(D_tall, ! condition %in% c(1:11,13:15,17:19,21:56))
subD8$condition = factor(subD8$condition)
lme.fit.bminusBB = lme(measure~condition, random=~1|ID, data=subD8)

# summary of model 
summary(lme.fit.bminusBB)

# get ANOVA results for lme model
anova.lme(lme.fit.bminusBB)

# run model
subD2 = subset(D_tall, ! condition %in% c(1:21,23:25,27:29,31:56))
lme.fit.aminpremidIS = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.aminpremidIS)

# get ANOVA results for lme model
anova.lme(lme.fit.aminpremidIS)



## BOOTSTRAPPED ANALYSES FOR THE B- EVENTS ##

# B- pre 
set.seed(2017)
dif.lm.fit.bminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==12,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminuspre = boot(D_tall, dif.lm.fit.bminuspre, R=5000) 
dif.lm.Bootobj.bminuspre
dif.lm.Bootobj.bminuspre.confint = 46.19048 + 1.96*c(-4.319935,4.319935)


# B- mid 
set.seed(2017)
dif.lm.fit.bminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==16,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminusmid = boot(D_tall, dif.lm.fit.bminusmid, R=5000) 
dif.lm.Bootobj.bminusmid
dif.lm.Bootobj.bminusmid.confint = 42.85714 + 1.96*c(-3.608798, 3.608798)

# B- post 
set.seed(2017)
dif.lm.fit.bminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==20,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.bminuspost = boot(D_tall, dif.lm.fit.bminuspost, R=5000) 
dif.lm.Bootobj.bminuspost
dif.lm.Bootobj.bminuspost.confint = 51.42857 + 1.96*c(-6.977627, 6.977627)



## PERMUTATION ANALYSIS FOR THE B+ EVENTS ##
# B- pre vs mid
# create a function that computes the difference between B- pre vs post
subBminuspremid = subset(D_tall, ! condition %in% c(1:11,13:15,17:56))
subBminuspremid$condition = factor(subBminuspremid$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminuspremid$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminuspremid) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminuspremid)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001



# B- pre vs post
# create a function that computes the difference between B- pre vs post
subBminusprepost = subset(D_tall, ! condition %in% c(1:11,13:19,21:56))
subBminusprepost$condition = factor(subBminusprepost$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusprepost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusprepost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusprepost)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001



# B- mid vs post
# create a function that computes the difference between B- pre vs post
subBminusmidpost = subset(D_tall, ! condition %in% c(1:15,17:19,21:56))
subBminusmidpost$condition = factor(subBminusmidpost$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusmidpost$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusmidpost) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusmidpost)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #0.0001





#######################################
# B+ post across BB and IS conditions #
#######################################

# B+ BB vs IS permutation test
# create a function that computes the difference between B- pre vs post
subBplusBBIS = subset(D_tall, ! condition %in% c(1:18,20:30,32:56))
subBplusBBIS$condition = factor(subBplusBBIS$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusBBIS$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusBBIS) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusBBIS)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000 #p-value, two-tailed #0.0001


# B- BB vs IS permutation test
# create a function that computes the difference between B- pre vs post
subBminusBBIS = subset(D_tall, ! condition %in% c(1:19,21:31,33:56))
subBminusBBIS$condition = factor(subBminusBBIS$condition)

# permutation test for B+ pre vs post ratings

b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusBBIS$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusBBIS) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusBBIS)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000 #p-value, two-tailed #0.0001





##################
## 1C CONDITION ##
##################
 
## A+ condition across the pre-, mid-, and post-rating phases

# run model
subD10 = subset(D_tall, ! condition %in% c(1:32,34:36,38:40,42:56))
subD10$condition = factor(subD10$condition)
lme.fit.aplus1C = lme(measure~condition, random=~1|ID, data=subD10)

# summary of model 
summary(lme.fit.aplus1C)

# get ANOVA results for lme model
anova.lme(lme.fit.aplus1C)

## BOOTSTRAPPED ANALYSES FOR THE A+ EVENTS ##

# A+ pre 
set.seed(2017)
dif.lm.fit.Apluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==33,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Apluspre = boot(D_tall, dif.lm.fit.Apluspre, R=5000) 
dif.lm.Bootobj.Apluspre
dif.lm.Bootobj.Apluspre.confint = 70 + 1.96*c(-4.909101,4.909101)

# A+ mid 
set.seed(2017)
dif.lm.fit.Aplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==37,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aplusmid = boot(D_tall, dif.lm.fit.Aplusmid, R=5000) 
dif.lm.Bootobj.Aplusmid
dif.lm.Bootobj.Aplusmid.confint = 87.14286 + 1.96*c(-4.06697,4.06697)

# A+ post 
set.seed(2017)
dif.lm.fit.Apluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==41,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Apluspost = boot(D_tall, dif.lm.fit.Apluspost, R=5000) 
dif.lm.Bootobj.Apluspost
dif.lm.Bootobj.Apluspost.confint = 90.2381 + 1.96*c(-3.372631,3.372631)

## PERMUTATION TESTS ##

# A+ pre vs mid 
subApluspremid1C = subset(D_tall, ! condition %in% c(1:32,34:36,38:56))
subApluspremid1C$condition = factor(subApluspremid1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subApluspremid1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subApluspremid1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000

# A+ mid vs post 
subAplusmidpost1C = subset(D_tall, ! condition %in% c(1:36,38:40,42:56))
subAplusmidpost1C$condition = factor(subAplusmidpost1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusmidpost1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusmidpost1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusmidpost1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000





## A- condition across the pre-, mid-, and post-rating phases

# run model
subD11 = subset(D_tall, ! condition %in% c(1:33,35:37,39:41,43:56))
subD11$condition = factor(subD11$condition)
lme.fit.aminus1C = lme(measure~condition, random=~1|ID, data=subD11)

# summary of model 
summary(lme.fit.aplus1C)

# get ANOVA results for lme model
anova.lme(lme.fit.aminus1C)

## BOOTSTRAPPED ANALYSES FOR THE A- EVENTS ##

# A- pre 
set.seed(2017)
dif.lm.fit.Aminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==34,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminuspre = boot(D_tall, dif.lm.fit.Aminuspre, R=5000) 
dif.lm.Bootobj.Aminuspre
dif.lm.Bootobj.Aminuspre.confint = 44.52381 + 1.96*c(-3.135468,3.135468)

# A- mid 
set.seed(2017)
dif.lm.fit.Aminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==38,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminusmid = boot(D_tall, dif.lm.fit.Aminusmid, R=5000) 
dif.lm.Bootobj.Aminusmid
dif.lm.Bootobj.Aminusmid.confint = 26.66667 + 1.96*c(-6.12378,6.12378)

# A- post 
set.seed(2017)
dif.lm.fit.Aminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==42,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminuspost = boot(D_tall, dif.lm.fit.Aminuspost, R=5000) 
dif.lm.Bootobj.Aminuspost
dif.lm.Bootobj.Aminuspost.confint = 22.61905 + 1.96*c(-5.778659,5.778659)
# 11.29288 33.94522

## PERMUTATION TESTS ##

# A- pre vs mid 
subAminuspremid1C = subset(D_tall, ! condition %in% c(1:33,35:37,39:56))
subAminuspremid1C$condition = factor(subAminuspremid1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminuspremid1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminuspremid1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminuspremid1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000

# A- mid vs post 
subAminusmidpost1C = subset(D_tall, ! condition %in% c(1:36,38:40,42:56))
subAminusmidpost1C$condition = factor(subAminusmidpost1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminusmidpost1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusmidpost1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusmidpost1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000





## B+ condition across the pre-, mid-, and post-rating phases

# run model
subD14 = subset(D_tall, ! condition %in% c(1:34,36:38,40:42,44:56))
subD14$condition = factor(subD14$condition)
lme.fit.bplus1C = lme(measure~condition, random=~1|ID, data=subD14)

# summary of model 
summary(lme.fit.bplus1C)

# get ANOVA results for lme model
anova.lme(lme.fit.bplus1C)

## BOOTSTRAPPED ANALYSES FOR THE B+ EVENTS ##

# B+ pre 
set.seed(2017)
dif.lm.fit.Bpluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==35,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bpluspre = boot(D_tall, dif.lm.fit.Bpluspre, R=5000) 
dif.lm.Bootobj.Bpluspre
dif.lm.Bootobj.Bpluspre.confint = 70 + 1.96*c(-4.897462,4.897462)
# 60.40097 79.59903

# B+ mid 
set.seed(2017)
dif.lm.fit.Bplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==39,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bplusmid = boot(D_tall, dif.lm.fit.Bplusmid, R=5000) 
dif.lm.Bootobj.Bplusmid
dif.lm.Bootobj.Bplusmid.confint = 25.04762 + 1.96*c(-5.257135,5.257135)
# 14.74364 35.35160

# B+ post 
set.seed(2017)
dif.lm.fit.Bpluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==43,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bpluspost = boot(D_tall, dif.lm.fit.Bpluspost, R=5000) 
dif.lm.Bootobj.Bpluspost
dif.lm.Bootobj.Bpluspost.confint = 23.2381 + 1.96*c(-5.356558,5.356558)
# 12.73925 33.73695

## PERMUTATION TESTS ##

# B+ pre vs mid 
subBpluspremid1C = subset(D_tall, ! condition %in% c(1:34,36:38,40:56))
subBpluspremid1C$condition = factor(subBpluspremid1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBpluspremid1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBpluspremid1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBpluspremid1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000

# B+ mid vs post 
subBplusmidpost1C = subset(D_tall, ! condition %in% c(1:38,40:42,44:56))
subBplusmidpost1C$condition = factor(subBplusmidpost1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusmidpost1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusmidpost1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusmidpost1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000



## B- condition across the pre-, mid-, and post-rating phases

# run model
subD15 = subset(D_tall, ! condition %in% c(1:35,37:39,41:43,45:56))
subD15$condition = factor(subD15$condition)
lme.fit.bminus1C = lme(measure~condition, random=~1|ID, data=subD15)

# summary of model 
summary(lme.fit.bminus1C)

# get ANOVA results for lme model
anova.lme(lme.fit.bminus1C)

## BOOTSTRAPPED ANALYSES FOR THE B- EVENTS ##

# B- pre 
set.seed(2017)
dif.lm.fit.Bminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==36,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminuspre = boot(D_tall, dif.lm.fit.Bminuspre, R=5000) 
dif.lm.Bootobj.Bminuspre
dif.lm.Bootobj.Bminuspre.confint = 44.52381 + 1.96*c(-3.136366,3.136366)
# 38.37653 50.67109

# B- mid 
set.seed(2017)
dif.lm.fit.Bminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==40,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminusmid = boot(D_tall, dif.lm.fit.Bminusmid, R=5000) 
dif.lm.Bootobj.Bminusmid
dif.lm.Bootobj.Bminusmid.confint = 73.42857 + 1.96*c(-6.565599,6.565599)
# 60.56000 86.29714

# B- post 
set.seed(2017)
dif.lm.fit.Bminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==44,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminuspost = boot(D_tall, dif.lm.fit.Bminuspost, R=5000) 
dif.lm.Bootobj.Bminuspost
dif.lm.Bootobj.Bminuspost.confint = 73.38095 + 1.96*c(-6.805077,6.805077)
# 60.0430 86.7189


## PERMUTATION TESTS ##

# B- pre vs mid 
subBpminuspremid1C = subset(D_tall, ! condition %in% c(1:35,37:39,41:56))
subBpminuspremid1C$condition = factor(subBpminuspremid1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBpminuspremid1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBpminuspremid1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBpminuspremid1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000

# B- mid vs post 
subBminusmidpost1C = subset(D_tall, ! condition %in% c(1:38,40:42,44:56))
subBminusmidpost1C$condition = factor(subBminusmidpost1C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusmidpost1C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusmidpost1C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusmidpost1C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000




##################
## 2C CONDITION ##
##################
## A+ condition across the pre-, mid-, and post-rating phases

# run model
subD18 = subset(D_tall, ! condition %in% c(1:44,46:48,50:52,54:56))
subD18$condition = factor(subD18$condition)
lme.fit.Aplus2c = lme(measure~condition, random=~1|ID, data=subD18)

# summary of model 
summary(lme.fit.Aplus2c)

# get ANOVA results for lme model
anova.lme(lme.fit.Aplus2c)

## BOOTSTRAPPED ANALYSES FOR THE A+ EVENTS ##

# A+ pre 
set.seed(2017)
dif.lm.fit.Apluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==45,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Apluspre = boot(D_tall, dif.lm.fit.Apluspre, R=5000) 
dif.lm.Bootobj.Apluspre
dif.lm.Bootobj.Apluspre.confint = 57.61905 + 1.96*c(-3.686809,3.686809)
# 50.3929 64.8452

# A+ mid 
set.seed(2017)
dif.lm.fit.Aplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==49,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aplusmid = boot(D_tall, dif.lm.fit.Aplusmid, R=5000) 
dif.lm.Bootobj.Aplusmid
dif.lm.Bootobj.Aplusmid.confint = 92.38095 + 1.96*c(-2.353957,2.353957)
# 87.76719 96.99471

# A+ post 
set.seed(2017)
dif.lm.fit.Apluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==53,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Apluspost = boot(D_tall, dif.lm.fit.Apluspost, R=5000) 
dif.lm.Bootobj.Apluspost
dif.lm.Bootobj.Apluspost.confint = 94.04762 + 1.96*c(-2.143268,2.143268)
# 89.84681 98.24843


## PERMUTATION TESTS ##

# A+ pre vs mid 
subApluspremid2C = subset(D_tall, ! condition %in% c(1:44,46:48,50:56))
subApluspremid2C$condition = factor(subApluspremid2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subApluspremid2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subApluspremid2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subApluspremid2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000


# A+ mid vs post 
subAplusmidpost2C = subset(D_tall, ! condition %in% c(1:48,50:52,54:56))
subAplusmidpost2C$condition = factor(subAplusmidpost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAplusmidpost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAplusmidpost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAplusmidpost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000






################################################################
## A- condition across the pre-, mid-, and post-rating phases ##
################################################################

# run model
subD19 = subset(D_tall, ! condition %in% c(1:45,47:49,51:53,55:56))
subD19$condition = factor(subD19$condition)
lme.fit.Aminus2c = lme(measure~condition, random=~1|ID, data=subD19)

# summary of model 
summary(lme.fit.Aminus2c)

# get ANOVA results for lme model
anova.lme(lme.fit.Aminus2c)


## BOOTSTRAPPED ANALYSES FOR THE A- EVENTS ##

# A- pre 
set.seed(2017)
dif.lm.fit.Aminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==46,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminuspre = boot(D_tall, dif.lm.fit.Aminuspre, R=5000) 
dif.lm.Bootobj.Aminuspre
dif.lm.Bootobj.Aminuspre.confint = 44.28571 + 1.96*c(-6.093333,6.093333)
# 32.34278 56.22864

# A- mid 
set.seed(2017)
dif.lm.fit.Aminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==50,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminusmid = boot(D_tall, dif.lm.fit.Aminusmid, R=5000) 
dif.lm.Bootobj.Aminusmid
dif.lm.Bootobj.Aminusmid.confint = 25.95238 + 1.96*c(-6.732869,6.732869)
# 12.75596 39.14880

# A- post 
set.seed(2017)
dif.lm.fit.Aminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==54,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Aminuspost = boot(D_tall, dif.lm.fit.Aminuspost, R=5000) 
dif.lm.Bootobj.Aminuspost
dif.lm.Bootobj.Aminuspost.confint = 16.90476 + 1.96*c(-6.070494,6.070494)
# 5.006592 28.802928


## PERMUTATION TESTS ##

# A- pre vs mid 
subAminuspremid2C = subset(D_tall, ! condition %in% c(1:45,47:49,51:56))
subAminuspremid2C$condition = factor(subAminuspremid2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminuspremid2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminuspremid2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminuspremid2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000


# A- mid vs post 
subAminusmidpost2C = subset(D_tall, ! condition %in% c(1:49,51:53,55:56))
subAminusmidpost2C$condition = factor(subAminusmidpost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subAminusmidpost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subAminusmidpost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subAminusmidpost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000


################################################################
## B+ condition across the pre-, mid-, and post-rating phases ##
################################################################

# run model
subD21 = subset(D_tall, ! condition %in% c(1:46,48:54,56))
subD21$condition = factor(subD21$condition)
lme.fit.Bplus2c = lme(measure~condition, random=~1|ID, data=subD21)

# summary of model 
summary(lme.fit.Bplus2c)

# get ANOVA results for lme model
anova.lme(lme.fit.Bplus2c)


## BOOTSTRAPPED ANALYSES FOR THE A- EVENTS ##

# B+ pre 
set.seed(2017)
dif.lm.fit.Bpluspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==47,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bpluspre = boot(D_tall, dif.lm.fit.Bpluspre, R=5000) 
dif.lm.Bootobj.Bpluspre
dif.lm.Bootobj.Bpluspre.confint = 57.61905 + 1.96*c(-3.642884,3.642884)
# 50.4790 64.7591

# B+ mid 
set.seed(2017)
dif.lm.fit.Bplusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==51,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bplusmid = boot(D_tall, dif.lm.fit.Bplusmid, R=5000) 
dif.lm.Bootobj.Bplusmid
dif.lm.Bootobj.Bplusmid.confint = 47.14286 + 1.96*c(-4.094394,4.094394)
# 39.11785 55.16787

# B+ post 
set.seed(2017)
dif.lm.fit.Bpluspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==55,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bpluspost = boot(D_tall, dif.lm.fit.Bpluspost, R=5000) 
dif.lm.Bootobj.Bpluspost
dif.lm.Bootobj.Bpluspost.confint = 65.90476 + 1.96*c(-2.304792,2.304792)
# 61.38737 70.42215


## PERMUTATION TESTS ##

# B+ pre vs mid 
subBpluspremid2C = subset(D_tall, ! condition %in% c(1:46,48:50,52:56))
subBpluspremid2C$condition = factor(subBpluspremid2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBpluspremid2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBpluspremid2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBpluspremid2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000
hist(b)
abline(v = beta_actual, col = "blue", lwd = 2)

# B+ mid vs post 
subBplusmidpost2C = subset(D_tall, ! condition %in% c(1:50,52:54,56))
subBplusmidpost2C$condition = factor(subBplusmidpost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusmidpost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusmidpost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusmidpost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000



# B+ pre vs post 
subBplusprepost2C = subset(D_tall, ! condition %in% c(1:46,48:54,56))
subBplusprepost2C$condition = factor(subBplusprepost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBplusprepost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBplusprepost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBplusprepost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000



################################################################
## B- condition across the pre-, mid-, and post-rating phases ##
################################################################

# run model
subD24 = subset(D_tall, ! condition %in% c(1:47,49:51,53:55))
subD24$condition = factor(subD24$condition)
lme.fit.Bminus2c = lme(measure~condition, random=~1|ID, data=subD24)

# summary of model 
summary(lme.fit.Bminus2c)

# get ANOVA results for lme model
anova.lme(lme.fit.Bminus2c)


## BOOTSTRAPPED ANALYSES FOR THE A- EVENTS ##

# B- pre 
set.seed(2017)
dif.lm.fit.Bminuspre = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==48,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminuspre = boot(D_tall, dif.lm.fit.Bminuspre, R=5000) 
dif.lm.Bootobj.Bminuspre
dif.lm.Bootobj.Bminuspre.confint = 44.28571 + 1.96*c(-6.177464,6.177464)
# 32.17788 56.39354

# B- mid 
set.seed(2017)
dif.lm.fit.Bminusmid = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==52,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminusmid = boot(D_tall, dif.lm.fit.Bminusmid, R=5000) 
dif.lm.Bootobj.Bminusmid
dif.lm.Bootobj.Bminusmid.confint = 45.85714 + 1.96*c(-5.460093,5.460093)
# 39.11785 55.16787

# B- post 
set.seed(2017)
dif.lm.fit.Bminuspost = function(data,b,formula){ 
  d= data[b,] 
  dif.1 =  mean(d[d$condition==56,4], data=D_tall) 
  return(dif.1)
}

dif.lm.Bootobj.Bminuspost = boot(D_tall, dif.lm.fit.Bminuspost, R=5000) 
dif.lm.Bootobj.Bminuspost
dif.lm.Bootobj.Bminuspost.confint = 41.71429 + 1.96*c(-3.839308,3.839308)
# 34.18925 49.23933


## PERMUTATION TESTS ##

# B- pre vs mid 
subBminuspremid2C = subset(D_tall, ! condition %in% c(1:47,49:51,53:56))
subBminuspremid2C$condition = factor(subBminuspremid2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminuspremid2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminuspremid2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminuspremid2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) > beta_actual)/5000


# B- mid vs post 
subBminusmidpost2C = subset(D_tall, ! condition %in% c(1:51,53:55))
subBminusmidpost2C$condition = factor(subBminusmidpost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusmidpost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusmidpost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusmidpost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000

# B- pre vs post 
subBminusprepost2C = subset(D_tall, ! condition %in% c(1:47,49:55))
subBminusprepost2C$condition = factor(subBminusprepost2C$condition)


b = rep(0,5000) 
for(i in 1:5000){
  y = sample(subBminusprepost2C$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subBminusprepost2C) 
  b[i] = fixed.effects(lm_1)[2]
}

lm.fit = lme(measure~condition, random=~1|ID, data=subBminusprepost2C)
beta_actual = fixed.effects(lm.fit)[2]

sum(abs(b) < beta_actual)/5000


#####
# Clustered Bar Graph
#####
F = D
F = as.data.frame(F)
F = F[,-c(2:10)]
fix(F)

# names: "ID"           "test.trial"   "measure"      "condition"    "period"       "test.trial.2" 

# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 2:49, v.names = "measure", timevar = "test.trial", idvar = "ID", new.row.names = 1:1008, direction = "long")
F_tall$condition = rep(c("BB","IS","1C","2C"), each = 252)
F_tall$condition = factor(F_tall$condition, levels = c("BB", "IS", "1C", "2C"))

F_tall$period = rep(c("pre","mid","post"), each = 84)
F_tall$period = factor(F_tall$period, levels = c("pre", "mid", "post")) # you did this, as opposed to "as.factor"
                                                                        # to ensure that the labels are ordered
                                                                        # in THAT particular way
F_tall$period = as.factor(F_tall$period)



F_tall$test.trial.2 = rep(c("A+", "A-", "B+", "B-"), each = 21, times = 4)
F_tall$test.trial.2 = as.factor(F_tall$test.trial.2)

F_tall = F_tall[, c("ID", "test.trial","condition", "test.trial.2","measure", "period")] # reoder column names


condition_barplot = ggplot(F_tall, aes(test.trial.2, measure, fill = period)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  facet_wrap(~condition) + # create as many separate graphs as there are conditions 
  scale_x_discrete(limits=c("A+","A-","B+", "B-")) + 
  ylab("causal (likelihood) ratings") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) + # ensure that bars hit the x-axis
  theme(strip.background =element_rect(fill='coral2')) +
  theme(strip.text = element_text(colour = 'white')) + # change facet labels color and text color
  labs(x = "Test trials") # change the main x-axis label
