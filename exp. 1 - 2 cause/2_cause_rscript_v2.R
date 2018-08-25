########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 1 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################

## INITIAL SET UP ##
# import "2_cause_CSV.csv"
D = read.csv(file.choose(), header = TRUE)

# convert data from "wide" to "tall" format
D_tall = reshape(D, varying = 3:18, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:960, direction = "long")

# order by ID
D_tall = D_tall[order(D_tall$ID),]

# remove scientific notation
options(scipen=999)

# add a binary column to the data frame
D_tall$measure.bin = rep(0, 960)
for(i in 1:960){
  bin.measure = ifelse(D_tall$measure[i]>=50,1,0)
  D_tall$measure.bin[i] = bin.measure
}

# install boostrap library
library(boot)


###########################################################################
#############   chi-square to determine if the choice (dependent) #########
############    is independent of the cause vs non.cause          #########     
###########################################################################

intro.chi = c(120,0)

intro.chi.test = chisq.test(intro.chi)

########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

##==> Normality check <==##

# plot the histograms
par(mfrow=c(4,4)) 
for (ii in 1:16)  hist(D_tall$measure[D_tall$condition==ii], breaks=5)
par(mfrow=c(1,1)) 

# formal test of normality
shapiro.ps = rep(0,16)
for(i in 1:16) {
  shap.calc = shapiro.test(D_tall$measure[D_tall$condition==i])
  shapiro.ps[i] = shap.calc$p.value
}

shapiro.mat = matrix(NA, nrow=4, ncol=4, byrow = TRUE)
for(ii in 1:16){
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
install.packages("car")
library(car)
leveneTest(D_tall$measure, as.factor(D_tall$condition), center=median) # used 'median' because it's a better measure of central tendency given the non-normality


##==> summary of equal variance check <==##
# Based on the analysis above, there is evidence of unequal variance.
# Violations were indicated by p-values that were less than .005.
# Corroborating the 'normality" analyses above, then, subsequent analyses
# will use bootstrapping and parametric tests.
###########################################


########################################################
########################################################
########################################################
#############                              #############
#############            Models            #############
#############                              #############
########################################################
########################################################
########################################################



########################################################
####          CONTROL CONDITITION ANALYSES          ####
########################################################

# install 'lme4' package to use the lmer function
install.packages("lme4")
library(lme4)
library(nlme)



          #######
  ######################
###########################
## A pre vs A post in 1C ##
###########################
  #######################
          #######
# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD = subset(D_tall, ! condition %in% c(1:8,10,12:16))
subD$condition = factor(subD$condition)

# run a lm and lmer model and compare the models using BIC

# partial model
lm.fit.1ca = lm(measure~factor(D_tall$condition, levels=c(9,11)), data=D_tall)
lmer.fit.1ca = lmer(measure~factor(D_tall$condition, levels=c(9,11))+(1|ID), data=D_tall)
lme.fit.1ca = lme(measure~condition, random=~1|ID, data=subD)

# summary of models
summary(lm.fit.1ca)
summary(lmer.fit.1ca)
summary(lme.fit.1ca)

# get ANOVA results for lme model
anova.lme(lme.fit.1ca)

# best fit models evaluations
BIC(lm.fit.1ca)
BIC(lmer.fit.1ca)


# obtained 95% Bootstrapped CI for pre- and post-ratings of A
bbB.lme.fit  <- function(d,i) {  # this worked for lme boot!
  d <- d[i,]
  # create new Subject id so each is unique
  #d$Subject <- 1:dim(d)[1]#
  thecoef <- tryCatch({bootfm.lme <- lme(measure ~ condition, data=d,
                                         random=~1|ID)
  fixef(bootfm.lme)},
  error = function(e) {print(e)
    return(rep(NA,length(fixef(lme.fit.1ca))))})
  #    browser()
  return(thecoef)
}

set.seed(2017)
bbB.lme.Bootobj = boot(subD, bbB.lme.fit , R=4000)
bbB.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Apre1c = 49.6  + 1.96*c(-2.027916, 2.027916)
ci.Apost1c = 45.21667 + 1.96*c(-2.970056,2.970056)

# derive a null distribution via permuation and obtain p-value for pre-
# and post-ratings of A
install.packages("pgirmess")
library(pgirmess)
PermTest(lme.fit.1ca, B=4000)

          #######
  ######################
###########################
## B pre vs B post in 1C ##
###########################
  #######################
          #######
# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD2 = subset(D_tall, ! condition %in% c(1:9,11,13:16))
subD2$condition = factor(subD2$condition)

# run model
lme.fit.1cb = lme(measure~condition, random=~1|ID, data=subD2)

# summary of model 
summary(lme.fit.1cb)

# get ANOVA results for lme model
anova.lme(lme.fit.1cb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD2, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Bpre1c = 49.56667  + 1.96*c(-2.198787, 2.198787)
ci.Bpost1c = -39.06667 + 1.96*c(-3.635452,3.635452)


### One-factor Repeated Measures ANOVA ###
install.packages("ez")
library(ez)

# subset data to include pre- and post-ratings of A and B in the 1C condition 
subD3 = subset(D_tall, ! condition %in% c(1:8,13:16))
subD3$condition = factor(subD3$condition)

# run repeated-measures ANOVA
repeat_AOV = ezANOVA(subD3, dv = measure, wid = ID, within = as.factor(condition))
print(repeat_AOV) # note that 'print' here is like 'summary' for lm, lmer, glm, glmer models

# test contrasts of Apre vs Apost and Bpre vs Bpost w/ Bonferroni correction
# but first determine whether the difference scores are normally distributed
# create new data frame 'F' which is just a copy of the 'D' data frame and then
# add a column that computes the difference between A pre and post in 1c and then another one
# for B pre and post for 1C

F = D
# adding a pre and post column for 1c
F$diffApreApost1c = rep(0, 60)
for(i in 1:length(F$diffApreApost1c)){
  calc = (F$oc.a.pre[i]-F$oc.a.post)
  F$diffApreApost1c[i] = calc
}

# adding b pre and post column for 1c
F$diffBpreBpost1c = rep(0, 60)
for(i in 1:length(F$diffBpreABpost1c)){
  calc = (F$oc.b.pre[i]-F$oc.b.post)
  F$diffBpreBpost1c[i] = calc
}

# transform F from wide to tall format
# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 3:20, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:1080, direction = "long")

# order by ID
F_tall = F_tall[order(F_tall$ID),]

# test to see if difference scores are normally distributed
apreapost_1c_norm_check = shapiro.test(F_tall$measure[F_tall$condition==17]) # clearly not normally distributed
bprebpost_1c_norm_check = shapiro.test(F_tall$measure[F_tall$condition==18]) # clearly not normally distributed

# because the dif scores above were not normally distributed, use wilcoxon sign rank test for paired t-test with
# Bonferroni correction 
apre_apost_1c_contrasts = wilcox.test(F_tall$measure[F_tall$condition==9],F_tall$measure[F_tall$condition==11],paired=TRUE)
bpre_bpost_1c_contrasts = wilcox.test(F_tall$measure[F_tall$condition==10],F_tall$measure[F_tall$condition==12],paired=TRUE)

# Bonferroni adjusted-alpha: .05/4 = .0125

F = Dr
F$diffApreApost1c = rep(0, 60)
for(i in 1:length(F$diffApreApost1c)){
  calc = (F$oc.a.pre[i]-F$oc.a.post)
  F$diffApreApost1c[i] = calc
}


          #######
  ######################
###########################
## A pre vs A post in 2C ##
###########################
  #######################
         #######
# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD3 = subset(D_tall, ! condition %in% c(1:12,14,16))
subD3$condition = factor(subD3$condition)

# run model
lme.fit.2ca = lme(measure~condition, random=~1|ID, data=subD3)

# summary of model 
summary(lme.fit.2ca)

# get ANOVA results for lme model
anova.lme(lme.fit.2ca)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
bbA2C.lme.fit  <- function(d,i) {  # this worked for lme boot!
  d <- d[i,]
  # create new Subject id so each is unique
  #d$Subject <- 1:dim(d)[1]#
  thecoef <- tryCatch({bootfm.lme <- lme(measure ~ condition, data=d,
                                         random=~1|ID)
  fixef(bootfm.lme)},
  error = function(e) {print(e)
    return(rep(NA,length(fixef(lme.fit.1ca))))})
  #    browser()
  return(thecoef)
}

set.seed(2017)
bbA2c.lme.Bootobj = boot(subD3, bbA2C.lme.fit , R=4000)
bbA2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Apre2c = 49.93333  + 1.96*c(-2.031376, 2.031376)
ci.Apost2c = 44.81667 + 1.96*c(-2.769960,2.769960)


         #######
  ######################
###########################
## B pre vs B post in 2C ##
###########################
  #######################
          #######
# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD4 = subset(D_tall, ! condition %in% c(1:13,15))
subD4$condition = factor(subD4$condition)

# run model
lme.fit.2cb = lme(measure~condition, random=~1|ID, data=subD4)

# summary of model 
summary(lme.fit.2cb)

# get ANOVA results for lme model
anova.lme(lme.fit.2cb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
bbB2C.lme.fit  <- function(d,i) {  # this worked for lme boot!
  d <- d[i,]
  # create new Subject id so each is unique
  #d$Subject <- 1:dim(d)[1]#
  thecoef <- tryCatch({bootfm.lme <- lme(measure ~ condition, data=d,
                                         random=~1|ID)
  fixef(bootfm.lme)},
  error = function(e) {print(e)
    return(rep(NA,length(fixef(lme.fit.1ca))))})
  #    browser()
  return(thecoef)
}

set.seed(2017)
bbB2c.lme.Bootobj = boot(subD4, bbB2C.lme.fit , R=4000)
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Bpre2c = 48.75  + 1.96*c(-1.982465, 1.982465)
ci.Bpost2c = 29.60 + 1.96*c(-3.102846,3.102846)


# one sample t-test to determine if post-rating differs
# from the rating that would be expected if participants
# were tracking the statistics
postB2c.ttest = t.test (D_tall$measure[D_tall$condition==16], mu=66.67)


### One-factor Repeated Measures ANOVA to test pre and post ratings of A and B in 2C ###
install.packages("ez")
library(ez)

# subset data to include pre- and post-ratings of A and B in the 1C condition 
subD6 = subset(D_tall, ! condition %in% c(1:12))
subD6$condition = factor(subD6$condition)

# because sphericity wasn't met, correct for degrees of freedom
3/.7
177/.7

# run repeated-measures ANOVA
repeat_AOV = ezANOVA(subD6, dv = measure, wid = ID, within = as.factor(condition))
print(repeat_AOV) # note that 'print' here is like 'summary' for lm, lmer, glm, glmer models


# adding a pre and post column for 2c
F$diffApreApost2c = rep(0, 60)
for(i in 1:length(F$diffApreApost2c)){
  calc = (F$tc.a.pre[i]-F$tc.a.post)
  F$diffApreApost2c[i] = calc
}

# adding b pre and post column for 2c
F$diffBpreBpost2c = rep(0, 60)
for(i in 1:length(F$diffBpreBpost2c)){
  calc = (F$tc.b.pre[i]-F$tc.b.post)
  F$diffBpreBpost2c[i] = calc
}


# transform F from wide to tall format
# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 3:22, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:1200, direction = "long")

# order by ID
F_tall = F_tall[order(F_tall$ID),]


# test to see if difference scores are normally distributed
apreapost_2c_norm_check = shapiro.test(F_tall$measure[F_tall$condition==19]) # clearly not normally distributed
bprebpost_2c_norm_check = shapiro.test(F_tall$measure[F_tall$condition==20]) # clearly not normally distributed

# because the dif scores above were not normally distributed, use wilcoxon sign rank test for paired t-test with
# Bonferroni correction 
apre_apost_2c_contrasts = wilcox.test(F_tall$measure[F_tall$condition==13],F_tall$measure[F_tall$condition==15],paired=TRUE)
bpre_bpost_2c_contrasts = wilcox.test(F_tall$measure[F_tall$condition==14],F_tall$measure[F_tall$condition==16],paired=TRUE)

# Bonferroni adjusted-alpha: .05/4 = .0125

# loop to obtain means of pre and post ratings of A and B in 2C
easy = rep(0,4)
for(i in 13:16){
  calc = mean(D_tall$measure[D_tall$condition==i])
  easy[i] = calc
}



########################################################
####          MAIN CONDITITION ANALYSES          ####
########################################################


          #######
  ######################
###########################
## A pre vs A post in IS ##
###########################
  #######################
          #######

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD7 = subset(D_tall, ! condition %in% c(1:4,6,8:16))
subD7$condition = factor(subD7$condition)

# run model
lme.fit.ISa = lme(measure~condition, random=~1|ID, data=subD7)


# summary of model 
summary(lme.fit.ISa)

# get ANOVA results for lme model
anova.lme(lme.fit.ISa)

set.seed(2017)
bbB2c.lme.Bootobj = boot(subD7, bbB2C.lme.fit , R=4000) # note that I didn't rename the object
                                                        # out of laziness, but clearly
                                                        # I'm using the pre- and post-ratings
                                                        # of A in the IS condition.
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreIS = 55.20000 + 1.96*c(-2.574410, 2.574410)
ci.ApostIS = -43.16667 + 1.96*c(-4.414155,4.414155)


          #######
  ######################
###########################
## B pre vs B post in IS ##
###########################
  #######################
          #######


# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD8 = subset(D_tall, ! condition %in% c(1:5,7,9:16))
subD8$condition = factor(subD8$condition)

# run model
lme.fit.ISb = lme(measure~condition, random=~1|ID, data=subD8)


# summary of model 
summary(lme.fit.ISb)

# get ANOVA results for lme model
anova.lme(lme.fit.ISb)

set.seed(2017)
bbB2c.lme.Bootobj = boot(subD8, bbB2C.lme.fit , R=4000) # note that I didn't rename the object
# out of laziness, but clearly
# I'm using the pre- and post-ratings
# of A in the IS condition.
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BpreIS = 49.50000 + 1.96*c(-2.085494, 2.085494)
ci.BpostIS = 40.81667 + 1.96*c(-3.111538,3.111538)

# Bayes' Factor analysis

# define the null and alternative models #
BIS.lme.null = lme(measure~1, random=~1|ID, data=subD8)
BIS.lme.alt = lme(measure~condition, random=~1|ID, data=subD8)



        #######
  ######################
###########################
## A pre vs A post in BB ##
###########################
  #######################
        #######

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD9 = subset(D_tall, ! condition %in% c(2,4:16))
subD9$condition = factor(subD9$condition)

# run model
lme.fit.BBa = lme(measure~condition, random=~1|ID, data=subD9)


# summary of model 
summary(lme.fit.BBa)

# get ANOVA results for lme model
anova.lme(lme.fit.BBa)

set.seed(2017)
bbB2c.lme.Bootobj = boot(subD9, bbB2C.lme.fit , R=4000) # note that I didn't rename the object
# out of laziness, but clearly
# I'm using the pre- and post-ratings
# of B in the BB condition.
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreBB = 51.25 + 1.96*c(-1.641334, 1.641334)
ci.ApostBB = 45.00 + 1.96*c(-2.344856,2.344856)


        #######
######################
###########################
## B pre vs B post in BB ##
###########################
#######################
        #######

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD10 = subset(D_tall, ! condition %in% c(1,3,5:16))
subD10$condition = factor(subD10$condition)

# run model
lme.fit.BBb = lme(measure~condition, random=~1|ID, data=subD10)


# summary of model 
summary(lme.fit.BBb)

# get ANOVA results for lme model
anova.lme(lme.fit.BBb)

# PERMUTATION TESTING
b = rep(0,4000) 

for(i in 1:4000){
  y = sample(subD10$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subD10) 
  b[i] = fixed.effects(lm_1)[2]
}

# run the model to see what the second coefficient actually is (or 1st if you're interested in the mean of reference group)
lm.fit = lme(measure~condition, random=~1|ID, data=subD10)
beta_actual = fixed.effects(lm.fit)[2]

# construct a histogram of using the 1000 values of b: 
hist(b)
abline(v = beta_actual, col = "blue", lwd = 2)

# place the beta_actual (the true beta) on the histogram to see how plausible it is under the H0:
segments(beta_actual,0,beta_actual,200, col="green")

# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b < beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/4000 #p-value, two-tailed
sum(b < beta_actual)/4000 #p-value, one-tailed


set.seed(2019)
bbB2c.lme.Bootobj = boot(subD10, bbB2C.lme.fit , R=4000) # note that I didn't rename the object
# out of laziness, but clearly
# I'm using the pre- and post-ratings
# of B in the BB condition.
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BpreBB = 46.50 + 1.96*c(-1.551589, 1.551589)
ci.BpostBB = -4.25 + 1.96*c(-2.596219,2.596219)




###############################################
###### Bayes Factor Analysis for B in BB ######
###############################################

# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=subD10)
lm.alt = lme(measure~condition, random=~1|ID, data=subD10)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01 
BF01 = exp((alt.bic - null.bic)/2) 
BF10 = 1/BF01

# transform BF01 into posterior probabilities - the percentage (i.e., posterior likelihood) of evidence in favor the null
BF01.posterior = BF01/(1+BF01)

# transform BF10 into posterior probabilities - the percentage (i.e., posterior likelihood) of evidence in favor the alternative (e.g., 4% likely under the null than the alternative)
BF01.posterior = BF10/(1+BF10)


        #######
  ######################
###########################
## B post BB vs B post in IS ##
###########################
  #######################
        #######


# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD11 = subset(D_tall, ! condition %in% c(1:3,5:7,9:16))
subD11$condition = factor(subD11$condition)

# run model
lme.fit.BBb = lme(measure~condition, random=~1|ID, data=subD11)


# summary of model 
summary(lme.fit.BBb)

# get ANOVA results for lme model
anova.lme(lme.fit.BBb)

# PERMUTATION TESTING
b = rep(0,4000) 

for(i in 1:4000){
  y = sample(subD11$measure, replace=TRUE)
  lm_1 = lme(y ~ condition, random=~1|ID, data=subD11) 
  b[i] = fixed.effects(lm_1)[2]
}

# run the model to see what the second coefficient actually is (or 1st if you're interested in the mean of reference group)
lm.fit = lme(measure~condition, random=~1|ID, data=subD10)
beta_actual = fixed.effects(lm.fit)[2]

# construct a histogram of using the 1000 values of b: 
hist(b)
abline(v = beta_actual, col = "blue", lwd = 2)

# place the beta_actual (the true beta) on the histogram to see how plausible it is under the H0:
segments(beta_actual,0,beta_actual,200, col="green")

# compute the p-value for your actual statistic based on the null distribution that was generated:
sum(b < beta_actual) #number of actual cases that are greater than the actual beta
sum(abs(b) < beta_actual)/4000 #p-value, two-tailed
sum(b < beta_actual)/4000 #p-value, one-tailed


set.seed(2017)
bbB2c.lme.Bootobj = boot(subD11, bbB2C.lme.fit , R=4000) # note that I didn't rename the object
# out of laziness, but clearly
# I'm using the pre- and post-ratings
# of B in the BB condition.
bbB2c.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BppostBB = 42.25000 + 1.96*c(-2.274312, 2.274312)
ci.BpostISvBB = 48.06667 + 1.96*c(-3.287434,3.287434)



#####
# Clustered Bar Graph
#####
F = D
F = as.data.frame(F)
fix(F)

# names: "ID"           "test.trial"   "measure"      "condition"    "period"       "test.trial.2" 

# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 3:18, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:960, direction = "long")
F_tall$condition = rep(c("BB","IS","1C","2C"), each = 240)
F_tall$condition = factor(F_tall$condition, levels = c("BB", "IS", "1C", "2C")) 

F_tall$period = rep(c("pre","post"), each = 120) # just take 240 and divide it by number of items in this vector, which is 2
F_tall$period = factor(F_tall$period, levels = c("pre", "post")) # you did this, as opposed to "as.factor"
# to ensure that the labels are ordered
# in THAT particular way


F_tall$test.trial.2 = rep(c("A", "B"), each = 60, times = 4)
F_tall$test.trial.2 = as.factor(F_tall$test.trial.2)

F_tall = F_tall[, c("ID", "condition", "test.trial.2","measure", "period")] # reoder column names


condition_barplot = ggplot(F_tall, aes(test.trial.2, measure, fill = period)) # create the bar graph with test.trial.2 on the x-axis and measure on the y-axis
condition_barplot + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + # add the bars, which represent the means and the place them side-by-side with 'dodge'
  stat_summary(fun.data=mean_cl_boot, geom = "errorbar", position = position_dodge(width=0.90), width = 0.2) + # add errors bars
  facet_wrap(~condition) + # create as many separate graphs as there are conditions 
  ylab("causal (likelihood) ratings") + # change the label of the y-axis
  theme_bw() + # remove the gray background
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + # remove the major and minor grids
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) + # ensure that bars hit the x-axis
  theme(strip.background =element_rect(fill='coral2')) +
  theme(strip.text = element_text(colour = 'white')) + # change facet labels color and text color
  labs(x = "Test trials") + # change the main x-axis label
  facet_wrap(~condition, scales = "free_x")
