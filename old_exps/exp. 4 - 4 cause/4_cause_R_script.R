########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 3 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################
## INITIAL SET UP ##
# import "4_cause_CSV.csv"
D = read.csv(file.choose(), header = TRUE)

D = D[-c(44:1013),]

# convert data from "wide" to "tall" format
D_tall = reshape(D, varying = 2:33, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:1920, direction = "long")

# order by ID
D_tall = D_tall[order(D_tall$ID),]

# remove scientific notation
options(scipen=999)


# libraries:
library(boot)
library(car)
library(lme4)
library(nlme)
library(pgirmess)
library(ez)

########################################################
#############                              #############
#############      Assumption Checks       #############
#############                              #############
########################################################

##==> Normality check <==##

# plot the histograms
par(mfrow=c(12,2)) 
for (ii in 1:16)  hist(D_tall$measure[D_tall$condition==ii], breaks=5)
par(mfrow=c(1,1)) 

# formal test of normality
shapiro.ps = rep(0,32)
for(i in 1:32) {
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

  ######################
###########################
## A pre vs A post in 1C ##
###########################
  #######################
# 17, 21
# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD = subset(D_tall, ! condition %in% c(1:16,18:20,22:32))
subD$condition = as.factor(subD$condition)
subD$condition = factor(subD$condition)

# run a lm and lmer model and compare the models using BIC

# partial model
lme.fit.1ca = lme(measure~condition, random=~1|ID, data=subD)

# summary of models
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
ci.Apre1c = 50.31667  + 1.96*c(-1.703637, 1.703637)
ci.Apost1c = 47.16667 + 1.96*c(-2.217616,2.217616)
ci.Apre1c
ci.Apost1c

  ######################
###########################
## B pre vs B post in 1C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD2 = subset(D_tall, ! condition %in% c(1:17,19:21,23:32))
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
ci.Bpre1c = 53.76667  + 1.96*c(-2.170496, 2.170496)
ci.Bpost1c = -45.45000 + 1.96*c(-3.257758,3.257758)
ci.Bpre1c
ci.Bpost1c

  ######################
###########################
## C pre vs C post in 1C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD3 = subset(D_tall, ! condition %in% c(1:18,20:22,24:32))
subD3$condition = factor(subD3$condition)

# run model
lme.fit.1cc = lme(measure~condition, random=~1|ID, data=subD3)

# summary of model 
summary(lme.fit.1cc)

# get ANOVA results for lme model
anova.lme(lme.fit.1cc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD3, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Cpre1c = 48.0500000  + 1.96*c(-2.048609, 2.048609)
ci.Cpost1c = -0.8833333 + 1.96*c(-3.139232,3.139232)
ci.Cpre1c
ci.Cpost1c

  ######################
###########################
## D pre vs D post in 1C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD4 = subset(D_tall, ! condition %in% c(1:19,21:23, 25:32))
subD4$condition = factor(subD4$condition)

# run model
lme.fit.1cd = lme(measure~condition, random=~1|ID, data=subD4)

# summary of model 
summary(lme.fit.1cd)

# get ANOVA results for lme model
anova.lme(lme.fit.1cd)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD4, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Dpre1c = 44.93333  + 1.96*c(-2.010829, 2.010829)
ci.Dpost1c = -5.10000 + 1.96*c(-3.005910,3.005910)
ci.Dpre1c
ci.Dpost1c


  ######################
###########################
## A pre vs A post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD5 = subset(D_tall, ! condition %in% c(1:24,26:28,30:32))
subD5$condition = factor(subD5$condition)

# run model
lme.fit.2ca = lme(measure~condition, random=~1|ID, data=subD5)

# summary of model 
summary(lme.fit.2ca)

# get ANOVA results for lme model
anova.lme(lme.fit.2ca)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD5, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Apre2c = 48.91667  + 1.96*c(-2.516262, 2.516262)
ci.Apost2c = 46.13333 + 1.96*c(-3.387557,3.387557)
ci.Apre2c 
ci.Apost2c
  ######################
###########################
## B pre vs B post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD6 = subset(D_tall, ! condition %in% c(1:25,27:29,31:32))
subD6$condition = factor(subD6$condition)

# run model
lme.fit.2cb = lme(measure~condition, random=~1|ID, data=subD6)

# summary of model 
summary(lme.fit.2cb)

# get ANOVA results for lme model
anova.lme(lme.fit.2cb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD6, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Apre2c = 44.37209  + 1.96*c(-2.175629, 2.175629)
ci.Apost2c = 24.17217 + 1.96*c(-3.458214,3.458214)
ci.Apre2c
ci.Apost2c
mean(D_tall$measure[D_tall$condition==26])
mean(D_tall$measure[D_tall$condition==30])

  ######################
###########################
## C pre vs C post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD7 = subset(D_tall, ! condition %in% c(1:26,28:30,32))
subD7$condition = factor(subD7$condition)

# run model
lme.fit.2cc = lme(measure~condition, random=~1|ID, data=subD7)

# summary of model 
summary(lme.fit.2cc)

# get ANOVA results for lme model
anova.lme(lme.fit.2cc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD7, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Cpre2c = 45.33333  + 1.96*c(-1.937039, 2.314974)
ci.Cpost2c = -2.86667 + 1.96*c(-2.693079,2.693079)
ci.Cpre2c
ci.Cpost2c 
mean(D_tall$measure[D_tall$condition==27])
mean(D_tall$measure[D_tall$condition==31])

  ######################
###########################
## D pre vs D post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD8 = subset(D_tall, ! condition %in% c(1:27,29:31))
subD8$condition = factor(subD8$condition)

# run model
lme.fit.2cd = lme(measure~condition, random=~1|ID, data=subD8)

# summary of model 
summary(lme.fit.2cd)

# get ANOVA results for lme model
anova.lme(lme.fit.2cd)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD8, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Dpre2c = 48.80000  + 1.96*c(-2.044887, 2.044887)
ci.Dpost2c = -2.13333 + 1.96*c(-3.199209,3.199209)
ci.Dpre2c
ci.Dpost2c
mean(D_tall$measure[D_tall$condition==28])
mean(D_tall$measure[D_tall$condition==32])

########################################################
####          MAIN CONDITITION ANALYSES             ####
########################################################
##################
## IS CONDITION ##
##################


  ######################
###########################
## A pre vs A post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD9 = subset(D_tall, ! condition %in% c(1:8,10:12,14:32))
subD9$condition = factor(subD9$condition)

# run model
lme.fit.ISa = lme(measure~condition, random=~1|ID, data=subD9)

# summary of model 
summary(lme.fit.ISa)

# get ANOVA results for lme model
anova.lme(lme.fit.ISa)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD9, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreIS = 50.33333   + 1.96*c(-2.2295, 2.2295)
ci.ApostIS = -39.71667 + 1.96*c(-3.4704,3.4704)
ci.ApreIS
ci.ApostIS
mean(D_tall$measure[D_tall$condition==9])
mean(D_tall$measure[D_tall$condition==13])
  ######################
###########################
## B pre vs B post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD10 = subset(D_tall, ! condition %in% c(1:9,11:13,15:32))
subD10$condition = factor(subD10$condition)

# run model
lme.fit.ISb = lme(measure~condition, random=~1|ID, data=subD10)

# summary of model 
summary(lme.fit.ISb)

# get ANOVA results for lme model
anova.lme(lme.fit.ISb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD10, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BpreIS = 48  + 1.96*c(-1.854653, 1.854653)
ci.BpostIS = 45.21667 + 1.96*c(-2.680070,2.680070)
ci.BpreIS
ci.BpostIS
mean(D_tall$measure[D_tall$condition==10])
mean(D_tall$measure[D_tall$condition==14])

  ######################
###########################
## C pre vs C post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD11 = subset(D_tall, ! condition %in% c(1:10,12:14,16:32))
subD11$condition = factor(subD11$condition)

# run model
lme.fit.ISc = lme(measure~condition, random=~1|ID, data=subD11)

# summary of model 
summary(lme.fit.ISc)

# get ANOVA results for lme model
anova.lme(lme.fit.ISc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD11, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpreIS = 50.250000   + 1.96*c(-1.929560, 1.929560)
ci.CpostIS = -5.516667 + 1.96*c(-2.876336,2.876336)
ci.CpreIS 
ci.CpostIS
  ######################
###########################
## D pre vs D post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD12 = subset(D_tall, ! condition %in% c(1:11,13:15,17:32))
subD12$condition = factor(subD12$condition)

# run model
lme.fit.ISd = lme(measure~condition, random=~1|ID, data=subD12)

# summary of model 
summary(lme.fit.ISd)

# get ANOVA results for lme model
anova.lme(lme.fit.ISd)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD12, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.DpreIS = 46   + 1.96*c(-1.436350, 1.436350)
ci.DpostIS = -0.6 + 1.96*c(-2.184581,2.184581)
ci.DpreIS
ci.DpostIS

##################
## BB CONDITION ##
##################

  ######################
###########################
## A pre vs A post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD13 = subset(D_tall, ! condition %in% c(2:4,6:32))
subD13$condition = factor(subD13$condition)

# run model
lme.fit.BBa = lme(measure~condition, random=~1|ID, data=subD13)

# summary of model 
summary(lme.fit.BBa)

# get ANOVA results for lme model
anova.lme(lme.fit.BBa)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD13, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreBB = 46.08333   + 1.96*c(-2.022054, 2.022054)
ci.ApostBB = 48.65000 + 1.96*c(-2.826058,2.826058)
ci.ApreBB
ci.ApostBB
  ######################
###########################
## B pre vs B post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD14 = subset(D_tall, ! condition %in% c(1,3:5,7:32))
subD14$condition = factor(subD14$condition)

# run model
lme.fit.BBb = lme(measure~condition, random=~1|ID, data=subD14)

# summary of model 
summary(lme.fit.BBb)

# get ANOVA results for lme model
anova.lme(lme.fit.BBb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD14, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BpreBB = 47.66667   + 1.96*c(-2.567940, 2.567940)
ci.BpostBB = -0.31667 + 1.96*c(-3.929385,3.929385)
ci.BpreBB
ci.BpostBB
  ######################
###########################
## C pre vs C post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD15 = subset(D_tall, ! condition %in% c(1:2,4:6,8:32))
subD15$condition = factor(subD15$condition)

# run model
lme.fit.BBc = lme(measure~condition, random=~1|ID, data=subD15)

# summary of model 
summary(lme.fit.BBc )

# get ANOVA results for lme model
anova.lme(lme.fit.BBc )

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD15, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpreBB = 44.26667   + 1.96*c(-1.969108, 1.969108)
ci.CpostBB = 0.50000 + 1.96*c(-2.997505,2.997505)
ci.CpreBB 
ci.CpostBB
  ######################
###########################
## D pre vs D post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD16 = subset(D_tall, ! condition %in% c(1:3,5:7,9:32))
subD16$condition = factor(subD16$condition)

# run model
lme.fit.BBd = lme(measure~condition, random=~1|ID, data=subD16)

# summary of model 
summary(lme.fit.BBd)

# get ANOVA results for lme model
anova.lme(lme.fit.BBd)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD16, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.DpreBB = 50.833333   + 1.96*c(-1.997886, 1.997886)
ci.DpostBB = -7.083333 + 1.96*c(-2.960815,2.960815)
ci.DpreBB 
ci.DpostBB
###############################################
###### Bayes Factor Analysis for D in BB ######
###############################################

# define the null and alternative models #
lm.null = lme(measure~1, random=~1|ID, data=subD16)
lm.alt = lme(measure~condition, random=~1|ID, data=subD16)

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


  ######################
###########################
## Bpost IS vs Bpost BB ##
###########################
#######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD17 = subset(D_tall, ! condition %in% c(1:5,7:13,15:32))
subD17$condition = factor(subD17$condition)

# run model
lme.fit.ISBBb = lme(measure~condition, random=~1|ID, data=subD17)

# summary of model 
summary(lme.fit.ISBBb)

# get ANOVA results for lme model
anova.lme(lme.fit.ISBBb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD17, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.DpreBB = 49.44186   + 1.96*c(-3.051812, 3.051812)
ci.DpostBB = 44.11628 + 1.96*c(-4.188242,4.188242)



#####
# Clustered Bar Graph
#####
F = D
F = F[-c(44:1013),]
F = as.data.frame(F)
fix(F)

# names: "ID"           "test.trial"   "measure"      "condition"    "period"       "test.trial.2" 

# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 2:33, v.names = "measure", timevar = "condition", new.row.names = 1:1920, idvar = "ID",  direction = "long")
F_tall$condition = rep(c("BB","IS","1C","2C"), each = 480)
F_tall$condition = factor(F_tall$condition, levels = c("BB", "IS", "1C", "2C")) 

F_tall$period = rep(c("pre","post"), each = 240) # just take 360 and divide it by number of items in this vector, which is 2
F_tall$period = factor(F_tall$period, levels = c("pre", "post")) # you did this, as opposed to "as.factor"
# to ensure that the labels are ordered
# in THAT particular way


F_tall$test.trial.2 = rep(c("A", "B", "C", "D"), each = 60, times = 4)
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




