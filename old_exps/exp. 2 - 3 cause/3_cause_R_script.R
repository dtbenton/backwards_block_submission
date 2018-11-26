########################################################
########################################################
########################################################
#############                              #############
#############      EXPERIMENT 2 SCRIPT     #############
#############                              #############
########################################################
########################################################
########################################################

# NOTE: Object 'C' serves the same role that B did in Experiment 1.

## INITIAL SET UP ##
# import "2_cause_CSV.csv"
D = read.csv(file.choose(), header = TRUE)

# convert data from "wide" to "tall" format
D_tall = reshape(D, varying = 3:26, v.names = "measure", timevar = "condition", idvar = "ID", new.row.names = 1:1440, direction = "long")

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
shapiro.ps = rep(0,24)
for(i in 1:24) {
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

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD = subset(D_tall, ! condition %in% c(1:12,14:15,17:24))
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
ci.Apre1c = 50.34883  + 1.96*c(-2.112401, 2.112401)
ci.Apost1c = 46.41117 + 1.96*c(-2.799613,2.799613)

  ######################
###########################
## B pre vs B post in 1C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD2 = subset(D_tall, ! condition %in% c(1:13,15:16,18:24))
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
ci.Bpre1c = 49.5655  + 1.96*c(-1.605429, 1.605429)
ci.Bpost1c = -3.2105 + 1.96*c(-2.418517,2.418517)


  ######################
###########################
## C pre vs C post in 1C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD3 = subset(D_tall, ! condition %in% c(1:14,16:17,19:24))
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
ci.Cpre1c = 50.04883  + 1.96*c(-2.225426, 2.225426)
ci.Cpost1c = -36.79883 + 1.96*c(-3.926308,3.926308)

  
  ######################
###########################
## A pre vs A post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD4 = subset(D_tall, ! condition %in% c(1:18,20:21,23:24))
subD4$condition = factor(subD4$condition)

# run model
lme.fit.2ca = lme(measure~condition, random=~1|ID, data=subD4)

# summary of model 
summary(lme.fit.2ca)

# get ANOVA results for lme model
anova.lme(lme.fit.2ca)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD4, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Apre2c = 49.11000   + 1.96*c(-2.091122, 2.091122)
ci.Apost2c = 47.01167 + 1.96*c(-2.830684,2.830684)

  ######################
###########################
## C pre vs C post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD5 = subset(D_tall, ! condition %in% c(1:20,22:23))
subD5$condition = factor(subD5$condition)

# run model
lme.fit.2cc = lme(measure~condition, random=~1|ID, data=subD5)

# summary of model 
summary(lme.fit.2cc)

# get ANOVA results for lme model
anova.lme(lme.fit.2cc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD5, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.Cpre2c = 49.02667   + 1.96*c(-2.281470, 2.281470)
ci.Cpost2c = 25.50667 + 1.96*c(-3.520436,3.520436)

  ######################
###########################
## B pre vs B post in 2C ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD6 = subset(D_tall, ! condition %in% c(1:19,21:22,24))
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
ci.Cpre2c = 46.710000   + 1.96*c(-1.684113, 1.684113)
ci.Cpost2c = -5.388333 + 1.96*c(-2.611612,2.611612)



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
subD7 = subset(D_tall, ! condition %in% c(1:6,8:9, 11:24))
subD7$condition = factor(subD7$condition)

# run model
lme.fit.ISa = lme(measure~condition, random=~1|ID, data=subD7)

# summary of model 
summary(lme.fit.ISa)

# get ANOVA results for lme model
anova.lme(lme.fit.ISa)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD7, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreIS = 57.36000   + 1.96*c(-2.345740, 2.345740)
ci.ApostIS = -44.82667 + 1.96*c(-3.813169,3.813169)

  ######################
###########################
## B pre vs B post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD9 = subset(D_tall, ! condition %in% c(1:7,9:10,12:24))
subD9$condition = factor(subD9$condition)

# run model
lme.fit.ISb = lme(measure~condition, random=~1|ID, data=subD9)

# summary of model 
summary(lme.fit.ISb)

# get ANOVA results for lme model
anova.lme(lme.fit.ISb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD9, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpreIS = 45.4100000   + 1.96*c(-1.456779, 1.456779)
ci.CpostIS = -0.4383333 + 1.96*c(-2.140663,2.140663)


  ######################
###########################
## C pre vs C post in IS ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD8 = subset(D_tall, ! condition %in% c(1:8,10:11,13:24))
subD8$condition = factor(subD8$condition)

# run model
lme.fit.ISc = lme(measure~condition, random=~1|ID, data=subD8)

# summary of model 
summary(lme.fit.ISc)

# get ANOVA results for lme model
anova.lme(lme.fit.ISc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD8, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.BpreIS = 47.77667   + 1.96*c(-1.636913, 1.636913)
ci.BpostIS = 46.41667 + 1.96*c(-2.359320,2.359320)


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
subD9 = subset(D_tall, ! condition %in% c(2:3,5:24))
subD9$condition = factor(subD9$condition)

# run model
lme.fit.BBa = lme(measure~condition, random=~1|ID, data=subD9)

# summary of model 
summary(lme.fit.BBa)

# get ANOVA results for lme model
anova.lme(lme.fit.BBa)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD9, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreBB = 49.66000   + 1.96*c(-2.073400, 2.073400)
ci.ApostBB = 45.47833 + 1.96*c(-2.928393,2.928393)


  ######################
###########################
## B pre vs B post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD10 = subset(D_tall, ! condition %in% c(1,3:4,6:24))
subD10$condition = factor(subD10$condition)

# run model
lme.fit.BBb = lme(measure~condition, random=~1|ID, data=subD10)

# summary of model 
summary(lme.fit.BBb)

# get ANOVA results for lme model
anova.lme(lme.fit.BBb)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD10, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.ApreBB = 51.84333   + 1.96*c(-1.995331, 1.995331)
ci.ApostBB = -9.18333 + 1.96*c(-2.939661,2.939661)

# one sample t-test to determine if post-rating differs
# from the rating that would be expected if participants
# were tracking the statistics
install.packages("ggpubr")
postBbb.ttest.wilcox = wilcox.test(D_tall$measure[D_tall$condition==6], mu = 50, alternative = "two.sided")
postBbb.ttest = t.test (D_tall$measure[D_tall$condition==6], mu=50)


  ######################
###########################
## C pre vs C post in BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD11 = subset(D_tall, ! condition %in% c(1:2,4:5,7:24))
subD11$condition = factor(subD11$condition)

# run model
lme.fit.BBc = lme(measure~condition, random=~1|ID, data=subD11)

# summary of model 
summary(lme.fit.BBc)

# get ANOVA results for lme model
anova.lme(lme.fit.BBc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD11, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpreBB = 44.02667   + 1.96*c(-1.839741, 1.839741)
ci.CpostBB = -1.10500 + 1.96*c(-2.957971,2.957971)


  ######################
###########################
## C post IS vs C post BB ##
###########################
  #######################

# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD11 = subset(D_tall, ! condition %in% c(1:2,4:5,7:24))
subD11$condition = factor(subD11$condition)

# run model
lme.fit.BBc = lme(measure~condition, random=~1|ID, data=subD11)

# summary of model 
summary(lme.fit.BBc)

# get ANOVA results for lme model
anova.lme(lme.fit.BBc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD11, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpreBB = 44.02667   + 1.96*c(-1.839741, 1.839741)
ci.CpostBB = -1.10500 + 1.96*c(-2.957971,2.957971)


# subset full data set to include only pre- and post-ratings 
# of A in the 1C condition
subD12 = subset(D_tall, ! condition %in% c(1:5,7:11,13:24))
subD12$condition = factor(subD12$condition)

# run model
lme.fit.BBISc = lme(measure~condition, random=~1|ID, data=subD12)

# summary of model 
summary(lme.fit.BBISc)

# get ANOVA results for lme model
anova.lme(lme.fit.BBISc)

# obtained 95% Bootstrapped CI for pre- and post-ratings of B
set.seed(2017)
bbB2.lme.Bootobj = boot(subD12, bbB.lme.fit , R=4000)
bbB2.lme.Bootobj

# obtain 95% CI from SE from the 'boot' analysis above
ci.CpostBB = 42.92167  + 1.96*c(-2.389536, 2.389536)
ci.CpostIS = 51.27167 + 1.96*c(-3.136548,3.136548)



#####
# Clustered Bar Graph
#####
F = D
F = as.data.frame(F)
fix(F)

# names: "ID"           "test.trial"   "measure"      "condition"    "period"       "test.trial.2" 

# convert data from "wide" to "tall" format
F_tall = reshape(F, varying = 3:26, v.names = "measure", timevar = "test.trial", idvar = "ID", new.row.names = 1:1440, direction = "long")
F_tall$condition = rep(c("BB","IS","1C","2C"), each = 360)
F_tall$condition = factor(F_tall$condition, levels = c("BB", "IS", "1C", "2C")) 

F_tall$period = rep(c("pre","post"), each = 180) # just take 360 and divide it by number of items in this vector, which is 2
F_tall$period = factor(F_tall$period, levels = c("pre", "post")) # you did this, as opposed to "as.factor"
# to ensure that the labels are ordered
# in THAT particular way


F_tall$test.trial.2 = rep(c("A", "B", "C"), each = 60, times = 4)
F_tall$test.trial.2 = as.factor(F_tall$test.trial.2)

F_tall = F_tall[, c("ID", "test.trial","condition", "test.trial.2","measure", "period")] # reoder column names


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