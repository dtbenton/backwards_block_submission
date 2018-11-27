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

D_tall = reshape(D, varying = 4:27, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

D_tall$measure = as.numeric(D_tall$measure)
D_tall$sex = as.factor(D_tall$sex)
D_tall$condition_order = as.factor(D_tall$condition_order)
D_tall$condition = as.factor(D_tall$condition)

D_tall = D_tall[order(D_tall$ID),]

# ADD A CONDITION NAME COLUMN
D_tall$condition_names = as.factor(rep(1:4, each = 6, times = 24))
D_tall$condition_names = revalue(x = as.factor(D_tall$condition_names), 
                           c("1" = "BB", "2"="IS", "3" = "1C", 
                             "4" = "2C"))

# ADD A 'PHASE' COLUMN
D_tall$phase = as.factor(rep(1:3, each = 2, times = 96))
D_tall$phase = revalue(x = as.factor(D_tall$phase), 
                                 c("1" = "Pre", "2"="Mid", "3" = "Post"))

# RENAME SEX COLUMN
D_tall$sex = revalue(x = as.factor(D_tall$sex), 
                       c("1" = "M", "2"="F"))

# REORDER COLUMNS'
D_tall = D_tall[,c(1,2,3,6,7,5,4)]
D_tall$condition = NULL
D_tall$row.names = NULL



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
