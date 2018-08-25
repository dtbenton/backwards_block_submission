########################################################
#############                              #############
#############           Initial SET UP     #############
#############                              #############
########################################################

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

shapiro.mat = matrix(NA, nrow=4, ncol=4)
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


