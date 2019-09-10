##############################################################
##############################################################
##############################################################
#############                                    #############
#############      CROSS EXPERIMENT ANALYSIS     #############
#############                                    #############
##############################################################
##############################################################
##############################################################

# load all relevant libraries
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


# LOAD TWO CAUSE DATA #
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

# OBJECT COLUMN
D_tall$objects = as.factor(rep(1:2, times = 288))
D_tall$objects = revalue(x = as.factor(D_tall$objects), 
                         c("1" = "A", "2"="B"))
# REORDER COLUMNS'
D_tall$condition = NULL
D_tall$row.names = NULL
D_tall = D_tall[,c(1,2,3,5,7,6,4)]


# ADD DIFFERENCE COLUMN BETWEEN BB-B-PRE AND BB-B-MID
D_tall$bprebmid = D_tall$

# LOAD THREE CAUSE DATA #
R = read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
R = R[c(1:24),]

R_tall = reshape(R, varying = 4:39, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

R_tall$measure = as.numeric(R_tall$measure)
R_tall$sex = as.factor(R_tall$sex)
R_tall$condition = as.factor(R_tall$condition)
R_tall$condition_order = as.factor(R_tall$condition_order)

R_tall = R_tall[order(R_tall$ID),]

# ADD A CONDITION NAME COLUMN
R_tall$condition_names = as.factor(rep(1:4, each = 9, times = 24))
R_tall$condition_names = revalue(x = as.factor(R_tall$condition_names), 
                                 c("1" = "BB", "2"="IS", "3" = "1C", 
                                   "4" = "2C"))

# ADD A 'PHASE' COLUMN
R_tall$phase = as.factor(rep(1:3, each = 3, times = 96))
R_tall$phase = revalue(x = as.factor(R_tall$phase), 
                       c("1" = "Pre", "2"="Mid", "3" = "Post"))

# RENAME SEX COLUMN
R_tall$sex = revalue(x = as.factor(R_tall$sex), 
                     c("1" = "M", "2"="F"))

# OBJECT COLUMN
R_tall$objects = as.factor(rep(1:3, times = 288))
R_tall$objects = revalue(x = as.factor(R_tall$objects), 
                         c("1" = "A", "2"="B", "3"="C"))
# REORDER COLUMNS'
R_tall$condition = NULL
R_tall$row.names = NULL
R_tall$condition_order = NULL
R_tall = R_tall[,c(1,2,4,5,6,3)]


# LOAD FOUR CAUSE DATA #
V = read.csv(file.choose(), header = TRUE, stringsAsFactors = FALSE)
V = V[c(1:20),]

V_tall = reshape(V, varying = 4:51, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")

V_tall$measure = as.numeric(V_tall$measure)
V_tall$sex = as.factor(V_tall$sex)
V_tall$condition = as.factor(V_tall$condition)


V_tall = V_tall[order(V_tall$ID),]

# ADD A CONDITION NAME COLUMN
V_tall$condition_names = as.factor(rep(1:4, each = 12, times = 20))
V_tall$condition_names = revalue(x = as.factor(V_tall$condition_names), 
                                 c("1" = "BB", "2"="IS", "3" = "1C", 
                                   "4" = "2C"))


# ADD A CONDITION ORDER COLUMN
V_tall$condition_order = revalue(x = as.factor(V_tall$condition), 
                                 c("1" = "1234", "2"="2413", "3"="3142", "4"="4321"))

# ADD A 'PHASE' COLUMN
V_tall$phase = as.factor(rep(1:3, each = 4, times = 80))
V_tall$phase = revalue(x = as.factor(V_tall$phase), 
                       c("1" = "Pre", "2"="Mid", "3" = "Post"))

# RENAME SEX COLUMN
V_tall$sex = revalue(x = as.factor(V_tall$sex), 
                     c("1" = "M", "2"="F"))

# OBJECT COLUMN
V_tall$objects = as.factor(rep(1:4, times = 240))
V_tall$objects = revalue(x = as.factor(V_tall$objects), 
                         c("1" = "A", "2"="B", "3"="C", "4"="D"))
# REORDER COLUMNS'
V_tall$condition = NULL
V_tall$row.names = NULL
V_tall = V_tall[,c(1,4,6,7,3)]



# T.Tests #
# BB CONDITION #
# two-cause pre-mid vs three-cause pre-mid
twoVthree_ttest = t.test(D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Pre" & D_tall$objects=="B"]-
         D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Mid" & D_tall$objects=="B"],
       R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Pre" & R_tall$objects=="B"]-
         R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Mid" & R_tall$objects=="B"],
       alternative = "two.sided", var.equal = FALSE)


# Bayes factor for twoVthree_ttest #
N2 = data.frame (ID = 1:24, two.bbpremid = as.numeric(D$PRE_B_BB)-as.numeric(D$MID_B_BB), 
                 three.bbpremid = as.numeric(R$PRE_B_BB)-as.numeric(R$MID_B_BB))
N2$ID = as.factor(N2$ID)

N2_tall = reshape(N2, varying = 2:3, v.names = "measure", timevar = "condition", 
                 idvar = "ID", 
                 direction = "long")
N2_tall$phase_type = as.factor(rep(1:2, each = 24, times = 1))
N2_tall$phase_type = revalue(x = as.factor(N2_tall$phase_type), 
                       c("1"="two.bbpremid", "2"="three.bbpremid"))
N2_tall = N2_tall[,c(1,2,4,3)]



lm.null = lme(measure~1, random=~1|ID, data=N2_tall)
lm.alt = lme(measure~phase_type, random=~1|ID, data=N2_tall)

#obtain BICs for the null and alternative models
null.bic = BIC(lm.null)
alt.bic = BIC(lm.alt)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01 = exp((alt.bic - null.bic)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10 = 1/BF01



# TWO V. FOUR BBBPRE BBBMID #

twoVfour_ttest = t.test(D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Pre" & D_tall$objects=="B"]-
                           D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Mid" & D_tall$objects=="B"],
                         V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Pre" & V_tall$objects=="B"]-
                           V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Mid" & V_tall$objects=="B"],
                         alternative = "two.sided", var.equal = FALSE)



# Bayes factor for twoVthree_ttest #
four_preBBB = append(V$PRE_B_BB, rep(NA,4))
four_midBBB = append(V$MID_B_BB, rep(NA,4))
N3 = data.frame (ID = 1:24, two.bbpremid = as.numeric(D$PRE_B_BB)-as.numeric(D$MID_B_BB), 
                 four.bbpremid = as.numeric(four_preBBB)-as.numeric(four_midBBB))
N3$ID = as.factor(N3$ID)

N3_tall = reshape(N3, varying = 2:3, v.names = "measure", timevar = "condition", 
                  idvar = "ID", 
                  direction = "long")
N3_tall$phase_type = as.factor(rep(1:2, each = 24, times = 1))
N3_tall$phase_type = revalue(x = as.factor(N3_tall$phase_type), 
                             c("1"="two.bbpremid", "2"="four.bbpremid"))
N3_tall = na.omit(N3_tall)
N3_tall = N3_tall[,c(1,2,4,3)]



lm.null.2 = lme(measure~1, random=~1|ID, data=N3_tall)
lm.alt.2 = lme(measure~phase_type, random=~1|ID, data=N3_tall)

#obtain BICs for the null and alternative models
null.bic.2 = BIC(lm.null.2)
alt.bic.2 = BIC(lm.alt.2)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01.2 = exp((alt.bic.2 - null.bic.2)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10.2 = 1/BF01.2



#############
# OBJECT A  #
#############
# two-cause pre-mid vs three-cause pre-mid
twoVthree_ttest.bb.a = t.test(D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Pre" & D_tall$objects=="A"]-
                           D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Mid" & D_tall$objects=="A"],
                         R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Pre" & R_tall$objects=="A"]-
                           R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Mid" & R_tall$objects=="A"],
                         alternative = "two.sided", var.equal = FALSE)


# Bayes factor for twoVthree_ttest object A for BB #
N4 = data.frame (ID = 1:24, two.bbpremid = as.numeric(D$PRE_A_BB)-as.numeric(D$MID_A_BB), 
                 three.bbpremid = as.numeric(R$PRE_A_BB)-as.numeric(R$MID_A_BB))
N4$ID = as.factor(N4$ID)

N4_tall = reshape(N4, varying = 2:3, v.names = "measure", timevar = "condition", 
                  idvar = "ID", 
                  direction = "long")
N4_tall$phase_type = as.factor(rep(1:2, each = 24, times = 1))
N4_tall$phase_type = revalue(x = as.factor(N4_tall$phase_type), 
                             c("1"="two.bbpremid", "2"="three.bbpremid"))
N4_tall = N4_tall[,c(1,2,4,3)]



lm.null.3 = lme(measure~1, random=~1|ID, data=N4_tall)
lm.alt.3 = lme(measure~phase_type, random=~1|ID, data=N4_tall)

#obtain BICs for the null and alternative models
null.bic.3 = BIC(lm.null.3)
alt.bic.3 = BIC(lm.alt.3)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01.3 = exp((alt.bic.3 - null.bic.3)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10.3 = 1/BF01.3





# TWO V. FOUR BBBPRE BBBMID #

twoVfour_ttest_a_bb = t.test(D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Pre" & D_tall$objects=="A"]-
                          D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Mid" & D_tall$objects=="A"],
                        V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Pre" & V_tall$objects=="A"]-
                          V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Mid" & V_tall$objects=="A"],
                        alternative = "two.sided", var.equal = FALSE)



# Bayes factor for twoVthree_ttest #
four_preABB = append(V$PRE_A_BB, rep(NA,4))
four_midABB = append(V$MID_A_BB, rep(NA,4))
N5 = data.frame (ID = 1:24, two.bbpremid = as.numeric(D$PRE_A_BB)-as.numeric(D$MID_A_BB), 
                 four.bbpremid = as.numeric(four_preABB)-as.numeric(four_midABB))
N5$ID = as.factor(N5$ID)

N5_tall = reshape(N5, varying = 2:3, v.names = "measure", timevar = "condition", 
                  idvar = "ID", 
                  direction = "long")
N5_tall$phase_type = as.factor(rep(1:2, each = 24, times = 1))
N5_tall$phase_type = revalue(x = as.factor(N5_tall$phase_type), 
                             c("1"="two.bbpremid", "2"="four.bbpremid"))
N5_tall = na.omit(N5_tall)
N5_tall = N5_tall[,c(1,2,4,3)]



lm.null.4 = lme(measure~1, random=~1|ID, data=N5_tall)
lm.alt.4 = lme(measure~phase_type, random=~1|ID, data=N5_tall)

#obtain BICs for the null and alternative models
null.bic.4 = BIC(lm.null.4)
alt.bic.4 = BIC(lm.alt.4)

# compute the BF01  - this is the BF whose value is interpreted as the evidence in favor of the null (e.g., if the BF01 = 2.6, this means that there is 2.6 times as much evidence for the null than for the alternative or the evidence is 2.6:1 in favor of the null)

BF01.4 = exp((alt.bic.4 - null.bic.4)/2) # this yields a BF that is interpreted as the evidence in favor of the null; it's critical that the alt.bic comes first otherwise your interpretation of the resulting BF value will be incorrect
BF10.4 = 1/BF01.4



#############################
# RECOMPUTED BAYES' FACTORS #
#############################
BBbpremidExp1 = D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Pre" & D_tall$objects=="B"]-
  D_tall$measure[D_tall$condition_names=="BB" & D_tall$phase=="Mid" & D_tall$objects=="B"]

BBbpremidExp2 = R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Pre" & R_tall$objects=="B"]-
  R_tall$measure[R_tall$condition_names=="BB" & R_tall$phase=="Mid" & R_tall$objects=="B"]

BBbpremidExp3 = V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Pre" & V_tall$objects=="B"]-
  V_tall$measure[V_tall$condition_names=="BB" & V_tall$phase=="Mid" & V_tall$objects=="B"]

# exp1 vs exp2
ttestBF(x=BBbpremidExp1,y=BBbpremidExp2,paired=FALSE)

# exp1 vs exp3
ttestBF(x=BBbpremidExp1,y=BBbpremidExp3,paired=FALSE)

