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
