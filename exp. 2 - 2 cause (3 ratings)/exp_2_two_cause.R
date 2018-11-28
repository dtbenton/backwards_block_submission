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

# OBJECT COLUMN
D_tall$objects = as.factor(rep(1:2, times = 288))
D_tall$objects = revalue(x = as.factor(D_tall$objects), 
                         c("1" = "A", "2"="B"))
# REORDER COLUMNS'
D_tall$condition = NULL
D_tall$row.names = NULL
D_tall = D_tall[,c(1,2,3,5,7,6,4)]




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
# create a data frame in which the 1C condition is subsetted
1C_subset = subset(D_tall, ! condition_names %in% c("BB","IS","2C"))

# 1C condition
lme.fit.1C = lme(measure~(condition_names+phase+objects)^3, 
                 random=~1|ID, 
                 data=D_tall)

lme.fit.1C = lme(measure[condition_names=="1C"]~condition_names=="1C", 
                 random=~1|ID, 
                 data=D_tall)

anova.lme(lme.fit.1C)



practice = lme(D_tall$measure[D_tall$condition_names=="1C"]~D_tall$condition_names[D_tall$condition_names=="1C"]+
                 D_tall$objects[D_tall$condition_names=="1C"],
               random=~1|ID,
               data=D_tall)
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
