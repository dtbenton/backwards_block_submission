# note: even though I computed for # and power, you'll only want to report on the # of participants that are required
# to reach power of 0.8; in other words, compute number when power is given.

###=============power analysis for demographic study=============###

# install the power analysis package
install.packages("pwr")
library(pwr)

########################################################
#############                              #############
#############     College Power analyses   #############
#############                              #############
########################################################

# run a power analysis to determine # of subjects to detect differences in rating (should one exist)
# given a comparison between two different majors
col_pwr_n.les = pwr.anova.test(k = 6, n = NULL , f = 0.5, sig.level = 0.05, power = 0.8) # if at fewest 1 argument is kept blank, R will solve for that argument given the values of the other arguments. "les" stands for large effect size
col_pwr_n.mes = pwr.anova.test(k = 6, n = NULL , f = 0.25, sig.level = 0.05, power = 0.8)

##### # of subjects pwr analysis summary ##################
Given the current design, I would have needed roughly 64  
subjects in each of the 6 colleges condition to have an 80% 
chance of detecting a reliable difference a reliable difference 
should one exist, where the effect size is 
conventionally 0.5 (considered large), at alpha = 0.05. Given that there are
7 reported colleges among the sample of students tested,
this analysis suggests that I would have needed to collect
384 subjects (64 per group, i.e., 6 x 64) to have enough
signal to detect an effect. Because of this, it does not 
make much sense to test the effect of "college" on causal
rating. Note that I would need to run just 35 subjects to detect
a difference if we use a medium-effect-size cut-off (of 0.25). Either way,
my study is underpowerd, as some groups (i.e., colleges) have as few as 
1 subject in them.
###########################################################




# run a power analysis to determine how much power the study currently has
# given the "college" parameter
col_pwr_pwr.les = pwr.anova.test(k = 6, n = 9, f = 0.5, sig.level = 0.05, power = )
col_pwr_pwr.mes = pwr.anova.test(k = 6, n = 9, f = 0.25, sig.level = 0.05, power = )
##### # of subjects pwr analysis summary ##################
Based on this particular power analysis and (roughly) on the number of
students assigned to a particular college in my data set, my study currently
has a 77% and 22% chance of detecting a large and medium sized effect for the 
effect of "college" on causal ratings (if this relationship is real).
More generally, given the parameters of the study, I have an 16% chance of detecting an effect. This
warrants the conclusion that my study is severely underpowered at least with regard to the factor of college. 
###########################################################



########################################################
#############                              #############
############# College Major Power analyses #############
#############                              #############
########################################################

# determine power given # of subjects
major_pwr_pwr.les = pwr.anova.test(k = 19, n = 3, f = 0.5, sig.level = 0.05, power = )
major_pwr_pwr.mes = pwr.anova.test(k = 19, n = 3, f = 0.25, sig.level = 0.05, power = )

##### # of subjects pwr analysis summary ##################
This analysis revealed that given my current study design, I have
roughly a 43% chance of detecting an effect should one exist, the
effect being that college major affects causal ratings.
###########################################################


# determine # of subjects
major_pwr_num.les = pwr.anova.test(k = 19, n = , f = 0.5, sig.level = 0.05, power = 0.8)
major_pwr_num.mes = pwr.anova.test(k = 19, n = , f = 0.25, sig.level = 0.05, power = 0.8)
##### # of subjects pwr analysis summary ##################
This analysis revealed that we would need at fewest 5
subjects per major to detect a reliable effect (should one exist)
between major and causal ratings. This is not true in the current case (because of majors 
                                                                        have as few as 1
                                                                        subject),
which supports the decision not to use this paramter when estimating a model
###########################################################



########################################################
#############                              #############
#############       Year in school         #############
#############                              #############
########################################################

year_pwr_num.les = pwr.anova.test(k = 4, n = , f = 0.5, sig.level = 0.05, power = 0.8)
year_pwr_num.mes = pwr.anova.test(k = 4, n = , f = 0.25, sig.level = 0.05, power = 0.8)

##### # of subjects pwr analysis summary ##################
This analysis revealed that I would need 11 and 44 subjects 
to detect a large and medium effect, respectively in which 
year in school affects causal rating. This is not the case 
here given that some groups (i.e., year) have as few as 3
participants (e.g., 3 juniors in the sample)
###########################################################