              #################################
###########################################################
        # BACKWARDS-BLOCKING NEURAL-NETWORK SCRIPT #
###########################################################
              #################################

# This script can be used to generate a pretraining, habituation, and
# testing set. Critically, the pretraining set for this script includes
# preventative, causal, AND ambiguous cases.


# INSTALL 'COMBOS' PACKAGE 
require(hier.part)


#############              
## OBJECTS ##
#############
objects_full = combos(12)$binary



# A CAUSES
a_objects_stim = as.data.frame(matrix(NA, nrow=100, ncol=12))
for(i in 1:nrow(objects_full)){
  if(objects_full[i,4]==1|objects_full[i,5]==1|objects_full[i,6]==1|objects_full[i,7]==1|objects_full[i,8]==1|objects_full[i,9]==1|objects_full[i,10]==1|objects_full[i,11]==1|objects_full[i,12]==1){
    a_objects_stim[i,] = NA
  } else {
    a_objects_stim[i,] = objects_full[i,]
  }
}
a_objects_stim = na.omit(a_objects_stim)
names(a_objects_stim) = NULL
rownames(a_objects_stim) = c(1:7)





# B CAUSES
b_objects_stim = as.data.frame(matrix(NA, nrow=100, ncol=12))
for(i in 1:nrow(objects_full)){
  if(objects_full[i,1]==1|objects_full[i,2]==1|objects_full[i,3]==1|objects_full[i,7]==1|objects_full[i,8]==1|objects_full[i,9]==1|objects_full[i,10]==1|objects_full[i,11]==1|objects_full[i,12]==1){
    b_objects_stim[i,] = NA
  } else {
    b_objects_stim[i,] = objects_full[i,]
  }
}

b_objects_stim = na.omit(b_objects_stim)
names(b_objects_stim) = NULL
rownames(b_objects_stim) = c(1:7)





# AMBIGUOUS CAUSES 1
ambig_objects_stim_1 = as.data.frame(matrix(NA, nrow=100, ncol=12))
for(i in 1:nrow(objects_full)){
  if(objects_full[i,1]==1|objects_full[i,2]==1|objects_full[i,3]==1|objects_full[i,4]==1|objects_full[i,5]==1|objects_full[i,6]==1|objects_full[i,10]==1|objects_full[i,11]==1|objects_full[i,12]==1){
    ambig_objects_stim_1[i,] = NA
  } else {
    ambig_objects_stim_1[i,] = objects_full[i,]
  }
}

ambig_objects_stim_1 = na.omit(ambig_objects_stim_1)
names(ambig_objects_stim_1) = NULL
rownames(ambig_objects_stim_1) = c(1:7)





# AMBIGUOUS CAUSES 2
ambig_objects_stim_2 = as.data.frame(matrix(NA, nrow=100, ncol=12))
for(i in 1:nrow(objects_full)){
  if(objects_full[i,1]==1|objects_full[i,2]==1|objects_full[i,3]==1|objects_full[i,4]==1|objects_full[i,5]==1|objects_full[i,6]==1|objects_full[i,7]==1|objects_full[i,8]==1|objects_full[i,9]==1){
    ambig_objects_stim_2[i,] = NA
  } else {
    ambig_objects_stim_2[i,] = objects_full[i,]
  }
}

ambig_objects_stim_2 = na.omit(ambig_objects_stim_2)
names(ambig_objects_stim_2) = NULL
rownames(ambig_objects_stim_2) = c(1:15)


# ABSENT CAUSE
absent_cause = data.frame(x = c('0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'))
names(absent_cause) = NULL
rownames(absent_cause) = NULL


# AMBIGUOUS OBJECTS
ambig_objects = data.frame(x = c('1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0', 
                                 '0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'))
rownames(ambig_objects) = NULL
names(ambig_objects) = NULL



# OUTCOME EVENTS 
outcomes = data.frame(x = c('1 0','0 1'))
names(outcomes) = NULL

# Detector ON: 1 0
# Detector OFF: 0 1



# POSITION
position = data.frame(x = c('1 0','0 1'))

# Position ON: 1 0
# Position OFF: 0 1