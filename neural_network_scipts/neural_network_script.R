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





# C CAUSES 
c_objects_stim = as.data.frame(matrix(NA, nrow=100, ncol=12))
for(i in 1:nrow(objects_full)){
  if(objects_full[i,1]==1|objects_full[i,2]==1|objects_full[i,3]==1|objects_full[i,4]==1|objects_full[i,5]==1|objects_full[i,6]==1|objects_full[i,10]==1|objects_full[i,11]==1|objects_full[i,12]==1){
    c_objects_stim[i,] = NA
  } else {
    c_objects_stim[i,] = objects_full[i,]
  }
}

c_objects_stim = na.omit(c_objects_stim)
names(c_objects_stim) = NULL
rownames(c_objects_stim) = c(1:7)



# ABSENT CAUSE
absent_cause = data.frame(x = c('0 0 0 0 0 0 0 0 0 0 0 0'))
names(absent_cause) = NULL
rownames(absent_cause) = NULL


# AMBIGUOUS OBJECTS
ambig_objects = data.frame(x = c('1 1 1 1 1 1 1 1 1 1 1 0', 
                                 '0 1 1 1 1 1 1 1 1 1 1 1',
                                 '1 1 1 1 1 1 1 1 1 1 1 1'))
rownames(ambig_objects) = NULL
names(ambig_objects) = NULL



# OUTCOME EVENTS 
outcomes = data.frame(x = c('1 0','0 1'))
names(outcomes) = NULL

# Detector ON: 1 0
# Detector OFF: 0 1



# POSITION
position = data.frame(x = c('1 0','0 1'))
names(position) = NULL

# Position (w.r.t. the machine) ON: 1 0
# Position (w.r.t. the machine) OFF: 0 1

# PRETRAINING
sink('pretrain_backwards_blocking_new.ex')
for(i in 1:nrow(a_objects_stim)){
  # A ON; B ON; C OFF #
  cat(paste("name:","AonBonCoff",rownames(a_objects_stim)[i], "\n", sep=""))
  cat(paste("I:", "\n", sep="\t"))
  # Object A
  cat(paste("(Object_A)", sep="\t"))
  print(a_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_A)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object B
  cat(paste("(Object_B)", sep="\t"))
  print(b_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_B)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object C
  cat(paste("(Object_C)", sep="\t"))
  print(c_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_C)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Input_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(paste("T:", "\n", sep="\t"))
  cat(paste("(Output_Activation)", sep="\t"))
  print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste(";", sep="\t"))
  cat("\n")
  
  
  # A ON; B OFF; C OFF #
  cat(paste("name:","AonBoffCoff",rownames(a_objects_stim)[i], "\n", sep=""))
  cat(paste("I:", "\n", sep="\t"))
  # Object A
  cat(paste("(Object_A)", sep="\t"))
  print(a_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_A)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object B
  cat(paste("(Object_B)", sep="\t"))
  print(b_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_B)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object C
  cat(paste("(Object_C)", sep="\t"))
  print(c_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_C)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Input_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(paste("T:", "\n", sep="\t"))
  cat(paste("(Output_Activation)", sep="\t"))
  print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste(";", sep="\t"))
  cat("\n")
  
  
  # A OFF; B ON; C OFF #
  cat(paste("name:","AoffBonCoff",rownames(a_objects_stim)[i], "\n", sep=""))
  cat(paste("I:", "\n", sep="\t"))
  # Object A
  cat(paste("(Object_A)", sep="\t"))
  print(a_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_A)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object B
  cat(paste("(Object_B)", sep="\t"))
  print(b_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_B)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object C
  cat(paste("(Object_C)", sep="\t"))
  print(c_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_C)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Input_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(paste("T:", "\n", sep="\t"))
  cat(paste("(Output_Activation)", sep="\t"))
  print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste(";", sep="\t"))
  cat("\n")

  
  
  # A OFF; B OFF; C ON #
  cat(paste("name:","AoffBoffCon",rownames(a_objects_stim)[i], "\n", sep=""))
  cat(paste("I:", "\n", sep="\t"))
  # Object A
  cat(paste("(Object_A)", sep="\t"))
  print(a_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_A)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object B
  cat(paste("(Object_B)", sep="\t"))
  print(b_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_B)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object C
  cat(paste("(Object_C)", sep="\t"))
  print(c_objects_stim[i,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_C)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Input_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(paste("T:", "\n", sep="\t"))
  cat(paste("(Output_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste(";", sep="\t"))
  cat("\n")
}
sink()




# AB+ #
sink('AB_backwards_blocking.ex')
  # A ON; B ON; C OFF #
  cat(paste("name:","AonBonCoff",rownames(a_objects_stim)[i], "\n", sep=""))
  cat(paste("I:", "\n", sep="\t"))
  # Object A
  cat(paste("(Object_A)", sep="\t"))
  print(ambig_objects[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_A)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object B
  cat(paste("(Object_B)", sep="\t"))
  print(ambig_objects[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_B)", sep="\t"))
  print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
  # Object C
  cat(paste("(Object_C)", sep="\t"))
  print(ambig_objects[3,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Position_C)", sep="\t"))
  print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("(Input_Activation)", sep="\t"))
  print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(paste("T:", "\n", sep="\t"))
  cat(paste("(Output_Activation)", sep="\t"))
  print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste(";", sep="\t"))
  cat("\n")
sink()



# A+ #
sink('A_backwards_blocking.ex')
# A ON; B OFF; C OFF #
cat(paste("name:","AonBonCoff",rownames(a_objects_stim)[i], "\n", sep=""))
cat(paste("I:", "\n", sep="\t"))
# Object A
cat(paste("(Object_A)", sep="\t"))
print(ambig_objects[1,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_A)", sep="\t"))
print(position[1,], sep = "\t", quote = FALSE, row.names = FALSE)
# Object B
cat(paste("(Object_B)", sep="\t"))
print(ambig_objects[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_B)", sep="\t"))
print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
# Object C
cat(paste("(Object_C)", sep="\t"))
print(ambig_objects[3,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_C)", sep="\t"))
print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Input_Activation)", sep="\t"))
print(outcomes[2,1], sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")

cat(paste("T:", "\n", sep="\t"))
cat(paste("(Output_Activation)", sep="\t"))
print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste(";", sep="\t"))
cat("\n")
sink()




# Test #
sink('test_backwards_blocking.ex')
cat(paste("defI:-", "\n", sep="\t"))
# A #
cat(paste("name:","A", "\n", sep=""))
cat(paste("I:", "\n", sep="\t"))
# Object A
cat(paste("(Object_A)", sep="\t"))
print(ambig_objects[1,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_A)", sep="\t"))
print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")

cat(paste("T:", "\n", sep="\t"))
cat(paste("(Output_Activation)", sep="\t"))
print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste(";", sep="\t"))
cat("\n")


# B #
cat(paste("name:","B", "\n", sep=""))
cat(paste("I:", "\n", sep="\t"))
# Object B
cat(paste("(Object_B)", sep="\t"))
print(ambig_objects[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_B)", sep="\t"))
print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")

cat(paste("T:", "\n", sep="\t"))
cat(paste("(Output_Activation)", sep="\t"))
print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste(";", sep="\t"))
cat("\n")


# C #
cat(paste("name:","C","\n", sep=""))
cat(paste("I:", "\n", sep="\t"))
# Object C
cat(paste("(Object_C)", sep="\t"))
print(ambig_objects[3,], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste("(Position_C)", sep="\t"))
print(position[2,], sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n")

cat(paste("T:", "\n", sep="\t"))
cat(paste("(Output_Activation)", sep="\t"))
print(outcomes[1,1], sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste(";", sep="\t"))
cat("\n")
sink()