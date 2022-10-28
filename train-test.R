#########################################
##
##    Script for splitting data into
##     Training and Testing Datasets
##
#########################################


# Load Libraries
library(caret)

base_seed <- 719921  # generated using random.org
set.seed(base_seed)  # set seed for reproducibility

# Separate data to training (75%) and testing (25%) sets 
SVAD_US <- read.csv('data/SVAD_US.csv')  # import US data
train_idx <- createDataPartition(SVAD_US$SEX, p=0.75, list=F)  # part data
train_US <- SVAD_US[train_idx, ]  # training set
test_US <- SVAD_US[-train_idx, ]  # testing set

write.csv(train_US, "data/SVAD_US_train.csv")  # save training set
write.csv(test_US, "data/SVAD_US_test.csv")  # save testing set
