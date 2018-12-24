library(tidyr)
library(dplyr)
library(ggplot2)
library(plyr)
library(zoo)
library(ROCR)
library(C50)

##Model Creation and Prediction##

#Create new column for prediction
ProteinComparison_df$Type2 <- ProteinComparison_df$Type #create new column Type2
ProteinComparison_df$Type2 <- gsub("Transmembrane", "0", ProteinComparison_df$Type2) #assign 0 to transmembrane proteins
ProteinComparison_df$Type2 <- gsub("Non-membrane", "1", ProteinComparison_df$Type2) #assign 1 to non-membrane proteins
ProteinComparison_df$Type2 <- factor(ProteinComparison_df$Type2) #make Type2 column factors

#Split final data frame into training and testing sets
index <- sample(1:nrow(ProteinComparison_df),36) #split final data frame randomly 
train_df <- ProteinComparison_df[index, ] #place ~70% of data frame into training set
test_df <- ProteinComparison_df[-index, ] #place ~30% of data frame into testing set

#Create logistic regression model
model1 <- glm(Type2 ~ Average_Hydrophobic + Lys + Met + Trp, family = "binomial", data = train_df)
summary(model1)
anova(model1, test = "Chisq") #measure of significance for each feature

#Model prediction with ROC curve and auc
prob <- predict(model1, newdata = test_df, type = "response")
pred <- prediction(prob, test_df$Type2)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf) #plots ROC curve
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #value represents percentage of accurately predicted proteins

#Create Decision Tree model
train_tree <- C5.0(Type2 ~ Average_Hydrophobic + Met + Trp, data = train_df)
plot(train_tree)

#Decision Tree model predicition 
results <- predict(object = train_tree, newdata = test_df, type = "class")
table(results, test_df$Type2)