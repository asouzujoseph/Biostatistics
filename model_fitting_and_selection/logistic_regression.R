## Logistic regression

# title: "Model fitting and selection using logistic regression"
# author: "Nnamdi Joseph Asouzu"
# date: "8/22/2021"

## Case study
# Low birth weight data: The goal of this study was to identify risk factors associated with
# giving birth to a low birth weight baby (weighing less than 2500 grams). Data were collected
# on 189 women, 59 of which had low birth weight babies and 130 of which had normal birth
# weight babies. Four variables which were thought to be of importance were age, weight of
# the subject at her last menstrual period, race, and the number of physician visits during the
# first trimester of pregnancy.

## ID and their meaning
# LOW;Low BirthWeight (0 = BirthWeight >= 2500g, 1 = Birth Weight < 2500g)
# AGE;Age of the Mother in Years 
# LWT;Weight in Pounds at the Last Menstrual Period 
# RACE;Race (1 = White, 2 = Black, 3 = Other) 
# SMOKE;Smoking Status During Pregnancy (1 = Yes, 0 = No) 
# PTL;History of Premature Labor (0 = None 1 = One, etc.) 
# HT;History of Hypertension (1 = Yes, 0 = No) 
# UI;Presence of Uterine Irritability (1 = Yes, 0 = No) 
# FTV;Number of Physician Visits During the First Trimester (0 = None, 1 = One, 2 = Two, etc.)
# BWT;Birth Weight in Grams

# Use the continuous variable LOW as the response variable.
# Analysis 1
# Considering the four most clinically important variables previously discussed, i.e., age,
# weight of the subject at her last menstrual period, race, and the number of physician
# visits during the first trimester of pregnancy, carry out a model building exercise to
# determine the best model to describe the data. Based your model building exercise on
# the AIC and AIC weights. Discuss your results carefully. 

# Use backward elimination procedure to decide which predictor variables can be dropped
# from the regression model. Control the type I error at alpha <= 0.10 at each stage. Which
# variables are retained? How does this compare to the results obtained with the AIC? 
  
## Analysis 2:
#Consider now the 8 covariates included in the data set (all except Birth Weight in Grams). 
# Repeat the previous model building exercise using the AIC and the AIC weights. 
# Compare the results of both analysis and make a careful discussion.

## Use backward elimination procedure to decide which predictor variables can be dropped
#from the regression model. Control the type I error at alpha <= 0.10 at each stage. Which
#variables are retained? How does this compare to the results obtained with the AIC? 
  
  
## Using the final model selected in the previous model building exercise (based on the AIC)
# answer the following scientific questions 
# What scientific insight does your final model offer? 
# Interpret the final model using the estimated coefficients and odd ratios. 

#load the necessary libraries
suppressWarnings(suppressMessages(library("pastecs")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("PerformanceAnalytics")))
suppressWarnings(suppressMessages(library("AICcmodavg")))
suppressWarnings(suppressMessages(library("ROCR")))
suppressWarnings(suppressMessages(library("caret")))

# import dataset
df <- read.table(file="C:/Users/Nnamdi/Desktop/Bioinformatics/Statistical_methods_for_bioinformatics/part_I/Projects/Project2/lowbwt.dat", header=T)
head(df)

# Descriptive statistics
desc.df <- stat.desc(df[,c("AGE","LWT","RACE","SMOKE","PTL","HT","UI","FTV")], basic =T, desc=T)
desc.df
df$RACE <- as.factor(df$RACE)
df$FTV <- as.factor(df$FTV)
df$SMOKE <- as.factor(df$SMOKE)
df$PTL <- as.factor(df$PTL)
df$HT <- as.factor(df$HT)
df$UI <- as.factor(df$UI)
df$LOW <- as.factor(df$LOW)

## ANALYSIS 1
# Exploratory data analysis
ggplot(df,aes(LOW,AGE,color = LOW))+geom_boxplot() + theme_classic()
ggplot(df,aes(LOW,LWT,color = LOW))+geom_boxplot() + theme_classic()
ggplot(df)+geom_bar(aes(x=RACE, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=FTV, fill=as.factor(LOW)))+ theme_classic()

# Forward selection using AICc score
mod.list <- list()
mod.list[[1]] <- glm(LOW ~ AGE+LWT+RACE+FTV, data=df, family=binomial(link="logit"))
mod.list[[2]]<- glm(LOW ~ AGE*FTV+AGE*LWT+AGE*RACE+LWT*RACE+LWT*FTV+RACE*FTV, data = df,family=binomial(link="logit"))
mod.list[[3]] <- glm(LOW ~ LWT+RACE, data=df, family=binomial(link="logit"))
mod.list[[4]] <- glm(LOW ~ AGE+LWT+RACE, data=df, family=binomial(link="logit"))
mod.list[[5]] <- glm(LOW ~ FTV+RACE+AGE, data=df, family=binomial(link="logit"))
mod.list[[6]] <- glm(LOW ~ AGE+LWT+RACE+LWT*RACE, data=df, family=binomial(link="logit"))
mod.aictab=aictab(cand.set = mod.list) 
mod.aictab

# The model fitted by forward selection technique using AICc score is LOW ~ AGE + LWT + RACE

## Backward selection using p-value<0.10
bk.mod <- glm(LOW~AGE+LWT+RACE+FTV, data = df, family=binomial(link="logit"))
summary(bk.mod)

#Remove variables with p-value above the cut-off score
bk.mod2 <- update(bk.mod, . ~ . - FTV, data = df)
summary(bk.mod2)

# The model fitted by backward selection technique using control p-value less than 0.10 is LOW ~ AGE + LWT + RACE

## Analysis 2
# Exploratory data analysis
ggplot(df,aes(LOW,AGE,color = LOW))+geom_boxplot() + theme_classic()
ggplot(df,aes(LOW,LWT,color = LOW))+geom_boxplot() + theme_classic()
ggplot(df)+geom_bar(aes(x=RACE, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=FTV, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=SMOKE, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=PTL, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=HT, fill=as.factor(LOW)))+ theme_classic()
ggplot(df)+geom_bar(aes(x=UI, fill=as.factor(LOW)))+ theme_classic()

## Model fitting using forward selection technique and AICc
mod.list2 <- list()

mod.list2[[1]] <- glm(LOW ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI+FTV, data=df, family=binomial(link="logit"))
mod.list2[[2]]<- glm(LOW ~ LWT+RACE+SMOKE+PTL+HT+UI, data = df, family=binomial(link="logit"))
mod.list2[[3]] <- glm(LOW ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI, data = df, family=binomial(link="logit"))
mod.list2[[4]] <- glm(LOW ~ LWT*RACE+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI, data=df, family=binomial(link="logit"))
mod.list2[[5]] <- glm(LOW ~ AGE+RACE*HT+RACE*UI+LWT*RACE+SMOKE*PTL+SMOKE*HT+SMOKE*RACE, data=df, family=binomial(link="logit"))
mod.list2[[6]] <- glm(LOW ~ LWT+RACE+SMOKE+PTL+HT+UI+RACE*UI, data=df, family=binomial(link="logit"))

mod.aictab2=aictab(cand.set = mod.list2) 
mod.aictab2

# The model fitted by forward selection technique using AICc score is LOW ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI

## Backward selection using p-value < 0.10
analysis2_mod <- glm(LOW ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI+FTV, data=df, family=binomial(link="logit")) # full model
summary(analysis2_mod)

# Drop variables with p-value > 0.10
clean.mod <- update(analysis2_mod, . ~ . - FTV, data = df)
summary(clean.mod)

# The model fitted by backward selection technique using control p-value less than 0.10 is: LOW ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI

## The final model using AICc criteria is LOW ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI.
final_model = glm(LOW ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI, data=df, family=binomial(link="logit"))
summary(final_model)

## Model is LOW ~ 0.85 - 0.04(Age) - 0.02(LWT) + 1.15(RACE2) + 0.74(RACE3) + 0.86(SMOKE1) + 1.58(PTL1) + 1.87 (HT1) - 0.86(UT1)

# Interpretation - odds ratio
exp(final_model$coefficients)

# For each unit increase in maternal age, the odds of giving birth to a baby with low birth weight decreases by 4%.
# For each unit increase in the weight of the mother before her last menstrual period, the odds of giving birth to a baby with low birth weight decreases by 1.6%. 
# The odds of giving birth to a baby with low birth weight is 3.15 times higher if the mother is of black ancestry. 
# The odds of giving birth to a baby with low birth weight is 2.09 times higher if the mother is of other ethnicity excluding black or white ancestry. 
# The odds of giving birth to a baby with low birth weight is 2.37 times higher if the mother smoked during pregnancy.
# The odds of giving birth to a baby with low birth weight is 4.84 times higher if the mother had one history of premature labor term. 
# The odds of giving birth to a baby with low birth weight is 6.47 times higher if the mother has hypertension . 
# The odds of giving birth to a baby with low birth weight is 2.35 times higher if the mother has uterine irritability.

## Model accuracy and predictions
# Create training and test data sets
set.seed(1110)
df_split = sort(sample(nrow(df), nrow(df)*0.8)) ## 80% of the dataset randomly selected
train<-df[df_split,]
test<-df[-df_split,]
train_mod <- glm(LOW ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI,family = binomial(link = "logit"), data = train)

# predictions
pred <- predict(train_mod, newdata = test, type = "response")
results <- ifelse(pred > 0.5,1,0)
table(results,test$LOW)
Accuracy.logistic <- round(mean(results == test$LOW), digits = 2)*100
print(paste('The accuracy of the model is ',Accuracy.logistic,"%"))

predict <- fitted(train_mod)
pred <- prediction(predict, train$LOW)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, main="sensitivity vs false positive rate",colorize=TRUE)
abline(a = 0, b = 1)

# ROCR and AUC
auc_ROCR <- performance(pred, measure = "auc")
auc_ROCR <- auc_ROCR@y.values[[1]]
auc_ROCR











