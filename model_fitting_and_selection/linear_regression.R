## Linear regression ###

#title: "Model fitting and selection using linear regression"
#author: "Nnamdi Joseph Asouzu"
#date: "8/22/2021"

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


# Use the continuous variable BWT as the response variable.
  
## Analysis 1:
# Considering the four most clinically important variables previously discussed, i.e., age,
# weight of the subject at her last menstrual period, race, and the number of physician
# visits during the first trimester of pregnancy, carry out a model building exercise to
# determine the best model to describe the data. Based your model building exercise on
# the AIC and AIC weights. Discuss your results carefully. 

# Use backward elimination procedure to decide which predictor variables can be dropped
#from the regression model. Control the type I error at alpha = 0.10 at each stage. Which
#variables are retained? How does this compare to the results obtained with the AIC? 
  
## Analysis 2:
# Consider now the 8 covariates included in the data set (all except Birth Weight in Grams). Repeat the previous model building exercise using the AIC and the AIC
# weights. Compare the results of both analysis and make a careful discussion.
# Use backward elimination procedure to decide which predictor variables can be dropped
# from the regression model. Control the type I error at alpha = 0.10 at each stage. Which
# variables are retained? How does this compare to the results obtained with the AIC? 

suppressWarnings(suppressMessages(library("pastecs")))
suppressWarnings(suppressMessages(library("tidyverse")))
suppressWarnings(suppressMessages(library("PerformanceAnalytics")))
suppressWarnings(suppressMessages(library("AICcmodavg")))
df <- read.table(file="C:/Users/Nnamdi/Desktop/Bioinformatics/Statistical_methods_for_bioinformatics/part_I/Projects/Project2/lowbwt.dat", header=T)
head(df)

## Analysis 1 
# descriptive statistics
desc.df <- stat.desc(df[,c("AGE","LWT","RACE","SMOKE","PTL","HT","UI","FTV")], basic =T, desc=T)
desc.df
df$RACE <- as.factor(df$RACE)
df$FTV <- as.factor(df$FTV)
df$SMOKE <- as.factor(df$SMOKE)
df$PTL <- as.factor(df$PTL)
df$HT <- as.factor(df$HT)
df$UI <- as.factor(df$UI)
df$LOW <- as.factor(df$LOW)

# exploratory data analysis
chart.Correlation(df[,c(3,4,11)], histogram=TRUE, pch=19)
ggplot(df,aes(FTV,BWT,color = FTV))+geom_boxplot() + theme_classic()
ggplot(df,aes(RACE,BWT,color = RACE))+geom_boxplot() + theme_classic()

## fitting models
mod <- lm(BWT ~ AGE+LWT+RACE+FTV, data=df)
summary (mod)

## Forward Selection using AICc score
mod.list <- list()
mod.list[[1]] <- lm(BWT ~ AGE+LWT+RACE+FTV, data=df)
mod.list[[2]]<- lm(BWT ~ AGE*FTV+AGE*LWT+AGE*RACE+LWT*RACE+LWT*FTV+RACE*FTV, data = df)
mod.list[[3]] <- lm(BWT ~ AGE+LWT+RACE, data=df)
mod.list[[4]] <- lm(BWT ~ AGE+FTV+RACE, data=df)
mod.list[[5]] <- lm(BWT ~ FTV+LWT+RACE, data=df)
mod.list[[6]] <- lm(BWT ~ AGE+LWT*RACE, data=df)
mod.list[[7]] <- lm(BWT ~ AGE+FTV*RACE, data=df)

mod.aictab=aictab(cand.set = mod.list) 
mod.aictab

# The model fitted by forward selection technique using the lowest AICc score is **BWT ~ AGE + LWT + RACE

## Backward selection using cut-off p-value < 0.10
df.mod <- lm(BWT~AGE+LWT+RACE+FTV, data = df)
summary(df.mod)

# Remove the variable with p-value above the cut off score
df.mod2 <- update(df.mod, . ~ . - FTV, data = df)
summary(df.mod2)

# Age has a p -value greater than 0.1 but we will include it in the model because it is known that maternal age affects birth weight of babies ( a condition known as Macrosomia) 
# The model fitted by backward selection technique using control p-value less than 0.10 is **BWT ~ AGE + LWT + RACE**. This model is same as the model fitted using the forward selection.  

### Analysis 2
# Exploratory data analysis
chart.Correlation(df[,-c(1,2,5:10)], histogram=TRUE, pch=19)
ggplot(df,aes(RACE,BWT,color = RACE))+geom_boxplot() + theme_classic()
ggplot(df,aes(SMOKE,BWT,color = SMOKE))+geom_boxplot() + theme_classic()
ggplot(df,aes(PTL,BWT,color = PTL))+geom_boxplot() + theme_classic()
ggplot(df,aes(HT,BWT,color = HT))+geom_boxplot() + theme_classic()
ggplot(df,aes(UI,BWT,color = UI))+geom_boxplot() + theme_classic()
ggplot(df,aes(FTV,BWT,color = FTV))+geom_boxplot() + theme_classic()


## Fitting Models
# Model building using forward selection technique and AICc score
mod.list2 <- list()
mod.list2[[1]] <- lm(BWT ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI+FTV, data=df)
mod.list2[[2]]<- lm(BWT ~ LWT+RACE+SMOKE+PTL+HT+UI, data = df)
mod.list2[[3]] <- lm(BWT ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI, data = df)
mod.list2[[4]] <- lm(BWT ~ LWT*RACE+LWT*SMOKE+LWT*PTL+LWT*HT+LWT*UI, data=df)
mod.list2[[5]] <- lm(BWT ~ AGE+RACE*HT+RACE*UI+LWT*RACE+SMOKE*PTL+SMOKE*HT+SMOKE*RACE, data=df)
mod.list2[[6]] <- lm(BWT ~ LWT+RACE+SMOKE+PTL+HT+UI+RACE*UI, data=df)

mod.aictab2=aictab(cand.set = mod.list2) 
mod.aictab2

# The model fitted by forward selection technique using AICc score is BWT ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI
# Although this model did not have the lowest AICc score, it is logical to include the mother's age as one of the predictors of the birth weight of a baby 
# because it is known that maternal age influences the birth weight of babies as seen in published articles.

## Backward selection using p-value < 0.10
mod2 <- lm(BWT ~ AGE+LWT+RACE+SMOKE+PTL+HT+UI+FTV, data=df)
summary (mod2)

# remove variables with p-value > 0.1
clean.mod2 <- update(mod2, . ~ . - FTV, data = df)
summary(clean.mod2)

## The model fitted by backward selection technique using control p-value less than 0.10 is
# BWT ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI. 

# The final model using AICc criteria is BWT ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI
final_model = lm(BWT ~ AGE + LWT + RACE + SMOKE + PTL + HT + UI, data=df)
summary(final_model)

# confidence interval
cbind(Estimate=final_model$coefficients,confint(final_model))

# The coefficients of the predictors are: 
# BWT = 2881.17 - 2.29(Age) + 4.38(LWT) - 454.40(RACE2) - 317.19(RACE3) - 331.65(SMOKE1) - 288.59(PTL1) + 1271.59(PTL3) - 576.64(HT1) - 543.57(UT1)
# The estimated values of the predictors fall within 95% confidence interval of the mean value of the population. Approximately 24 % of the variance in birth weight of newly born babies can be explained by the predictors in the model
  
#### Interpretation of the model
# A unit increase in maternal age will likely decrease the birth weight of a newborn baby by 2.29 grams, if all other variables are kept constant <br>
# A unit increase in weight of the woman at her last menstrual period will likely increase the birth weight of a newborn baby by 4.38 grams, if all other variables are kept constant <br>
# The birth weight of  a baby will likely decrease by 454.40 grams if the mother has black ancestry, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely decrease by 317.19 grams if the mother has neither black or white ancestry, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely decrease by 331.65 grams if the mother smoked cigarettes during pregnancy, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely decrease by 228.59 grams if the mother had one history of premature labor term, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely increase by 1271.89 grams if the mother had three history of premature labor term, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely decrease by 576.64 grams if the mother has a history of hypertension, while keeping all other variables constant. <br>
# The birth weight of  a baby will likely decrease by 543.57 grams if the mother has uterine irritability, while keeping all other variables constant. <br>







