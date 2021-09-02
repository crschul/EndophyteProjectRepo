# Testing Mixed Models for confounding variables 

getwd()
setwd("/home/coreyschultz/1.Projects/2.Maize.Endophyte.Project/Final_Data_and_Figs")
getwd()
library(tidyr)
library(tidyverse)
library(broom)
library(ggplot2)
library(car)
library(data.table)
library(ggpubr)
library(lme4)
library(lmerTest)

Herb_data = read.csv("Final_Herb_Data.csv", sep =",")
Burk_data = read.csv("Final_Burk_Data.csv", sep =",")
Serendip_data = read.csv("Final_Sbecsii_Data.csv", sep =",")

############## Do a variance decomposition loop with the mixed modeol no rep or residuals
HerbResultsT1 <- data.frame(Chlorophyll1=numeric(4),Chlorophyll2=numeric(4),Chlorophyll3=numeric(4),PlantHeight=numeric(4),LeafArea=numeric(4),RootLength=numeric(4),RootVolume=numeric(4))
rownames(HerbResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

# Sum of Squares analysis loop

myvars <- names(Herb_data[5:11]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  #linmod <- lmer(Herb_data[[m]] ~ (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Herb_data, REML = TRUE,
   #             control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  # linmod <- lmer(Herb_data[[m]] ~ (1 | Group_or_Date), data = Herb_data, REML = TRUE,
  #                control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  # linmod <- lmer(Herb_data[[m]] ~ (1 | Genotype) + Condition + Genotype:Condition, data = Herb_data, REML = TRUE,
  #                control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  linmod <- lmer(Herb_data[[m]] ~ (1 | Genotype) + (1 | Group_or_Date), data = Herb_data, REML = TRUE,
                 control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  rtab <- ranova(linmod)
  print(rtab)
  relgrad <- with(linmod@optinfo$derivs,solve(Hessian,gradient))
  print(max(abs(relgrad)))
  #print(summary(linmod))
  print("####################################################")
}

for( m in myvars){
  print(m)
  print(as.name(m))
  #linmod <- lmer(Herb_data[[m]] ~  Genotype + Condition + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Herb_data, REML = FALSE)
  linmod <- lmer(Herb_data[[m]] ~  (1 | Genotype) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Herb_data, REML = FALSE)
  #HPH1 <- anova(linmod, type = c("I"))
  #print(HPH1) #anova
  rtab <- ranova(linmod)
  print(rtab)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  print("####################################################")
  #print(summary(linmod)) # then print mm summary for residuals
}

# convergence test: https://github.com/lme4/lme4/issues/120

SereResultsT1 <- data.frame(PlantHeight=numeric(4),RootLength=numeric(4),RootMass=numeric(4),ShootMass=numeric(4))
rownames(SereResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

myvars <- names(Serendip_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  #linmod <- lmer(Serendip_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype), data = Serendip_data, REML = FALSE)
  #linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = TRUE)
  linmod <- lmer(Serendip_data[[m]] ~  (1 | Genotype) + (1 | Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  HPH1 <- anova(linmod, type = c("I"))
  print(HPH1) #anova
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  #print(summary(linmod)) # then print mm summary for residuals
  # SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  # itlist = HPH1[,6]
  # Signif_list <- append(Signif_list, itlist)
}

# All random
for( m in myvars){
  print(m)
  print(as.name(m))
  #linmod <- lmer(Serendip_data[[m]] ~  (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  linmod <- lmer(Serendip_data[[m]] ~  (1 | Group_or_Date) + (1 | Genotype), data = Serendip_data, REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=200000) ))
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
  rtab <- ranova(linmod)
  print(rtab)
  #print(summary(rtab))
  #print("/n")
  #print(summary(linmod))
  #print(deviance(linmod))
  print("################################################################")
}

# All random - try bayesian ANOVA
library(BayesFactor)

Serendip_data$Group_or_Date <- as.factor(Serendip_data$Group_or_Date)
banov <- anovaBF(PlantHeight ~ Genotype + Condition + Group_or_Date + Genotype:Condition, 
                 data = Serendip_data, whichRandom = "Group_or_Date", progress = FALSE, iterations = 20000)
print(banov)

# Anova summary no model just to see
summary(anova(PlantHeight ~ Genotype + Condition + Genotype:Condition + Error(Group_or_Date), 
            data = Serendip_data, type=("I")))

# Only group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  print(varcomp, comp = "Variance")
  print("/n")
  print(summary(linmod))
  print(deviance(linmod))
  print("################################################################")
  
}

# Anova just group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  HPH1 <- anova(linmod, type=c("I"))
  print(HPH1)
  print(summary(linmod))
}

# Anova plain lm
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ (1 |Group_or_Date) + Genotype + Condition + Genotype:Condition, data = Serendip_data)
  HPH1 <- anova(linmod)
  print(HPH1)
}

# Interaction +  group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Genotype + Condition + (1 |Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}

# Interaction + genotype + group random
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ (1 |Genotype) + Condition + (1 |Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}

# Just Genotype and group
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ (1 |Genotype) + Condition + Genotype:Condition + (1 | Group_or_Date), data = Serendip_data, REML = FALSE)
  varcomp <- VarCorr(linmod)
  print(varcomp, comp = "Variance")
}


library(pander)

linmod <- lmer(PlantHeight ~  (1 | Genotype) + (1 | Condition) + (1 | Genotype:Condition) + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("I"))
pander(HPH1)
varcomp <- VarCorr(linmod)
print(varcomp, comp = "Variance")
get_variance(linmod)


linmod <- lmer(PlantHeight ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("I"))
pander(HPH1)

linmod <- lmer(PlantHeight ~ Genotype + Condition + Genotype:Condition + Group_or_Date + (1 | Group_or_Date), data = Serendip_data)
HPH1 <- anova(linmod, type = c("III"))
pander(HPH1)

print(aov(PlantHeight ~ Genotype + Condition + Genotype:Condition + Error(Group_or_Date), data = Serendip_data))

library(nlme)
m1 <- lme(PlantHeight ~ Genotype + Condition + Genotype:Condition, random=~1|Group_or_Date, data = Serendip_data)
anova(m1)
summary(m1)


############################################# Test to see how confounding group vs. genotype is
# linear models

myvars <- names(Herb_data[5:11]) # create a list of traits

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Group_or_Date + Genotype:Condition, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  #rtab <- ranova(linmod)
  #print(rtab)
  #varcomp <- VarCorr(linmod)
  #print(varcomp, comp = "Variance")
  print("####################################################")
  #print(summary(linmod)) # then print mm summary for residuals
}


for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Group_or_Date*Genotype, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  #rtab <- ranova(linmod)
  #print(rtab)
  #varcomp <- VarCorr(linmod)
  #print(varcomp, comp = "Variance")
  print("####################################################")
}


################## Test that group is a counfounding variable

# ANCOVA

myvars <- names(Herb_data[5:11])

model1lm = lm(PlantHeight ~ Genotype + Group_or_Date, data = Herb_data)

model1 = aov(PlantHeight ~ Genotype + Group_or_Date, data = Herb_data)
Anova(model1, type="III") # throws an error: there is perfect multicollinearity
alias(PlantHeight ~ Genotype + Group_or_Date, data = Herb_data) 
# This shows that genotype terms are linearly dependent upon Group or date. 

# search for large variance inflation factors 
car::vif(model1lm)
# unable to run do to aliased coefficients in the model. 

summary.lm(model1)
posth=glht(model1, linfct=mcp(factorvariable="Tukey"))  ##gives the post-hoc Tukey analysis
summary(posth)

#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/alias
#https://stats.stackexchange.com/questions/51780/how-to-perform-an-ancova-in-r


# Loops to find dependencies in a Model
myvars <- names(Herb_data[5:11])
for( m in myvars){
  print(m)
  print(as.name(m))
  print(alias(Herb_data[[m]] ~ Genotype + Group_or_Date, data = Herb_data)) 
  print("####################################################")
}

myvars <- names(Serendip_data[5:8]) 
for( m in myvars){
  print(m)
  print(as.name(m))
  print(alias(Serendip_data[[m]] ~ Genotype + Group_or_Date, data = Serendip_data)) 
  print("####################################################")
}

# Loops for type I ANOVA test how switching variable order?

myvars <- names(Herb_data[5:11])
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Group_or_Date + Genotype + Group_or_Date:Genotype, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  print("####################################################")
}

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Genotype + Group_or_Date + Genotype:Group_or_Date, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  print("####################################################")
}

# Include Condition and Interaction. 
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Genotype + Group_or_Date + Condition + Genotype:Condition, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  print("####################################################")
}

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Condition, data = Herb_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  print("####################################################")
}

myvars <- names(Serendip_data[5:8]) 
SereResultsT1 <- data.frame(PlantHeight=numeric(3),RootLength=numeric(3),RootMass=numeric(3),ShootMass=numeric(3))
rownames(SereResultsT1) <- c("Group","Genotype","Combined")

dfcol = 1
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ Group_or_Date + Genotype + Group_or_Date:Genotype, data = Serendip_data)  
  HPH1 <- anova(linmod)
  #SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2])
  print(HPH1) #anova
  SereResultsT1[1,dfcol] <- HPH1[1,2]
  SereResultsT1[2,dfcol] <- HPH1[2,2]
  dfcol = dfcol + 1
  print("####################################################")
}

dfcol = 1
for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ Genotype + Group_or_Date + Genotype:Group_or_Date, data = Serendip_data)  
  HPH1 <- anova(linmod)
  print(HPH1) #anova
  SereResultsT1[3,dfcol] <- HPH1[1,2]
  dfcol = dfcol + 1
  print("####################################################")
}

SereResultsT1



