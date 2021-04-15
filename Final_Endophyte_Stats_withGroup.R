# This script will combine the analysis and create final figures for my final combined data. 

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

Herb_data = read.csv("Final_Herb_Data.csv", sep =",")
Burk_data = read.csv("Final_Burk_Data.csv", sep =",")
Serendip_data = read.csv("Final_Sbecsii_Data.csv", sep =",")


#get sd for error bars
data_summary <- function(data, varname, groupnames){
  require(plyr)
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
##########################################################################################
## Herbaspirillium Stuff
##########################################################################################

names(Herb_data) <- c("Genotype", "Condition", "Rep", "Endophyte","Chlorophyll1","Chlorophyll2","Chlorophyll3",
                      "PlantHeight","LeafArea","RootLength","RootVolume","Group_or_Date")


##### Visualize Phenotypes - use these two lines to make any histograms
dfsum <- data_summary(Herb_data, varname = "PlantHeight", groupnames = c("Genotype","Condition"))

ggplot(dfsum, aes(Genotype, PlantHeight, fill = Condition)) + geom_bar(stat = "summary", fun.y = "mean", position="dodge") + theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c("slateblue4", "slategrey")) + ggtitle("Herb. Inoculation and Plant Height") + geom_errorbar(aes(ymin=PlantHeight-sd, ymax=PlantHeight+sd), position=position_dodge())

##### T tests for all phenotypes
#T Tests for significance
genotypes <- unique(Herb_data$Genotype)
length(genotypes)

#Height
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$PlantHeight~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Leaf Area
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$LeafArea~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Length
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Volume
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$RootVolume~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

# Chlorophyll1
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$Chlorophyll1~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

# Chlorophyll2
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$Chlorophyll2~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

# Chlorophyll3
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$Chlorophyll3~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
}

##### Heritability                          
myvars <- as.list(colnames(Herb_data[5:11]))  # create a list of traits


# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition')), data = Herb_data) 
                 %>% anova()) 

calcs <- lapply(linreg, function(j)(j[5,2]/(j[5,2] + j[6,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
final.table <- cbind(((colnames(Herb_data[5:11]))),as.data.frame(do.call(rbind, calcs)))
names(final.table) <- c("Trait", "Heritability")

print(final.table)

# linear regression but drop rep as a variable - doesnt really change anything
linreg <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Group_or_Date + Genotype + Condition + Genotype:Condition')), data = Herb_data) 
                 %>% anova()) 

calcs <- lapply(linreg, function(j)(j[4,2]/(j[4,2] + j[5,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
final.table <- cbind(((colnames(Herb_data[5:11]))),as.data.frame(do.call(rbind, calcs)))
names(final.table) <- c("Trait", "Heritability")

print(final.table)

# non loop as an example:
#Height
height3model <- lm(Herb_data$PlantHeight ~ Herb_data$Genotype*Herb_data$Condition, data=Herb_data)
height3tab <- anova(height3model)
height3.het <- height3tab[3,2]/(height3tab[3,2] + height3tab[4,2])  #Heritability Calculation
height3.het

# Include random effects with genotype nested into group and vice versa
library(lme4)
library(lmerTest)
library(tidyverse)

# Condition + Genotype:Condition

ranef <- lapply(myvars, function(x)lmer(eval(paste0(x, '~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype) ')), data = Herb_data) 
                          %>% anova()) 

testmm <- lmer(PlantHeight ~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype), data = Herb_data, REML = FALSE)
anova(testmm)


##### Variance Decomposition Analysis - using type 1 ANOVA - loop below to handle signifigance values. 

HerbResultsT1 <- data.frame(Chlorophyll1=numeric(6),Chlorophyll2=numeric(6),Chlorophyll3=numeric(6),PlantHeight=numeric(6),LeafArea=numeric(6),RootLength=numeric(6),RootVolume=numeric(6))
rownames(HerbResultsT1) <- c("Rep","Group","Genotype","Inoculation","Genotype:Inoculation","Residuals")

HerbPlantHeightLM <- lm(PlantHeight ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(HerbPlantHeightLM)
HerbResultsT1$PlantHeight <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

HerbChl1LM <- lm(Chlorophyll1 ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(HerbChl1LM)
HerbResultsT1$Chlorophyll1 <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

linereg <- lm(Chlorophyll2 ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(linereg)
HerbResultsT1$Chlorophyll2 <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

linereg <- lm(Chlorophyll3 ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(linereg)
HerbResultsT1$Chlorophyll3 <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

linereg <- lm(LeafArea ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(linereg)
HerbResultsT1$LeafArea <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

linereg <- lm(RootLength ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(linereg)
HerbResultsT1$RootLength <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

linereg <- lm(RootVolume ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data=Herb_data)
HPH1 <- anova(linereg)
HerbResultsT1$RootVolume <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])

# Sum of Squares analysis loop
HerbResultsT1 <- data.frame(Chlorophyll1=numeric(6),Chlorophyll2=numeric(6),Chlorophyll3=numeric(6),PlantHeight=numeric(6),LeafArea=numeric(6),RootLength=numeric(6),RootVolume=numeric(6))
rownames(HerbResultsT1) <- c("Rep","Group","Genotype","Inoculation","Genotype:Inoculation","Residuals")

myvars <- names(Herb_data[5:11]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Herb_data[[m]] ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data = Herb_data)
  HPH1 <- anova(linmod)
  HerbResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])
  itlist = HPH1[,5]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
HerbT1 <- HerbResultsT1

HerbT1[] <- lapply(HerbT1[], function(x) x/sum(x))
setDT(HerbT1, keep.rownames = TRUE)[]

HerbT1_long <- HerbT1 %>%
  gather(HerbT1, value,Chlorophyll1:RootVolume)

HerbT1_long$Signif <- Signif_list

#HerbT1_long$Signif <- c("*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", "*", " ", " ", " ","*", " ", "*", " ", " ","*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", " ", " ", "*", " ")


#F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn))) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")
F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5)) + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")

F1


############## Do a variance decomposition loop with the mixed modeol no rep or residuals
HerbResultsT1 <- data.frame(Chlorophyll1=numeric(4),Chlorophyll2=numeric(4),Chlorophyll3=numeric(4),PlantHeight=numeric(4),LeafArea=numeric(4),RootLength=numeric(4),RootVolume=numeric(4))
rownames(HerbResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

# Sum of Squares analysis loop

myvars <- names(Herb_data[5:11]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Herb_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype), data = Herb_data, REML = FALSE)
  #print(summary(linmod))
  HPH1 <- anova(linmod)
  #print(HPH1)
  HerbResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,6]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
HerbT1 <- HerbResultsT1

HerbT1[] <- lapply(HerbT1[], function(x) x/sum(x))
setDT(HerbT1, keep.rownames = TRUE)[]

HerbT1_long <- HerbT1 %>%
  gather(HerbT1, value,Chlorophyll1:RootVolume)

HerbT1_long$Signif <- Signif_list

#HerbT1_long$Signif <- c("*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", "*", " ", " ", " ","*", " ", "*", " ", " ","*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", " ", " ", "*", " ")


#F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn))) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")
F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5)) + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")

F1






##### Variance Decomposition Analysis - using type 2 ANOVA
# signifigance is just the left column of the ANOVA duh you big dummy. 
# empty table
HerbResultsT2 <- data.frame(Chlorophyll1=numeric(5),Chlorophyll2=numeric(5),Chlorophyll3=numeric(5),PlantHeight=numeric(5),LeafArea=numeric(5),RootLength=numeric(5),RootVolume=numeric(5))
rownames(HerbResultsT2) <- c("Genotype","Condition","Rep","Genotype:Condition","Residuals")

#place ANOVA results into that empty table
HerbPlantHeightLM <- lm(PlantHeight ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(HerbPlantHeightLM, type = "II")
HerbResultsT2$PlantHeight=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

HerbChl1LM <- lm(Chlorophyll1 ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(HerbChl1LM, type = "II")
HerbResultsT2$Chlorophyll1=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

linereg <- lm(Chlorophyll2 ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(linereg, type = "II")
HerbResultsT2$Chlorophyll2=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

linereg <- lm(Chlorophyll3 ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(linereg, type = "II")
HerbResultsT2$Chlorophyll3=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

linereg <- lm(LeafArea ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(linereg, type = "II")
HerbResultsT2$LeafArea=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

linereg <- lm(RootLenth ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(linereg, type = "II")
HerbResultsT2$RootLength=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])

linereg <- lm(RootVolume ~ Genotype +Condition + Rep + Genotype*Condition, data=Herb_data)
HPH2 <- Anova(linereg, type = "II")
HerbResultsT2$RootVolume=c(HPH2[1,1],HPH2[2,1],HPH2[3,1],HPH2[4,1],HPH2[5,1])
HerbResultsT2

# Now Normalize all the columns
HerbT2 <- HerbResultsT2
HerbT2[] <- lapply(HerbT2[], function(x) x/sum(x))
setDT(HerbT2, keep.rownames = TRUE)[]
HerbT2_long <- HerbT2 %>%
  gather(HerbT2, value,Chlorophyll1:RootVolume)
HerbT2_long$Signif <- c("*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", "*", " ", " ", " ","*", " ", "*", " ", " ","*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", " ", " ", "*", " ")

ggplot(HerbT2_long, aes(x = HerbT2, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5)) + ggtitle("Herbaspirillum Variance Decomposition") + xlab("Phenotypes")





##########################################################################################
## Burkholderia Stuff
##########################################################################################

names(Burk_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                      "PlantHeight","LeafArea","RootLength","RootVolume")


##### Visualize Phenotypes - use these two lines to make any histograms
dfsum <- data_summary(Burk_data, varname = "PlantHeight", groupnames = c("Genotype","Condition"))

ggplot(dfsum, aes(Genotype, PlantHeight, fill = Condition)) + geom_bar(stat = "summary", fun.y = "mean", position="dodge") + theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c("slateblue4", "slategrey")) + ggtitle("Burkholderia Inoculation and Plant Height") + geom_errorbar(aes(ymin=PlantHeight-sd, ymax=PlantHeight+sd), position=position_dodge())

##### T tests for all phenotypes
#T Tests for significance
genotypes <- unique(Burk_data$Genotype)
length(genotypes)

#Height
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$PlantHeight~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Leaf Area
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$LeafArea~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Length
Burk_data$RootLength=as.numeric(as.character(Burk_data$RootLength))
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Volume
Burk_data$RootVolume=as.numeric(as.character(Burk_data$RootVolume))
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootVolume~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
}

##### Heritability                              
myvars <- as.list(colnames(Burk_data[5:8]))  # create a list of traits


# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Burk_data) 
                 %>% anova()) 

calcs <- lapply(linreg, function(j)(j[4,2]/(j[4,2] + j[5,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
final.table <- cbind(((colnames(Burk_data[5:8]))),as.data.frame(do.call(rbind, calcs)))
names(final.table) <- c("Trait", "Heritability")

print(final.table)
# non loop as an example:
#Height
heightmodel <- lm(Burk_data$PlantHeight ~ Burk_data$Genotype*Burk_data$Condition, data=Burk_data)
height3tab <- anova(heightmodel)
height3.het <- height3tab[3,2]/(height3tab[3,2] + height3tab[4,2])  #Heritability Calculation
height3.het

##### Variance Decomposition Analysis - using type 2 ANOVA





##########################################################################################
## Serendipita Stuff
##########################################################################################

names(Serendip_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                          "PlantHeight","RootLength","RootMass","ShootMass","Group_or_Date")


##### Visualize Phenotypes - use these two lines to make any histograms
dfsum <- data_summary(Serendip_data, varname = "PlantHeight", groupnames = c("Genotype","Condition"))
ggplot(dfsum, aes(Genotype, PlantHeight, fill = Condition)) + geom_bar(stat = "summary", fun.y = "mean", position="dodge") + theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c("slateblue4", "slategrey")) + ggtitle("Serendipita Inoculation and Plant Height") + geom_errorbar(aes(ymin=PlantHeight-sd, ymax=PlantHeight+sd), position=position_dodge())


# Root Mass
dfsum <- data_summary(Serendip_data, varname = "RootMass", groupnames = c("Genotype","Condition"))
dfsum['Significance'] <- c(" ", " ", " ", " "," ", " ", " ", " "," ", " ", "*", "*"," ", " ", " ", " "," ", " ", " ", " ","*", "*", " ", " ")

ggplot(dfsum, aes(Genotype, RootMass, fill = Condition)) + geom_bar(stat = "summary", fun = "mean", position="dodge") + 
  theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c( "tan3", "forestgreen")) + 
  ggtitle("Serendipita Inoculation and Root Mass") + 
  geom_errorbar(aes(ymin=RootMass-sd, ymax=RootMass+sd), position=position_dodge()) +
  geom_text(aes(label = Significance), vjust = 4.5, size = 5, position = position_dodge(width =1)) + 
  annotate("text", x = 10, y = 500, label = " * represent t-test pvalue < .05") + ylab(" Root Mass (mg)")

#Shoot Mass
dfsum <- data_summary(Serendip_data, varname = "ShootMass", groupnames = c("Genotype","Condition"))
dfsum['Significance'] <- c(" ", " ", " ", " "," ", " ", " ", " "," ", " ", " ", " "," ", " ", " ", " "," ", " ", "*", "*","*", "*", "*", "*")

ggplot(dfsum, aes(Genotype, ShootMass, fill = Condition, labels = Significance)) + geom_bar(stat = "summary", fun = "mean", position="dodge") + 
  theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c( "tan3", "forestgreen")) + 
  ggtitle("Serendipita Inoculation and Shoot Mass") + 
  geom_errorbar(aes(ymin=ShootMass-sd, ymax=ShootMass+sd), position=position_dodge()) +
  geom_text(aes(label = Significance), size = 5,vjust = 2,  position = position_dodge(width = 1)) + 
  annotate("text", x = 10, y = 1000, label = " * represent t-test pvalue < .05")+ ylab(" Shoot Mass (mg)")

################################################ Poster pics for root and shoot mass
keep_list <- c("TX303", "P39", "NC350")

dot_df <- Serendip_data[(Serendip_data$Genotype=="TX303"),]
dot_df <- rbind(Serendip_data[(Serendip_data$Genotype=="P39"),],dot_df)
dot_df <- rbind(Serendip_data[(Serendip_data$Genotype=="NC350"),],dot_df)
dot_df <- rbind(Serendip_data[(Serendip_data$Genotype=="MS71"),],dot_df)
dot_df <- rbind(Serendip_data[(Serendip_data$Genotype=="CML52"),],dot_df)
dot_df <- rbind(Serendip_data[(Serendip_data$Genotype=="B73"),],dot_df)

dot_df

dot_df$Condition <- gsub('C', 'Control', dot_df$Condition)
dot_df$Condition <- gsub('I', 'Inoculated', dot_df$Condition)

# Root Mass  
ggplot(dot_df, aes(x =Genotype, y=RootMass, group = Condition, fill = Condition)) + 
  geom_point(aes(colour = Condition),size=3, position = position_jitterdodge(jitter.width = .01)) + 
  scale_color_manual(values = c( "tan3", "forestgreen")) + ylab(" Root Mass (mg)") + theme(axis.text.x = element_text( size = 20), axis.title.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 16), plot.title = element_text(size=18)) + ggtitle("Serendipita Root Mass")


# Shoot Mass  
ggplot(dot_df, aes(x =Genotype, y=ShootMass, group = Condition, fill = Condition)) + 
  geom_point(aes(colour = Condition),size=3, position = position_jitterdodge(jitter.width = .01)) + 
  scale_color_manual(values = c( "tan3", "forestgreen")) + ylab(" Shoot Mass (mg)") + theme(axis.text.x = element_text( size = 20), axis.title.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 16), plot.title = element_text(size=18)) + ggtitle("Serendipita Shoot Mass")

##########################################################
##### T tests for all phenotypes
#T Tests for significance
genotypes <- unique(Serendip_data$Genotype)
length(genotypes)

#Height
for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$PlantHeight~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Length

for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Root Mass
for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootMass~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
}

#Shoot Mass
for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$ShootMass~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
}

##### Heritability  - Test with other ANOVA's       
myvars <- as.list(colnames(Serendip_data[5:8]))  # create a list of traits
# Type 1
# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition')), data = Serendip_data) 
                 %>% anova()) 
calcs <- lapply(linreg, function(j)(j[5,2]/(j[5,2] + j[6,2])))#extract and calculate the actual heritability for each table
#make a final table with the traits and calculated heritability
final.table <- cbind(((colnames(Serendip_data[5:8]))),as.data.frame(do.call(rbind, calcs)))
names(final.table) <- c("Trait", "Heritability")

print(final.table)

# Type II
linreg2 <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Serendip_data) 
                 %>% Anova(type="II")) 
lapply(linreg2, function(j)print(j))
calcs2 <- lapply(linreg2, function(j)(j[4,1]/(j[4,1] + j[5,1])))#extract and calculate the actual heritability for each table
#make a final table with the traits and calculated heritability
final.table2 <- cbind(((colnames(Serendip_data[5:8]))),as.data.frame(do.call(rbind, calcs2)))
names(final.table2) <- c("Trait", "Heritability")

print(final.table2)

# Type III
linreg3 <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Serendip_data) 
                 %>% Anova(type="III")) 
calcs3 <- lapply(linreg3, function(j)(j[5,1]/(j[5,1] + j[6,1])))#extract and calculate the actual heritability for each table
#make a final table with the traits and calculated heritability
final.table3 <- cbind(((colnames(Serendip_data[5:8]))),as.data.frame(do.call(rbind, calcs3)))
names(final.table3) <- c("Trait", "Heritability")

print(final.table3)
##### Variance decomposition analysis - ANOVA type I
# Sum of Squares analysis loop
SereResultsT1 <- data.frame(PlantHeight=numeric(6),RootLength=numeric(6),RootMass=numeric(6),ShootMass=numeric(6))
rownames(SereResultsT1) <- c("Rep", "Group","Genotype","Inoculation","Genotype:Inoculation", "Residuals")

myvars <- names(Serendip_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lm(Serendip_data[[m]] ~ Rep + Group_or_Date + Genotype + Condition + Genotype:Condition, data = Serendip_data)
  HPH1 <- anova(linmod)
  SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2],HPH1[6,2])
  itlist = HPH1[,5]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
HerbT1 <- SereResultsT1

HerbT1[] <- lapply(HerbT1[], function(x) x/sum(x))
setDT(HerbT1, keep.rownames = TRUE)[]

HerbT1_long <- HerbT1 %>%
  gather(HerbT1, value,PlantHeight:ShootMass)

HerbT1_long$Signif <- Signif_list

#HerbT1_long$Signif <- c("*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", "*", " ", " ", " ","*", " ", "*", " ", " ","*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", " ", " ", "*", " ")


#F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn))) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")
F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5)) + ggtitle("Variance Decomposition Analysis: Serendipita") + xlab("Growth Promoted Phenotypes")

F1





# Sum of Squares analysis loop - mixed models with no rep and no residuals include date
SereResultsT1 <- data.frame(PlantHeight=numeric(4),RootLength=numeric(4),RootMass=numeric(4),ShootMass=numeric(4))
rownames(SereResultsT1) <- c("Group","Genotype","Inoculation","Genotype:Inoculation")

myvars <- names(Serendip_data[5:8]) # create a list of traits
Signif_list <- list()

for( m in myvars){
  print(m)
  print(as.name(m))
  linmod <- lmer(Serendip_data[[m]] ~ Group_or_Date + Genotype + Condition + Genotype:Condition + (1 | Group_or_Date/Genotype), data = Serendip_data, REML = FALSE)
  HPH1 <- anova(linmod)
  print(HPH1)
  SereResultsT1[[m]] <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2])
  itlist = HPH1[,6]
  Signif_list <- append(Signif_list, itlist)
}

#Get rid of the NAs from residuals
Signif_list[is.na(Signif_list)] = 1
Signif_list

# Convert these p values to *
for( s in 1:length(Signif_list)){
  #print(s)
  if (Signif_list[s] < .05){
    Signif_list[s] <- '*'
    #print("True")
  } else {
    Signif_list[s] <- " "
    #print("False")
  }
}
Signif_list

# Now Normalize all the columns
HerbT1 <- SereResultsT1

HerbT1[] <- lapply(HerbT1[], function(x) x/sum(x))
setDT(HerbT1, keep.rownames = TRUE)[]

HerbT1_long <- HerbT1 %>%
  gather(HerbT1, value,PlantHeight:ShootMass)

HerbT1_long$Signif <- Signif_list

#HerbT1_long$Signif <- c("*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", "*", " ", " ", " ","*", " ", "*", " ", " ","*", " ", " ", " ", " ","*", " ", " ", " ", " ","*", " ", " ", "*", " ")


#F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn))) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + ggtitle("Variance Decomposition Analysis: Herbaspirillum") + xlab("Growth Promoted Phenotypes")
F1 <- ggplot(HerbT1_long, aes(x = HerbT1, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5)) + ggtitle("Variance Decomposition Analysis: Serendipita") + xlab("Growth Promoted Phenotypes")

F1



############ Non looped way without date
SereResultsT1 <- data.frame(PlantHeight=numeric(5),RootLength=numeric(5),RootMass=numeric(5),ShootMass=numeric(5))
rownames(SereResultsT1) <- c("Genotype","Inoculation","Rep","Genotype:Inoculation","Residuals")

SerePlantHeightLM <- lm(PlantHeight ~ Genotype +Condition + Rep + Genotype*Condition, data=Serendip_data)
HPH1 <- anova(SerePlantHeightLM)
SereResultsT1$PlantHeight <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2])

linereg <- lm(RootLength ~ Genotype +Condition + Rep + Genotype*Condition, data=Serendip_data)
HPH1 <- anova(linereg)
SereResultsT1$RootLength <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2])

linereg <- lm(RootMass ~ Genotype +Condition + Rep + Genotype*Condition, data=Serendip_data)
HPH1 <- anova(linereg)
SereResultsT1$RootMass <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2])

linereg <- lm(ShootMass ~ Genotype +Condition + Rep + Genotype*Condition, data=Serendip_data)
HPH1 <- anova(linereg)
SereResultsT1$ShootMass <- c(HPH1[1,2],HPH1[2,2],HPH1[3,2],HPH1[4,2],HPH1[5,2])

# Now Normalize all the columns
SereT1 <- SereResultsT1

SereT1[] <- lapply(SereT1[], function(x) x/sum(x))

setDT(SereT1, keep.rownames = TRUE)[]

SereT1_long <- SereT1 %>%
  gather(SereT1, value,PlantHeight:ShootMass)

SereT1_long$Signif <- c("*", "*", " ", " "," ", "*", "*", " "," ", " ", "*", "*"," ", " ", " ", "*","*", " ", " ", " ")

F1C <- ggplot(SereT1_long, aes(x = SereT1, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text(angle = 90)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5))+ ggtitle("Variance Decomposition Analysis: Serendipita") + xlab("Growth Promoted Phenotypes")

#################### Poster figure: combine herb and serendip variance

SereT1_long["Endophyte"] <- "Serendipita"
HerbT1_long["Endophyte"] <- "Herbaspirillum"

SereT1_long$SereT1 <- paste("Serendipita", SereT1_long$SereT1, sep="_")
HerbT1_long$HerbT1 <- paste("Herbaspirillum", HerbT1_long$HerbT1, sep="_")

names(SereT1_long)[names(SereT1_long) == "SereT1"] <- "Phenotype"
names(HerbT1_long)[names(HerbT1_long) == "HerbT1"] <- "Phenotype"

SandH <- rbind(SereT1_long, HerbT1_long)
ggplot(SandH, aes(x = Phenotype, y = value, fill = forcats::fct_rev(rn), label = Signif)) + geom_col(position=position_stack()) + theme(axis.text.x = element_text( size = 12)) + theme(axis.text.y = element_text(size = 14)) + labs(fill = "Variables") + geom_text(aes(label = Signif), size = 5, position = position_stack(vjust = .5))+ ggtitle("Variance Decomposition Analysis") + xlab("Growth Promoted Phenotypes") + coord_flip() 

###################################################################################################
# Real Time - qPCR    
###################################################################################################
# Note: each sample had 2 technical reps! Experiment 3 had a different dna concentration due to differences in 
# extracted DNA concnetrations. (they were too low compared to exps 1 and 2)

qpcr_data <- read.csv("CombinedRaw_qPCR_clean.csv", header = TRUE, sep = ",")
# Toss any rows with missing values
qpcr_data = qpcr_data[complete.cases(qpcr_data),]

#add empty column
qpcr_data['calcs'] <- NA

# delta ct: experimental gene - housekeeping gene - this is for a single sample
qpcr_delta1 <- transform(qpcr_data, calcs = Cp_ITS - Cp_CDK)

# average all the delta ct values by genotype and condition so you can do the second delta from the combined samples
qpcr_sum <- summarySE(data = qpcr_delta1, measurevar = "calcs", groupvars = c("Genotype","Condition"), na.rm = TRUE)

qgenotypes <- unique(qpcr_sum$Genotype)
length(qgenotypes)

# find the double delta for each genotype
doubledelta_df <- data.frame("Genotype" = qgenotypes, "Double_Delta" = NA)

# Inoculated - Control 
for( i in 1:length(qgenotypes)){
  df <-qpcr_sum[which (qpcr_sum$Genotype == (unique(qpcr_sum$Genotype)[i])),c(1:7)]
  doubled <- df[2,4]-df[1,4]
  doubledelta_df[i,2] = doubled
}

#Add experiments so we can fill by them
doubledelta_df['Experiment_Number'] <- c('3','3','3','1','3','2','1','2','1','1','1','3')

# Drop CML69 and MO18W - we have no phenotype data due to greenhouse fertilization by a worker
doubledelta_df <- doubledelta_df[!(doubledelta_df$Genotype=="CML69"),]
doubledelta_df <- doubledelta_df[!(doubledelta_df$Genotype=="MO18W"),]

# Graph double delta
ggplot(doubledelta_df, aes(Genotype, Double_Delta, fill = Experiment_Number)) + geom_bar(stat = "summary", fun = "mean", position="dodge") + theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c("slateblue4", "slategrey", "indianred3")) + ggtitle("Average Double Delta Value of Serendipita Inoculation by Genotype")

# Transform to fold change
foldchange_df = transform(doubledelta_df, Double_Delta = log2(-Double_Delta))

# Add phenotype significance from t - tests .1 - .05
foldchange_df['Significance'] <- c('*',' ',' ',' ','***','*',' ','***','***','***')

#Graph fold change
#ggplot(foldchange_df, aes(Genotype, Double_Delta, fill = Experiment_Number)) + geom_bar(stat = "summary", fun = "mean", position="dodge") + theme(axis.text.x = element_text(size=12)) + scale_fill_manual(values = c("slateblue4", "slategrey", "indianred3")) + ggtitle("Expression Change of Serendipita Inoculation by Genotype")

ggplot(foldchange_df, aes(x = reorder(Genotype, -Double_Delta), Double_Delta, fill = Experiment_Number, label = Significance)) + geom_bar(stat = "summary", fun = "mean", position="dodge") + theme(plot.title = element_text(size=18),axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), axis.title = element_text(size = 18)) + scale_fill_manual(values = c("dodgerblue", "orange")) + ggtitle("Expression Change of Serendipita Inoculation by Genotype") + xlab("Genotype") + ylab("2^Double_Delta") + geom_text(aes(label = Significance), size = 5, position = position_stack(vjust = .5)) + annotate("text", x = 8, y = 3.5, size=5, label = " P-Values for a phenotype less than:\n * = 0.1,    *** = 0.05 ") 

# T Test comparing inoculated and control single delta values. 
genotypes <- unique(qpcr_delta1$Genotype)
length(genotypes)

for( i in 1:length(genotypes)){
  df2 <-qpcr_delta1[which (qpcr_delta1$Genotype == (unique(qpcr_delta1$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df2$calcs~df2$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(qpcr_delta1$Genotype)[i])
  print(df2)
  print(pvalue)
}

#Heritability
qpcrFinModel <- lm(qpcr_delta1$calcs~ qpcr_delta1$Genotype*qpcr_delta1$Condition, data=qpcr_delta1)
qPCRFintab <- anova(qpcrFinModel)
colonizationFinal.het <- qPCRFintab[3,2]/(qPCRFintab[3,2] + qPCRFintab[4,2])
