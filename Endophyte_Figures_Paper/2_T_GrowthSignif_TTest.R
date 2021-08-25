# This script will combine the analysis and create final figures for my final combined data. 

getwd()                      # PLEASE CHANGE WD
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

names(Herb_data) <- c("Genotype", "Condition", "Rep", "Endophyte","Chlorophyll1","Chlorophyll2","Chlorophyll3",
                      "PlantHeight","LeafArea","RootLength","RootVolume","Group_or_Date")
names(Burk_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                      "PlantHeight","LeafArea","RootLength","RootVolume")
names(Serendip_data) <- c("Genotype", "Condition", "Rep", "Endophyte",
                          "PlantHeight","RootLength","RootMass","ShootMass")

Serendip_data$Condition <- gsub('C', 'Control', Serendip_data$Condition)
Serendip_data$Condition <- gsub('I', 'Inoculated', Serendip_data$Condition)
Herb_data$Condition <- gsub('C', 'Control', Herb_data$Condition)
Herb_data$Condition <- gsub('I', 'Inoculated', Herb_data$Condition)
Burk_data$Condition <- gsub('C', 'Control', Burk_data$Condition)
Burk_data$Condition <- gsub('I', 'Inoculated', Burk_data$Condition)

# Average the chlorophyll together
Herb_data$Chlorophyll <- rowMeans(Herb_data[ , c("Chlorophyll1","Chlorophyll2","Chlorophyll3")], na.rm=TRUE)

##### T tests for all phenotypes
#T Tests for significance
genotypes <- unique(Herb_data$Genotype)
length(genotypes)

ttest_df <- data.frame(Experiment=character(0),
                       Trait=character(0),
                       Genotype=character(0))

#Height
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$PlantHeight~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Herb","Plant Height",as.character(unique(Herb_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Leaf Area
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$LeafArea~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Herb","Leaf Area",as.character(unique(Herb_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Length
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Herb","Root Length",as.character(unique(Herb_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Volume
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:11)]
  pvalue <- t.test(df$RootVolume~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Herb","Root Volume",as.character(unique(Herb_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}


# Chlorophyll Average
for( i in 1:length(genotypes)){
  df <-Herb_data[which (Herb_data$Genotype == (unique(Herb_data$Genotype)[i])),c(1:13)]
  pvalue <- t.test(df$Chlorophyll~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Herb_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Herb","Chlorophyll",as.character(unique(Herb_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}



##### Burkholderia
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
  if(pvalue$p.value < .1){
    adddf<-data.frame("Burk","Plant Height",as.character(unique(Burk_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Leaf Area
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$LeafArea~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Burk","Leaf Area",as.character(unique(Burk_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Length
Burk_data$RootLength=as.numeric(as.character(Burk_data$RootLength))
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Burk","Root Length",as.character(unique(Burk_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Volume
Burk_data$RootVolume=as.numeric(as.character(Burk_data$RootVolume))
for( i in 1:length(genotypes)){
  df <-Burk_data[which (Burk_data$Genotype == (unique(Burk_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootVolume~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Burk_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Burk","Root Volume",as.character(unique(Burk_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

##### Serendipita
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
  if(pvalue$p.value < .1){
    adddf<-data.frame("Serendip","Plant Height",as.character(unique(Serendip_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Length

for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootLength~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Serendip","Root Length",as.character(unique(Serendip_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Root Mass
for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$RootMass~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)  
  if(pvalue$p.value < .1){
    adddf<-data.frame("Serendip","Root Mass",as.character(unique(Serendip_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

#Shoot Mass
for( i in 1:length(genotypes)){
  df <-Serendip_data[which (Serendip_data$Genotype == (unique(Serendip_data$Genotype)[i])),c(1:8)]
  pvalue <- t.test(df$ShootMass~df$Condition,mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)
  print(unique(Serendip_data$Genotype)[i])
  print(df)
  print(pvalue)
  if(pvalue$p.value < .1){
    adddf<-data.frame("Serendip","Shoot Mass",as.character(unique(Serendip_data$Genotype)[i]))
    names(adddf)<-c("Experiment","Trait","Genotype")
    ttest_df <- rbind(ttest_df,adddf)
  }
}

# I am not sure how to parse control and inoculated values in pvalue$estimate to make a pretty method for each experiment
# So you have to run each for loop individually to identify the direction. Circle back to explore some more. 

print(" Use this table to look at scale and direction of printed t test values ")
ttest_df

write.csv(ttest_df, "Ttest_significance.csv")
