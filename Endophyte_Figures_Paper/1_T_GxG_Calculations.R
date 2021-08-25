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


##### Heritability                          
herb_cols <- c(13, 8:11)
myvars <- names(Herb_data[,herb_cols]) # create a list of traits


# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Herb_data) 
                 %>% anova()) 

calcs <- lapply(linreg, function(j)(j[4,2]/(j[4,2] + j[5,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
end <- data.frame(endo=c("Herb","Herb","Herb","Herb","Herb"))

final.table1 <- cbind((names(Herb_data[,herb_cols])),as.data.frame(do.call(rbind, calcs)),end)

names(final.table1) <- c("Trait", "Heritability","Endo")

print(final.table1)


##### Serendipita                          
myvars <- as.list(colnames(Serendip_data[5:8]))  # create a list of traits


# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg2 <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Serendip_data) 
                 %>% anova()) 

calcs2 <- lapply(linreg2, function(j)(j[4,2]/(j[4,2] + j[5,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
end <- data.frame(endo=c("Serend","Serend","Serend","Serend"))

final.table2 <- cbind(((colnames(Serendip_data[5:8]))),as.data.frame(do.call(rbind, calcs2)),end)
names(final.table2) <- c("Trait", "Heritability","Endo")

print(final.table2)


### Burkholderia
# keeps defaulting to cols of factors, need to fix it or ANOVA wont works
Burk_data[,5:8] <- sapply(Burk_data[,5:8], as.numeric)

myvars <- as.list(colnames(Burk_data[5:8]))  # create a list of traits

# use lapply to loop through columns and create linear models. evale paste0 makes the variable work somehow
linreg3 <- lapply(myvars, function(x)lm(eval(paste0(x, '~ Rep + Genotype + Condition + Genotype:Condition')), data = Burk_data) 
                  %>% anova()) 

calcs3 <- lapply(linreg3, function(j)(j[4,2]/(j[4,2] + j[5,2])))#extract and calculate the actual heritability for each table

#make a final table with the traits and calculated heritability
end <- data.frame(endo=c("Burk","Burk","Burk","Burk"))

final.table3 <- cbind(((colnames(Burk_data[5:8]))),as.data.frame(do.call(rbind, calcs3)),end)
names(final.table3) <- c("Trait", "Heritability","Endo")

print(final.table3)

herit_table <- do.call("rbind", list(final.table1,final.table2,final.table3))
print(herit_table)


write_csv(herit_table,"GxG_Table.csv")
