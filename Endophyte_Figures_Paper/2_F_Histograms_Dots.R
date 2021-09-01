# This script will combine the analysis and create final figures for my final combined data. 

getwd()                        # PLEASE CHANGE YOUR WD
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
colnames(Serendip_data)

Fig2A_above <- ggplot(Serendip_data, aes(x =Genotype, y=ShootMass, group = Condition, fill = Condition)) + 
  geom_point(aes(colour = Condition),size=3, position = position_dodge(width = .75)) + 
  scale_color_manual(values = c( "tan3", "forestgreen")) + ylab(" Shoot Mass (mg)") + 
  theme(axis.text.x = element_text(angle = 90, size = 24), axis.title.x = element_text(size = 24)) + 
  theme(axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), plot.title = element_text(size=24)) + 
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24)) + 
  ggtitle("Serendipita Shoot Mass") + theme(panel.background = element_rect(fill = "white",
                                                                            colour = "white"), 
                                            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                            colour = "light grey"), 
                                            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                            colour = "light grey")) +
  annotate("text", x = 15, y = 250, label = " * p < .05", size = 8) + coord_cartesian(xlim = c(0,12), clip = 'off')+
  annotate("text", x = 10, y = 400, label = " * ", size = 18) + 
  annotate("text", x = 11, y = 400, label = " * ", size = 18) +
  annotate("text", x = 12, y = 750, label = " * ", size = 18) + theme(legend.position = "right")
Fig2A_above


Fig2B_above <- ggplot(Serendip_data, aes(x =Genotype, y=RootMass, group = Condition, fill = Condition)) + 
  geom_point(aes(colour = Condition),size=3, position = position_dodge(width = .75)) + 
  scale_color_manual(values = c( "tan3", "forestgreen")) + ylab(" Root Mass (mg)") + 
  theme(axis.text.x = element_text(angle = 90, size = 24), axis.title.x = element_text(size = 24)) + 
  theme(axis.text.y = element_text(size = 24), axis.title.y = element_text(size = 24), plot.title = element_text(size=24)) + 
  theme(legend.title = element_text(size=24), legend.text = element_text(size=24)) + 
  ggtitle("Serendipita Root Mass") + theme(panel.background = element_rect(fill = "white",
                                                                           colour = "white"), 
                                           panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                           colour = "light grey"), 
                                           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                           colour = "light grey")) +
  annotate("text", x = 15, y = 100, label = " * p < .05", size = 8) + coord_cartesian(xlim = c(0,12), ylim = c(0,600), clip = 'off') +
  annotate("text", x = 6, y = 500, label = " * ", size = 18) + 
  annotate("text", x = 11, y = 250, label = " * ", size = 18) 
# labs(caption=str_wrap(fig2_cap, 100)) + 
# theme(plot.caption = element_text(hjust = 0.5,size = 22)) + theme(legend.position = "right")

Fig2B_above



Figure2 <- ggarrange(Fig2A_above, Fig2B_above, ncol =1, nrow = 2, labels = c("A", "B"))
Figure2


