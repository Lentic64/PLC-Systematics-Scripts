
#MANOVA

#set working directory
setwd("/Volumes/Macintosh HD/Users/forrestfreund/Desktop/Berkeley/PopGen/")

#libraries needed
library(dplyr)
library(readr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(car)
library(broom)

#importing data
#Averaged leaf values
Heterozygocity <- read.csv("Heterozygosity_modified2.csv")


#Heterozygocity  plot

p <-ggdotplot(Heterozygocity, palette = "jco",
              x = "pop",
              y = "HetRatio",
              xlab = "Population",
              ylab = "Heterozygocity Ratio",
              binwidth = 0.01,
              yticks.by = 0.5,
              add = "boxplot",
              merge = TRUE,
              repel = TRUE,
              size = 6
              ) 

p + theme(axis.text.x=element_text(angle = 90, hjust = 0))

p + rotate_x_text(90)
