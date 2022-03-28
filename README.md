# Project Description
### Project 1 - Hedgehog ###  
This project focuses on reproducing the results from Marisa et al., focusing on the comparisions of the C3 and C4 tumor subtypes. This project uses the combined two-phase design single dataset (134 samples in total), where an initial set of “discovery” samples was used to identify patterns among the samples, and a separate set of “validation” samples was used to test if the results from the discovery set were robust.
This repository contains files from Programmer, Analyst, and Biologist

# Contributors
### Qinrui Wu - Data Curator ###
### Dylan Beeber - Programmer. ###
### Rojashree Jayakumar - Analyst ###
### Merai Dandouch  - Biologist ###

# Repository Contents
main.R - PCA and Explained Variables 

Analyst code and deliverables/Analysis-Final.R - Noise filtering & dimensionality reduction and Hierarchical clustering & subtype discovery

Biol/biol.R - Load differential expression results from 5.6 and map gene symbols to probe ID, generates top 10 and top 1000 up/down regulated genes, converting GMT to tibble report and fishers exact test on gene sets 
