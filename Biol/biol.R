#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133plus2.db")
BiocManager::install("GSEABase")
BiocManager::install("GSEAmining")
BiocManager::install("clusterProfiler")


# load bioconductor packages
library(tidyverse)
library(biomaRt)
library(AnnotationDbi)
library(hgu133plus2.db)
library(GSEABase)
library(GSEAmining)
library(clusterProfiler)


#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
load_expression <- function(filepath){
  data <- readr::read_delim(filepath, delim=",", col_names = TRUE)
  master_table <- tibble::tibble(data)
  return(master_table)
}

#### Mapping Gene Symbols to Probe ID ####
#' Add Gene Symbol Col to T_test_expr.csv
#'
#' @param affy A csv file with differential expression data.
#'
#' @return A tibble containing the gene symbols loaded from hgu133plus2.db. 
#' 
#' @details There are many identifier systems in hgu133plus2.db such as 
#' hgu133plus2ALIAS2PROBE, hgu133plus2ENSEMBL, hgu133plus2ENZYME, 
#' hgu133plus2GENENAME, hgu133plus2GO. The one used here is hgu133plus2ALIAS2PROBE.
#' 
#' PLEASE NOTE: 
#'   This mapping includes ALL gene symbols including those which are already listed 
#'   in the SYMBOL map. The SYMBOL map is meant to only list official gene symbols, 
#'   while the ALIAS maps are meant to store all used symbols. Therefore, the SYMBOL
#'   map is used here. 
#'
affy_to_hgnc <- function(affy) {
  data <- affy
  affy_vector <- dplyr::pull(affy, "Probeset_expr")
  gene_map <- AnnotationDbi::select(hgu133plus2.db, 
                                    keys=affy_vector, 
                                    columns=c("SYMBOL"))
  new_matrix <- merge(gene_map, data, by.x = "PROBEID", by.y = "Probeset_expr")
  new_matrix <- drop_na(new_matrix) #question here 
  return(new_matrix)
}

#### Top 10 Differentially Expressed Genes####
#'
#' @param mapped A tibble with differential expression data and hgnc_symbols.
#'
#' @return A tibble containing top 10 up and top 10 down regulated genes. 
#'
top_up_down_10 <- function(mapped) {
  top_1000 <- mapped %>%                                      
    arrange(desc(stat_expr))
  bot_1000 <- mapped %>%                                      
    arrange(stat_expr)
  top_10 <- top_1000[1:10,]
  bot_10 <- bot_1000[1:10,]
  top_10['reg'] <- rep("up")
  bot_10['reg'] <- rep("down")
  new_matrix_10 <- rbind(top_10, bot_10)
  return(new_matrix_10)
}

#### Top 1000 Up Regulated Genes####
#'
#' @param mapped A tibble with differential expression data and hgnc_symbols.
#'
#' @return A tibble containing top 1000 up regulated genes. 
#'
top_up_1000 <- function(mapped) {
  top_1000 <- mapped %>%                                      
    arrange(desc(stat_expr))
  top_1000 <- top_1000[1:1000,]
  top_1000['reg'] <- rep("up")
  return(top_1000)
}

#### Top 1000 Down Regulated Genes####
#'
#' @param mapped A tibble with differential expression data and hgnc_symbols.
#'
#' @return A tibble containing top 1000 down regulated genes. 
#'
top_down_1000 <- function(mapped) {
  bot_1000 <- mapped %>%                                      
    arrange(stat_expr)
  bot_1000 <- bot_1000[1:1000,]
  bot_1000['reg'] <- rep("down")
  return(bot_1000)
}

#### Fisher’s Exact Test on Differentially Expressed Genes and Gene Sets####
#'
#' @param affy_genes A vector with differentially expressed hgnc_symbols.
#' @param gsets_genes A vector with gene set symbols.
#' @param notreg_gens A list of hgnc_symbols not differentially expressed 
#'
#' @return A list of fisher exact parameters such as odds.ratio and p-value
#'
fishertable <- function(affy_genes, gsets_genes, notreg_genes){
  in_gset <- length(intersect(affy_genes, gsets_genes))
  not_in_gset <- length(affy_genes) - in_gset
  not_de_in_gset <- length(intersect(notreg_genes, gsets_genes))
  not_de_not_in_gset <- length(notreg_genes) - not_de_in_gset
  x <- fisher.test(matrix(c(in_gset, not_in_gset, not_de_in_gset, not_de_not_in_gset), nrow=2))
  print(x)
  return(x)
}

#### Map GMT Gene Symbols to HGNC Gene Symbols and Compute Fisher’s Exact Test####
#'
#' @param gset A vector with with gene set symbols.
#' @param report_top_1000 A vector with gene set symbols.
#' @param report_bot_1000 A list of hgnc_symbols not differentially expressed 
#' @param notup_reg A vector of genes that was not found in report_top_1000
#' @param notdown_reg A vector of genes that was not found in report_bot_1000
#'
#' @return A atomic vector data frame of gene sets and fisher exact parameters such as odds.ratio and p-value
#'
gmt2affy <- function(gset, report_top_1000, report_bot_1000, notup_reg, notdown_reg){
  data <- data.frame()
  for (i in 1:length(gset))
  {
    g <- geneIds(gset[i])
    f <- fishertable(report_top_1000$SYMBOL,g[[names(g)]],notup_reg$SYMBOL)
    fd <- fishertable(report_bot_1000$SYMBOL,g[[names(g)]],notdown_reg$SYMBOL)
    go_up <- c(pathway = names(g), p_value = f$p.value, f$estimate, exp = 'UP')
    go_down <- c(pathway = names(g), p_value = fd$p.value, fd$estimate, exp = 'DOWN')
    temp <- rbind(go_up, go_down)
    data <- rbind(data, temp)
  }
  return (data)
}

#PART 1: Load differential expression results from 5.6 and map gene symbols to probe ID
result_tib <- load_expression("Analyst code and deliverables/T_test_expr.csv")
print(nrow(result_tib))
print(result_tib)
mapped <- affy_to_hgnc(result_tib)
mapped

#PART 2: Generating top 10 and top 1000 up/down regulated gene report
report_10 <- top_up_down_10(mapped)
report_top_1000 <- top_up_1000(mapped)
report_bot_1000 <- top_down_1000(mapped)
write.csv(report_10, file = "report_top10.csv")
notup_reg <- subset(mapped, !mapped$SYMBOL %in% report_top_1000$SYMBOL) #extract genes not up regulated 
notdown_reg <- subset(mapped, !mapped$SYMBOL %in% report_bot_1000$SYMBOL) #extract genes not down regulated 

#PART 3: converting GMT to tibble - getGMT() IS DIFFICULT TO USE AND UNDERSTAND...!!!
kf <- read.gmt("Biol/c2.cp.kegg.v7.5.1.symbols.gmt")
kegg <- getGmt("Biol/c2.cp.kegg.v7.5.1.symbols.gmt",collectionType = BroadCollection(category = 'c5'), geneIdType = SymbolIdentifier())


gf <- read.gmt("Biol/c5.go.bp.v7.5.1.symbols.gmt")
go <- getGmt("Biol/c5.go.bp.v7.5.1.symbols.gmt",collectionType = BroadCollection(category = 'c5'), geneIdType = SymbolIdentifier())

hf <- read.gmt("Biol/h.all.v7.5.1.symbols.gmt")
hall <- getGmt("Biol/h.all.v7.5.1.symbols.gmt",collectionType = BroadCollection(category = 'c5'), geneIdType = SymbolIdentifier())

all <- rbind(kf, gf, hf)
all <- tibble::tibble(all)

#PART 4: Fisher Test  
#4.1: all gene collection fisher test results 
x <- fishertable(report_top_1000$SYMBOL, all$gene, notup_reg$SYMBOL)
y <- fishertable(report_bot_1000$SYMBOL, all$gene, notup_reg$SYMBOL)
allgenes_up <- c(p_value = x$p.value, x$estimate, exp = 'UP')
allgenes_down <- c(p_value = y$p.value, y$estimate, exp = 'DOWN')
allgenes_ft <- rbind(allgenes_up, allgenes_down)
allgenes_ft
write.csv(allgenes_ft, file = "all_ft.csv")

#4.2: kegg gene collection fisher test results 
x <- fishertable(report_top_1000$SYMBOL, kf$gene, notup_reg$SYMBOL)
y <- fishertable(report_bot_1000$SYMBOL, kf$gene, notup_reg$SYMBOL)
kegg_genes_up <- c(p_value = x$p.value, x$estimate, exp = 'UP')
kegg_genes_down <- c(p_value = y$p.value, y$estimate, exp = 'DOWN')
kegg_genes_ft <- rbind(kegg_genes_up, kegg_genes_down)
write.csv(kegg_genes_ft, file = "kegg_ft.csv")

#4.3: go gene collection fisher test results 
x <- fishertable(report_top_1000$SYMBOL, gf$gene, notup_reg$SYMBOL)
y <- fishertable(report_bot_1000$SYMBOL, gf$gene, notup_reg$SYMBOL)
go_genes_up <- c(p_value = x$p.value, x$estimate, exp = 'UP')
go_genes_down <- c(p_value = y$p.value, y$estimate, exp = 'DOWN')
go_genes_ft <- rbind(go_genes_up, go_genes_down)
write.csv(go_genes_ft, file = "go_ft.csv")

#4.4: hf gene collection fisher test results 
x <- fishertable(report_top_1000$SYMBOL, hf$gene, notup_reg$SYMBOL)
y <- fishertable(report_bot_1000$SYMBOL, hf$gene, notup_reg$SYMBOL)
hallmark_genes_up <- c(p_value = x$p.value, x$estimate, exp = 'UP')
hallmark_genes_down <- c(p_value = y$p.value, y$estimate, exp = 'DOWN')
hallmark_genes_ft <- rbind(hallmark_genes_up, hallmark_genes_down)
write.csv(hallmark_genes_ft, file = "hallmark_ft.csv")

#PART 5: Adjust p-values using BH-FDR and create master tibble of fisher test results 
master <- rbind(allgenes_ft, kegg_genes_ft, go_genes_ft, hallmark_genes_ft)
master<- data.frame(master)
master <- master %>%
  add_column(adj_pvalue = p.adjust(master$p_value, method = "BH", n = length(master$p_value)),.after = "p_value")
master
write.csv(master, file = "master_BH_FDR.csv")

#OR 

#PART 4: Fisher Test  
#4.1: kegg collection fisher test results 
k_data <- gmt2affy(kegg, report_top_1000, report_bot_1000, notup_reg, notdown_reg)
k_data <- data.frame(k_data)
k_data <- tibble::tibble(k_data)
k_data
k_data <- unique(k_data)
indexes <- which(k_data$p_value < 0.05)
k_data <- k_data[indexes,]
nrow(k_data) # the number of significantly enriched gene sets at adjusted 𝑝<0.05
k_data <- k_data %>%
  add_column(adj_pvalue = p.adjust(k_data$p_value, method = "BH", n = length(k_data$p_value)),.after = "p_value")
k_data
top_3_kegg <- k_data %>%                                      
  arrange(p_value)
top_3_kegg <- top_3_kegg[1:3,]
top_3_kegg
write.csv(top_3_kegg, file = "top_3_kegg.csv")

#4.2: go gene collection fisher test results 
go_data <- gmt2affy(go, report_top_1000, report_bot_1000, notup_reg, notdown_reg)
go_data <- data.frame(go_data)
go_data <- tibble::tibble(go_data)
go_data <- unique(go_data)
indexes <- which(go_data$p_value< 0.05)
go_data <- go_data[indexes,]
nrow(go_data) # the number of significantly enriched gene sets at adjusted 𝑝<0.05
go_data <- go_data %>%
  add_column(adj_pvalue = p.adjust(go_data$p_value, method = "BH", n = length(go_data$p_value)),.after = "p_value")
go_data
top_3_go <- go_data %>%                                      
  arrange(p_value)
top_3_go <- top_3_go[1:3,]
top_3_go
write.csv(top_3_go, file = "top_3_go.csv")

#4.3: hallmark gene collection fisher test results 
hall_data <- gmt2affy(hall, report_top_1000, report_bot_1000, notup_reg, notdown_reg)
hall_data <- data.frame(hall_data)
hall_data <- tibble::tibble(hall_data)
hall_data <- unique(hall_data)
indexes <- which(hall_data$p_value< 0.05)
hall_data <- hall_data[indexes,]
nrow(hall_data) # the number of significantly enriched gene sets at adjusted 𝑝<0.05
hall_data <- hall_data %>%
  add_column(adj_pvalue = p.adjust(hall_data$p_value, method = "BH", n = length(hall_data$p_value)),.after = "p_value")
hall_data
top_3_hall <- hall_data %>%                                      
  arrange(p_value)
top_3_hall <- top_3_hall[1:3,]
top_3_hall
write.csv(top_3_hall, file = "top_3_hall.csv")


