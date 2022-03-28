#importing libraries

library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(tidyverse)
library(gplots)
library(stats)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
(library(biomaRt))


#Reading the file, using the transposed file where the samples are rows and probes are columns
read_expression_table <- function(filename) {
  expression_data <- readr::read_delim(filename) #read file
  transposed_data <- t(expression_data) #transpose the file, now the data is in matrix, checked with is.matrix(transposed_data)
  colnames(transposed_data)<-transposed_data[1,]#first row of matrix becomes the name of columns
  GEO_column<-colnames(expression_data) #taking column names of original file
  data<-transposed_data[-1,]%>% #removing the first row as they are column names
    as_tibble(data) %>% #converting to tibble
    add_column(GEO_column[-1], .before = "1007_s_at") #adding GEO Id as column
  exp_data <- rename(data,'subject_id'='GEO_column[-1]')#renaming the first column name alone
  return(exp_data)
}
Inten_Data<-read_expression_table("/projectnb/bf528/users/hedgehog_2022/project_1/outputs/hedgehog_intensity.csv")
Inten_numeric<-type_convert(Inten_Data)

#4. Noise filtering & dimensionality reduction
#4.1 Filter1: Expressed with a value greater than log2(15) in at least 20% of samples
logical_data<-apply(Inten_numeric>log2(15),2,sum) 
Filter1<-select_if(Inten_numeric,logical_data>0.2*134)
rownames(Filter1) <- Filter1$subject_id

#4.2 Filter2: Have a variance significantly different from the median variance of all probe sets using a threshold of  𝑝<0.01
df = nrow(Filter1)-1 #degree of freedom = n-1
chilower = qchisq((0.01)/2, df) #pvalue = 0.01, lower tail
chiupper = qchisq((1 - 0.01)/2, df, lower.tail = FALSE) #upper tail
Variance<-apply(Filter1, 2, var) #Variance of genes
Variance[is.na(Variance)] = 0
Test_statistic <- (df*Variance/median(Variance[2:39662]))
Filter2 <-select_if(Filter1, Test_statistic > chiupper) #one tailed upper test


#Filter2 for the entire data, to check how many genes are filtering
Variance1<-apply(Inten_numeric, 2, var)
Variance[is.na(Variance)] = 0
Test_statistic1 <- (133*Variance1/median(Variance1[2:54675]))
Whole_Filter <-select_if(Inten_numeric, Test_statistic1 > chiupper)

#4.3 Filter3: Have a coefficient of variation > 0.186
CV<-function(x){
  sd(x)/mean(x)
}
coeff_var<-apply(Filter2[-1],2,CV)
Filter3<-select_if(Filter2[-1], coeff_var > 0.186)
Filter3<-Filter3%>%add_column(Subjectid=rownames(Filter1), .before = '1552283_s_at')
Filter2<-Filter2%>%add_column(Subjectid=rownames(Filter1), .before = '1053_at')

#writing the files
write.csv(Filter3, "gene_expression_matrix.csv", row.names = FALSE)

write.csv(Filter2, "expression_matrix.csv", row.names = FALSE)

write.csv(Filter3, "/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/gene_expression_matrix.csv", row.names = FALSE)
write.csv(Filter2, "/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/expression_matrix.csv", row.names = FALSE)

#Number of genes passing all 3 filters

print(paste0('Number of genes passing all thresholds is ', length(Filter3)-1))


#5. Hierarchical clustering & subtype discovery

#5.1 Heirarchial clustering

Gene_expression_matrix<- readr::read_csv("/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/gene_expression_matrix.csv")
distance<-dist(Gene_expression_matrix) #distance matrix
Cluster<-hclust(distance) #clustering
#5.2 Cutting the clusters
clusters<-cutree(Cluster, k = 2)
print (paste0('No of samples in 1st cluster ', sum(clusters==1))) #1st cluster
print (paste0('No of samples in 2nd cluster ', sum(clusters==2))) #2nd cluster

#5.3 Heatmap

Y<-data.matrix(Gene_expression_matrix[-1])
rownames(Y)<-Gene_expression_matrix[[1]]
metadata<-read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
condition_colors <-
  transmute(
    metadata,
    color=if_else(SixSubtypesClassification == "C3","red","blue")
  ) #coloring by subtypes
heatmap.2(t(Y),ColSideColors=condition_colors$color, xlab='Patient Samples', ylab='Genes', 
          main='Gene Expression Across Samples',trace='none', density.info = 'none',
          key.xlab='Expression Level', scale='row', margins=c(5,5), cexRow = 0.5, cexCol = 0.2, key = 'TRUE') #heatmap
#legend(x=0.1, y=0., legend=c("C3", "C4"),fill = c('red', 'blue'))

dev.off()
heatmap.2(t(Y),ColSideColors=condition_colors$color,trace='none', density.info = 'none',scale='row',col = cm.colors(4))
#5.4 Welch T-test
Obs1<-Gene_expression_matrix[clusters==1,] #1st cluster expression values for all genes
Obs2<-Gene_expression_matrix[clusters==2,] #2nd cluster expression values for all genes
#looping through each gene and storing the t-statistic and p-value in vectors
stat<-vector('list', 1532)
p_val<-vector('list', 1532)
for (col in 2:1532) {
  x<-t.test(Obs1[,col], Obs2[,col])
  stat[col]<-x$statistic
  p_val[col]<-x$p.value
}
#p-adjust
padjust<-as.list(p.adjust(as.numeric(unlist(p_val)),method = 'fdr'))
Probeset<- colnames(Gene_expression_matrix)
diff_exp<-cbind(Probeset,stat,p_val)
Differential<-diff_exp[-1,]
Differential_expression<-cbind(Differential,padjust)
write.csv(Differential_expression, "/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/T_test.csv", row.names = FALSE)

#Number of differentially expressed genes
xyz<-abs(p.adjust(as.numeric(unlist(p_val)),method = 'fdr'))
print (paste0('Number of differentially expressed genes is ',  sum(xyz<0.05)))

#5.5  most differentially expressed genes, based on mean expression
MeanObs1<-cbind(colnames(Obs1[-1]),as_tibble(apply(Obs1[,-1],2,mean)))
MeanObs2<-cbind(colnames(Obs2[-1]),as_tibble(apply(Obs2[,-1],2,mean)))

#5.6 Expression matrix t test
expression_matrix<- readr::read_csv("/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/expression_matrix.csv")
distance_expr<-dist(expression_matrix)
Cluster_expr<-hclust(distance_expr)
clusters_expr<-cutree(Cluster_expr, k = 2)
Obs1_expr<-expression_matrix[clusters_expr==1,]
Obs2_expr<-expression_matrix[clusters_expr==2,]
stat_expr<-vector('list', length(expression_matrix))
p_val_expr<-vector('list', length(expression_matrix))
# For every gene.
for (col in 2:length(expression_matrix)) {
  y<-t.test(Obs1_expr[,col], Obs2_expr[,col])
  stat_expr[col]<-y$statistic
  p_val_expr[col]<-y$p.value
}
padjust_expr<-as.list(p.adjust(as.numeric(unlist(p_val_expr)),method = 'fdr'))
Probeset_expr<- colnames(expression_matrix)
Expr1<-cbind(Probeset_expr,stat_expr,p_val_expr)
Expr2<-Expr1[-1,]
Expr3<-cbind(Expr2,padjust_expr)
write.csv(Expr3, "/projectnb/bf528/users/hedgehog_2022/project_1/Analyst_deliverables/T_test_expr.csv", row.names = FALSE)
abc<-abs(p.adjust(as.numeric(unlist(p_val_expr)),method='fdr'))
print (paste0('Number of differentially expressed genes is ',  sum(abc<0.05)))

#write a file mapping affy ids and HGNC symbols

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_map <- as_tibble(
  getBM(
    attributes=c("affy_hg_u133_plus_2", "hgnc_symbol"),
    mart=ensembl
  )
)
#write.csv(gene_map, "gene.csv", row.names = FALSE)

