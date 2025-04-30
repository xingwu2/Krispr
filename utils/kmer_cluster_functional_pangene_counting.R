rm(list = ls())

library(dplyr)
library(data.table)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input"),  type = "character", help = "dosage matrix file"),
  make_option(c("-e", "--expression"),  type = "character", help = "expression file"),
  make_option(c("-k", "--key"),  type = "character", help = "meta file"),
  make_option(c("-b", "--beta"),  type = "character", help = "beta file"),
  make_option(c("-o", "--output"), type = "character",help = "Output file")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

meta <- read.delim(opt$key,sep=",")

expression <- read.delim(opt$expression,header=F)

dm <- fread(opt$input,sep=",")
dm <- as.data.frame(dm)

beta <- read.delim(opt$beta)

cluster_names <- colnames(dm)

beta$pangene_occurance <- 0
beta$functional_pangene_count <- 0
beta$functional_pangene_ratio <- 0
beta$functional_pangene_list <- NA

for (i in 1:nrow(beta)){
  cluster_beta <- beta$kmer_effect[i]
  cluster_index <- which(beta$kmer_name[i] == cluster_names)
  non0_index <- which(dm[,cluster_index]!=0)
  pangenes <- unique(meta$Pangene_ID[non0_index])
  beta$pangene_occurance[i] <- length(pangenes)
  
  functional_pangene_count <- 0
  p_values <- c()
  pangene_beta <- c()
  for (j in 1:length(pangenes)){
    pangene_index <- which(meta$Pangene_ID==pangenes[j])
    pangene_expression <- expression$V1[pangene_index]
    cluster_i_pangene_j_dosage <- dm[pangene_index,cluster_index]
    if(var(cluster_i_pangene_j_dosage) == 0){
      p_values <- c(p_values,1)
      pangene_beta <- c(pangene_beta,0)
    } else {
      m <- summary(lm(pangene_expression~cluster_i_pangene_j_dosage))
      p_values <- c(p_values,m$coefficients[2,4])
      pangene_beta <- c(pangene_beta,m$coefficients[2,1])
    }
  }
  
  fdr_values <- p.adjust(p_values,"fdr")
  functional_pangene_index <- which(fdr_values<0.05 & cluster_beta *pangene_beta > 0)
  beta$functional_pangene_count[i] <- length(functional_pangene_index)
  beta$functional_pangene_ratio[i] <- beta$functional_pangene_count[i] / beta$pangene_occurance[i]
  beta$functional_pangene_list[i] <- paste(pangenes[functional_pangene_index],collapse=",")

}

write.table(beta,file =  paste0(opt$output,".txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
