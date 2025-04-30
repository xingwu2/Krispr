rm(list = ls())

## This script will generate dataset for 10-fold cross validation

library(optparse)
library(dplyr)
library(data.table)

option_list <- list(
  make_option(c("-i", "--input"),  type = "character", help = "Input file"),
  make_option(c("-e", "--expression"),  type = "character", help = "expression file"),
  make_option(c("-k", "--key"),  type = "character", help = "key file"),
  make_option(c("-o", "--output"), type = "character",help = "Output file")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

dm <- fread(opt$input,sep=",")
expression <- read.delim(opt$expression,header=F)
key <- read.delim(opt$key,header=T,sep=",")
  

## randomly make the 10-fold cross validation

r <- nrow(dm)

set.seed(123)
fold_id <- sample(rep(1:10, length.out=r))

for (i in 1:10) {
  test_idx  <- which(fold_id == i)
  train_idx <- setdiff(1:r, test_idx)
  
  train_dm <- dm[train_idx, ]
  test_dm  <- dm[test_idx, ]
  
  train_expression <- expression$V1[train_idx]
  test_expression <- expression$V1[test_idx]
  
  fwrite(train_dm,file = paste0(opt$output,"_DM_training_random_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  fwrite(test_dm,file = paste0(opt$output,"_DM_testing_random_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  write.table(train_expression,file = paste0(opt$output,"_expression_training_random_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
  write.table(test_expression,file = paste0(opt$output,"_expression_testing_random_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
}


## make the 10-fold cross validation by the pan-gene

#head(key)
pangene_ids <- unique(key$Pangene_ID)
l <- length(pangene_ids)

set.seed(123)
fold_id <- sample(rep(1:10, length.out=l))

for (i in 1:10) {
  test_idx  <- which(key$Pangene_ID %in% pangene_ids[which(fold_id == i)])
  train_idx <- setdiff(1:r, test_idx)
  
  train_dm <- dm[train_idx, ]
  test_dm  <- dm[test_idx, ]
  
  train_expression <- expression$V1[train_idx]
  test_expression <- expression$V1[test_idx]
  
  fwrite(train_dm,file = paste0(opt$output,"_DM_training_pangene_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  fwrite(test_dm,file = paste0(opt$output,"_DM_testing_pangene_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  write.table(train_expression,file = paste0(opt$output,"_expression_training_pangene_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
  write.table(test_expression,file = paste0(opt$output,"_expression_testing_pangene_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
}

## make the 10-fold cross validation by the individuals

#head(key)
ind_ids <- unique(key$Accession)
m <- length(ind_ids)

set.seed(123)
fold_id <- sample(rep(1:10, length.out=m))

for (i in 1:10) {
  test_idx  <- which(key$Accession %in% ind_ids[which(fold_id == i)])
  train_idx <- setdiff(1:r, test_idx)
  
  train_dm <- dm[train_idx, ]
  test_dm  <- dm[test_idx, ]
  
  train_expression <- expression$V1[train_idx]
  test_expression <- expression$V1[test_idx]
  
  fwrite(train_dm,file = paste0(opt$output,"_DM_training_accession_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  fwrite(test_dm,file = paste0(opt$output,"_DM_testing_accession_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = T)
  write.table(train_expression,file = paste0(opt$output,"_expression_training_accession_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
  write.table(test_expression,file = paste0(opt$output,"_expression_testing_accession_",i,".csv"),append = F,quote = F,sep = ",",row.names = F,col.names = F)
}

