library(MMUPHin)
library(dplyr)
library(magrittr)

# Dados
data <- read.csv("../02 - processamento/amplicon/intestinal/taxonomia_meta-analise-int.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("../02 - processamento/amplicon/intestinal/meta-analise-int_metadata.tsv", sep = '\t', check.names = FALSE)

colnames(metadata)[colnames(metadata) == "sample-id"] <- "sample_id"

data <- data[, !(colnames(data) %in% c("study_id", "disease_group", "group"))]

data[] <- lapply(data, function(x) as.numeric(as.character(x)))

metadata <- metadata[match(rownames(data), metadata$sample_id), ]
rownames(metadata) <- metadata$sample_id  

all(rownames(data) == metadata$sample_id)
                 
# Ajuste                 
fit_adjust <- adjust_batch(
    feature_abd = t(data),
    batch = "study_id", # coluna que identifica os diferentes estudos
    covariates = "group", # coluna que identifica os diferentes grupos
    data = metadata
)
         
# Salvando a tabela com ajuste                 
adjusted_table <- fit_adjust$feature_abd_adj

adjusted_table_transposed <- t(adjusted_table)

write.csv(adjusted_table_transposed, "taxonomia_amplicon_corrigida.csv", quote = FALSE)