library(MMUPHin)
library(dplyr)
library(magrittr)

# Lendo a tabela de abundância e os metadados
data <- read.csv("taxonomia_meta-analise-int-agrupada.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("meta-analise-int-agrupada_metadata.tsv", sep = '\t', check.names = FALSE)

# Renomeando a coluna 'sample-id' (esse passo só é necessário se o nome tiver algum caractere que o R não gosta)
colnames(metadata)[colnames(metadata) == "sample-id"] <- "sample_id"

# Removendo colunas do 'data'
data <- data[, !(colnames(data) %in% c("study_id", "disease_group", "group"))]

# Convertendo os dados de string para número
data[] <- lapply(data, function(x) as.numeric(as.character(x)))

# Dando match entre as amostras presentes no 'data' e no 'metadata'
metadata <- metadata[match(rownames(data), metadata$sample_id), ]
rownames(metadata) <- metadata$sample_id  

# Verifica o match (se TRUE, deu certo)
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

write.csv(adjusted_table_transposed, 
          "taxonomia_meta-analise-int-agrupada_ajustada.csv", 
          quote = FALSE)