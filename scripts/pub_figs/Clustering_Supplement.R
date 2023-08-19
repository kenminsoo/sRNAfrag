library(umap)
library(tidyverse)

# Load in counts
smallRNAs <- read.delim("output_v2.tsv", sep = "\t", skip = 1)

smallRNAs <- smallRNAs[,c(1, 6:length(smallRNAs))]

total_reads <- colSums(smallRNAs[,c(3:length(smallRNAs))])

samples <- colnames(smallRNAs[3:length(smallRNAs)])

storesums <- total_reads

rows <- nrow(smallRNAs)

rpkm_transformed <- data.frame(matrix(ncol = 0, nrow = rows))
rpkm_transformed$Geneid <- smallRNAs$Geneid

for (sample in samples){
  rpkm_transformed_col <- (10^9 * smallRNAs[sample]) / (smallRNAs['Length'] * as.numeric(storesums[sample]))
  rpkm_transformed[ , ncol(rpkm_transformed) + 1] <- rpkm_transformed_col
  
}
#write.table(rpkm_transformed, "jijiwa_042523_RPKM.tsv",sep="\t",quote=F,row.names=F)

# convert the rpkm table to tpm
## tp

storeRPKMsums <- colSums(rpkm_transformed[2:length(rpkm_transformed)]) #sum of RPKM

rows <- nrow(rpkm_transformed)

tpm_transformed <- data.frame(matrix(ncol = 0, nrow = rows))
tpm_transformed$Geneid <- rpkm_transformed$Geneid

for (sample in samples){
  tpm_transformed_col <- (10^6 * rpkm_transformed[sample]) / (storeRPKMsums[sample])
  tpm_transformed[ , ncol(tpm_transformed) + 1] <- tpm_transformed_col
}

tpm_filtered <- tpm_transformed %>%
  column_to_rownames("Geneid") %>%
  mutate(total = rowSums(.)) %>%
  filter(total > 250) %>%
  dplyr::select(!total)

class_vector <- c("Neuroblastoma", "Neuroblastoma", "Neuroblastoma", "Neuroblastoma",
                   "Breast Cancer", "Breast Cancer", "Breast Cancer", "Breast Cancer",
                   "Healthy Epithelial")

cell_name_vector <- c("SK-N-SH", "SK-N-SH", "IMR32", "IMR32", "T47D", "T47D", "HCC1143", "HCC1143", "HMEC")

umap_iso <- umap(t(tpm_filtered), n_neighbors= 4)

PCA_data <- prcomp(t(tpm_filtered), scale = T)
PCA_plot <-PCA_data$x
PCA_plot <- data.frame(PCA_plot)
PCA_plot$class <- class_vector
PCA_plot$cell <- cell_name_vector
miRNA_clusters <- PCA_plot %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "miRNA")
umap_data <- data.frame(umap_iso$layout)





tpm_snoRNA <- read.csv("tpm_norm_snoRNA.csv")

umap_sno <- tpm_snoRNA %>%
  column_to_rownames("new_id") %>%
  t() %>%
  umap(n_neighbors=4)

tpm_snoRNA <- tpm_snoRNA %>%
  column_to_rownames("new_id")

PCA_data <- prcomp(t(tpm_snoRNA), scale = T)
PCA_plot <-PCA_data$x
PCA_plot <- data.frame(PCA_plot)
PCA_plot$class <- class_vector
PCA_plot$cell <- cell_name_vector
summary(PCA_data)

miRNA_clusters <- PCA_plot %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "TPM Normalization for Fragment Type + Mature miRNA",
       x = "PC1 (32.4 %)", y = "PC2 (19.5%)")

ggsave("TPM_norm.png", miRNA_clusters, height =5, width = 9, dpi = 500)

umap_data <- data.frame(umap_sno$layout)

umap_data$class <- class_vector
umap_data$cell <- cell_name_vector

snoRNA_clusters <- umap_data %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA - TPM")

ratio_snoRNA <- read.csv("ratio_norm.csv")

ratio_snoRNA <- ratio_snoRNA %>%
  dplyr::select(!c(original_id_set,outside_maps)) %>%
  column_to_rownames("new_id")

PCA_data <- prcomp(t(ratio_snoRNA), scale = T)
PCA_plot <-PCA_data$x

summary(PCA_data)


PCA_plot <- data.frame(PCA_plot)
PCA_plot$class <- class_vector
PCA_plot$cell <- cell_name_vector
miRNA_clusters <- PCA_plot %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "PCA: merged snoRNA fragment / total miRNA reads",
       x = "PC1 (32.3%)", y = "PC2 (19.5%)")

ggsave("merg_frag_ratio.jpeg",miRNA_clusters, height =5, width = 9, dpi = 500)

umap_rat <- ratio_snoRNA %>%
  t() %>%
  umap(n_neighbors=4)

umap_data <- data.frame(umap_rat$layout)

umap_data$class <- class_vector
umap_data$cell <- cell_name_vector

snoRNA_rat_clusters <- umap_data %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA Fragments with Ratio Normalization")

cpm_snoRNA <- read.csv("cpm_norm.csv")

cpm_snoRNA <- cpm_snoRNA %>%
  column_to_rownames("new_id")

PCA_data <- prcomp(t(cpm_snoRNA), scale = T)
summary(PCA_data)
PCA_plot <-PCA_data$x
PCA_plot <- data.frame(PCA_plot)
PCA_plot$class <- class_vector
PCA_plot$cell <- cell_name_vector
miRNA_clusters <- PCA_plot %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA Fragments - CPM Normalization",
       x = "PC1 (27.5%)", y = "PC2 (24.4%)")

ggsave("CPM_norm.jpeg", miRNA_clusters, height = 5, width = 9)



umap_rat <- cpm_snoRNA %>%
  t() %>%
  umap(n_neighbors=4)

umap_data <- data.frame(umap_rat$layout)

umap_data$class <- class_vector
umap_data$cell <- cell_name_vector

snoRNA_cpm_clusters <- umap_data %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA - CPM")

#filt cor count clusters

filt_corr <-  read.csv("filtered_corrected_counts.csv")
filt_corr <- filt_corr[1:10]

filt_corr <- filt_corr %>%
  column_to_rownames("sequence")

PCA_data <- prcomp(t(filt_corr), scale = T)
PCA_plot <-PCA_data$x
PCA_plot <- data.frame(PCA_plot)
PCA_plot$class <- class_vector
PCA_plot$cell <- cell_name_vector
summary(PCA_data)
miRNA_clusters <- PCA_plot %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA Fragments - Raw Adjusted Filtered Counts",
       x = "PC1 (28.2%)", y = "PC2 (22.9%)")

ggsave("raw_counts.jpeg",miRNA_clusters, height =5, width = 9, dpi = 500)


umap_corr <- filt_corr %>%
  t() %>%
  umap(n_neighbors=4)

umap_data <- data.frame(umap_corr$layout)

umap_data$class <- class_vector
umap_data$cell <- cell_name_vector

snoRNA_corr_clusters <- umap_data %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(stroke = 1.5, size = 3, aes(color = class, shape = cell)) +
  labs(title = "snoRNA - Corrected Counts")



