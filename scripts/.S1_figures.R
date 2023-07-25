library(ggpmisc)
library(tidyverse)
library(ggsci)
library(patchwork)
library(ggseqlogo)
library(gt)

# Import data
counts <- read.csv("sRNA_frag_counts.csv")
sample_info <- read.csv("num_reads_bam.csv")

# Clean Reads for Unlabeled Transcripts
n_reads <- sample_info$n
names(n_reads) <- sample_info$sample

counts_clean <- counts %>%
  drop_na()

# Remove the strandedness flag
counts_clean <- counts_clean %>%
  dplyr::select(!flag) %>%
  distinct()

# Get sample index & ncols
sample_index <- grep("sources", colnames(counts_clean))
num_cols <- ncol(counts)

# Group together sequences that were and were not reversed comped
counts_clean <- counts_clean %>%
  group_by(sequence) %>%
  summarise(across(1:(sample_index-2), sum)) %>%
  left_join((counts %>%
               dplyr::select("sequence","sources":as.numeric(num_cols)) %>%
               distinct()), by = "sequence")

# >20 reads
counts_filtered_plotting <- counts_clean %>%
  mutate(sums=rowSums(.[2:(sample_index - 1)])) %>%
  filter(sums > 20)

# Divide by number of loci
counts_adjusted <- counts_clean
counts_adjusted[2:(sample_index-1)] <- counts_clean[2:(sample_index-1)] / counts_clean$num

counts_filtered <- counts_adjusted %>%
  mutate(sums=rowSums(.[2:(sample_index - 1)])) %>%
  filter(sums > 20)

# LOCI BIAS #
# Try to plot with a scatter
plot_log_loci_pre_adj <- pivot_longer(counts_filtered_plotting, cols = 2:(sample_index - 1), values_to = "counts", names_to = "samples") %>%
  filter(counts < 100, counts > 0) %>%
  ggplot(aes(x = num, y = counts)) +
  stat_poly_line(formula = y~x, method = "lm") +
  stat_poly_eq(use_label("eq"), label.y = 0.2, label.x = 0.9) +
  stat_poly_eq(use_label("R2"), label.y = 0.1, label.x = .9) +
  geom_point() +
  labs(x = "Number of Loci", 
       y = "Counts",
       title = "Counts vs. Number of loci",
       subtitle = "Counts > 0")

plot_log_loci_post_adj <- pivot_longer(counts_filtered, cols = 2:(sample_index - 1), values_to = "counts", names_to = "samples") %>%
  filter(counts < 100, counts > 0) %>%
  ggplot(aes(x = num, y = counts)) +
  stat_poly_line(formula = y~x, method = "lm") +
  stat_poly_eq(use_label("eq"), label.y = 0.2, label.x = 0.9) +
  stat_poly_eq(use_label("R2"), label.y = 0.1, label.x = .9) +
  geom_point() +
  labs(x = "Number of Loci", 
       y = "Counts",
       title = "Counts vs. Number of loci Correction",
       subtitle = "Counts > 0")

# Count Density
plot_count_density <- counts_adjusted %>%
  mutate(sumsv2 = rowSums(.[2:(sample_index - 1)])) %>%
  mutate(density = sumsv2 / sum(sumsv2)) %>%
  group_by(length) %>%
  summarise(total = sum(density)) %>%
  ggplot(aes(x = length, y = total)) +
  geom_line() +
  geom_point() +
  labs(title = "Count Density for All Fragments by Length",
       y = "Proportion of Counts",
       x = "Length of Fragment")

fig1 <- (plot_log_loci_pre_adj + plot_log_loci_post_adj) / (plot_count_density) + 
  plot_annotation(title = "Loci Correction & Count Density")

ggsave("S1_correction_density.jpeg", plot = fig1, dpi = 500, width = 8, height = 7)

# Figure 1 complete #

# Motif Analysis #

max_length <- max(counts_adjusted$length)

all_pos_motifs <- list()

k <- 1
for (i in 15:(max_length)){
  all_pos_motifs[[k]] <- counts_adjusted %>%
    filter(length == i) %>%
    dplyr::mutate(sums = rowSums(.[2:(sample_index - 2)])) %>%
    dplyr::select(sums, sequence) %>%
    separate(col = sequence, sep = "", into = as.character(c(1:(max(i) + 1)))) %>%
    dplyr::select(!2) %>%
    pivot_longer(cols = c(2:(max(i) + 1)), values_to = "Base", names_to = "Position") %>%
    mutate(across(2, as.integer)) %>%
    drop_na() %>%
    arrange(Position) %>%
    group_by(Position, Base) %>%
    dplyr::count() %>%
    ungroup() %>%
    group_by(Position) %>%
    mutate(total_num = sum(n)) %>%
    mutate(prop = n/total_num) %>%
    dplyr::select(prop, Position, Base) %>%
    pivot_wider(names_from = Position, values_from = prop) %>%
    column_to_rownames("Base") %>%
    replace(is.na(.), 0) %>%
    as.matrix() %>%
    ggseqlogo() +
    ylim(0,2) + labs(title = paste("Length:", i,sep = " "),
                     subtitle = "Passed Filtering")
  
  k <- k + 1
  
}

ggsave("S1_freq_motifs.jpeg", (wrap_plots(all_pos_motifs) + plot_annotation(title = "Motifs by Frequency")), width = 21, height = 12, dpi = 500)

all_pos_motifs_counts <- list()

counts_adjusted <- counts_adjusted %>%
mutate(adjusted_sums = rowSums(.[2:(sample_index - 1)]))

k <- 1
for (i in 15:(max_length)){
  all_pos_motifs_counts[[k]]<- counts_adjusted %>%
    filter(length == i) %>%
    dplyr::mutate(sums = rowSums(.[2:(sample_index - 2)])) %>%
    dplyr::select(adjusted_sums, sequence) %>%
    separate(col = sequence, sep = "", into = as.character(c(1:(max(i) + 1)))) %>%
    dplyr::select(!2) %>%
    pivot_longer(cols = c(2:(max(i) + 1)), values_to = "Base", names_to = "Position") %>%
    group_by(Position, Base) %>%
    summarise(total = sum(adjusted_sums)) %>%
    ungroup() %>%
    mutate(across(1, as.integer)) %>%
    drop_na() %>%
    arrange(Position) %>%
    group_by(Position) %>%
    mutate(total_pos = sum(total)) %>%
    mutate(prop = total / total_pos) %>%
    dplyr::select(prop, Position, Base) %>%
    pivot_wider(names_from = Position, values_from = prop) %>%
    column_to_rownames("Base") %>%
    replace(is.na(.), 0) %>%
    as.matrix() %>%
    ggseqlogo() +
    ylim(0,2) + labs(title = paste("Length:", i,sep = " "))
  
  k <- k + 1
}

print("motif_counts")

ggsave("S1_motifs_counts.jpeg", (wrap_plots(all_pos_motifs_counts) + plot_annotation(title = "Motifs by Counts")), width = 21, height = 12, dpi = 500)

# Motif Module Done #

# Start Loci Module #
std1_loci <- counts_adjusted %>%
  filter(num > 1) %>%
  ggplot(aes(x = std)) +
  geom_density() +
  labs(title = "Distribution of Start Loci Standard Deviation")

std2_loci <- counts_adjusted %>%
  ggplot(aes(x = avg_start)) +
  geom_density() +
  labs(title = "Distribution of Start Loci")

std3_loci <- counts_adjusted %>% 
  mutate(new_bin = cut(avg_start, breaks=c(-0.01,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))) %>%
  pivot_longer(cols = 2:(sample_index - 1), values_to = "counts", names_to = "sample") %>%
  group_by(new_bin) %>%
  summarise(total = sum(counts)) %>%
  ggplot(aes(x = new_bin, y = total)) +
  geom_col() +
  labs(x = "Distance from 5' End",
       y = "Counts",
       title = "Counts vs. Fragment Distance From 5' End")

std_loci <- (std1_loci + std2_loci) / std3_loci

ggsave("S1_std_loci.jpeg", std_loci, dpi = 500, width = 10)
# std loci module done #

total_counts <- counts_adjusted %>%
  group_by(flanking5p) %>%
  count()

total_counts <- sum(total_counts$n)

counts_adjusted %>%
  group_by(flanking5p) %>%
  count() %>%
  arrange(desc(n)) %>%
  dplyr::filter(flanking5p != "") %>%
  ungroup() %>%
  slice_max(order_by = n, n = 10) %>%
  mutate(prop = n / total_counts) %>%
  gt::gt() %>%
  tab_header(title = "Top 10 5' Flanking 5-mers") %>%
  gtsave("S1_five-mer_5p.html")

total_counts_3p <- counts_adjusted %>%
  group_by(flanking3p) %>%
  count()

total_counts_3p <- sum(total_counts_3p$n)

counts_adjusted %>%
  group_by(flanking3p) %>%
  count() %>%
  arrange(desc(n)) %>%
  mutate(length = str_length(flanking3p)) %>%
  dplyr::filter(flanking3p != "", length == 5) %>%
  ungroup() %>%
  slice_max(order_by = n, n = 10) %>%
  mutate(prop = n / total_counts_3p) %>%
  dplyr::select(!length) %>%
  gt::gt() %>%
  tab_header(title = "Top 10 3' Flanking 5-mers") %>%
  gtsave("S1_five-mer_3p.html")

write.csv(counts_filtered, file = "filtered_corrected_counts.csv", row.names = F)
write.csv(counts_filtered_plotting, file = "filtered_counts.csv", row.names = F)
