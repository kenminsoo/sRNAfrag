library(tidyverse)
library(patchwork)
library(ggpmisc)

# A more clean way of annotating our transcripts
mature <- read.csv("mature.csv")
pre <- read.csv("pre.csv")

joined <- mature %>%
  left_join(pre, by = c("id" = "name"), suffix = c(".mature",".pre"))

joined_plus <- joined %>%
  filter(strand.mature == "+")

joined_minus <- joined %>%
  filter(strand.mature == "-")

joined_plus$start <- joined_plus$start.mature - joined_plus$start.pre + 1
joined_plus$end <- joined_plus$end.mature - joined_plus$start.pre + 1

joined_minus$start <- joined_minus$end.pre - joined_minus$end.mature + 1
joined_minus$end <- joined_minus$end.pre - joined_minus$start.mature + 1

joined_long_plus <- joined_plus %>%
  pivot_longer(cols = c(start,end), values_to = "loci", names_to = "type")

joined_long_minus <- joined_minus %>%
  pivot_longer(cols = c(start,end), values_to = "loci", names_to = "type")

joined_small_plus <- joined_long_plus %>%
  select(derived_id, type, loci)

joined_small_minus <- joined_long_minus %>%
  select(derived_id, type, loci)

all_loci <- rbind(joined_small_plus, joined_small_minus)

for (i in -6:6){
  temp_col <- (all_loci$loci + i) %>%
    as.data.frame()
  
  colnames(temp_col) <- i
  
  all_loci <- cbind(all_loci, temp_col)
}

# build offsets to -5, +5
even_longer_df <- all_loci %>%
  pivot_longer(cols = c(4:16), values_to = "location", names_to = "offset")

even_longer_df <- even_longer_df %>%
  select(-loci)

# Import peaks
pot_peaks <- read.csv("cluster_peak_relationship_table.csv")

# Select for primary clusters & find expected counts
pot_peaks_filtered <- pot_peaks %>%
  filter(primary_cluster_peak == "True")

labeled_hits <- pot_peaks_filtered %>%
  left_join(even_longer_df, by = c("loci" = "location", "type", "source" = "derived_id"))

results_flai <- read.delim("results_flaimapper.tsv", sep = "\t")

# WE now filter the flai mapper for reads >10 
results_flai_filtered <- results_flai %>%
  filter(Corresponding.reads..total. > 10)

# Now we want to make the dataframe longer
results_flai_filtered <- results_flai_filtered %>%
  pivot_longer(cols = c("Start", "End"), names_to = "type", values_to = "loci")

results_flai_filtered$loci <- results_flai_filtered$loci + 1

results_flai_filtered$type <- str_to_lower(results_flai_filtered$type)

labeled_hits_flai <- results_flai_filtered %>%
  left_join(even_longer_df, by = c("loci" = "location", "type", "Reference.sequence" = "derived_id"))

labeled_hits %>%
  drop_na()

labeled_hits_flai %>%
  drop_na()

counts_df <- read.delim("output_v2.tsv", sep = "\t", skip = 1)
counts_df_lab <- read.delim("output.tsv", sep = "\t", skip = 1)

counts_df<- counts_df %>%
  mutate(sums = rowSums(.[,7:15]))

counts_df <- counts_df %>%
  filter(sums > 0)

counts_df_lab <- counts_df_lab %>%
  separate(col = Geneid, into = c("Derived", "ID"), sep = "_")

counts_df_lab <- counts_df_lab %>%
  select(Derived, ID)

counts_df <- counts_df %>%
  left_join(counts_df_lab, by = c("Geneid" = "ID"))

labler_df <- pre %>%
  select(name, derived_id) %>%
  distinct()

counts_df_labeled <- counts_df %>%
  left_join(labler_df, by = c("Derived" = "name"))

counts_sumar <- counts_df_labeled %>%
  select(sums, derived_id) %>%
  group_by(derived_id) %>%
  summarise(total=sum(sums))

counts_sumar_exp <- pot_peaks %>%
  select(counts,source) %>%
  group_by(source) %>%
  summarise(total=sum(counts))

counts_metrics <- counts_sumar_exp %>%
  left_join(counts_sumar, by = c("source" = "derived_id"), suffix = c("exp", "true")) %>%
  ggplot(aes(x = totaltrue, y = totalexp)) +
  geom_point() +stat_poly_line(formula = y~x, method = "lm") +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(use_label("R2"), label.y = 0.8) + 
  labs(x = "True Counts", y = "Peal Calling Counts",
       title = "Peak Calling vs. True Counts",
       caption = "Counts generated with AASRA")
  
# OFFSET COUNTS
joined_for_metrics <- pot_peaks_filtered %>%
  inner_join(even_longer_df, by = c("source" = "derived_id", "type", "loci" = "location"))

joined_for_metrics_flai <- results_flai_filtered %>%
  drop_na() %>%
  inner_join(even_longer_df, by = c("Reference.sequence" = "derived_id", "type", "loci" = "location"))

found_mirnas <- pot_peaks_filtered$source %>%
  unique() %>%
  as.data.frame()

found_mirnas_flai <- joined_for_metrics_flai$Reference.sequence %>%
  unique() %>%
  as.data.frame()

colnames(found_mirnas) <- "source"

all <- even_longer_df %>%
  filter(offset == 0) %>%
  inner_join(found_mirnas, by = c("derived_id" = "source"))

combined_metrics <- joined_for_metrics %>%
  group_by(type, offset) %>%
  count() %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(ratio = n / sum(n)) %>%
  mutate(pipeline = "sRNAfrag")

combined_metrics <- rbind(combined_metrics, joined_for_metrics_flai %>%
                            group_by(type, offset) %>%
                            count() %>%
                            ungroup() %>%
                            group_by(type) %>%
                            mutate(ratio = n / sum(n)) %>%
                            mutate(pipeline = "Flaimapper"))

offset_graph <- combined_metrics %>%
  ggplot(aes(x = as.integer(offset), y = ratio, fill = pipeline)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(data = . %>% filter(pipeline == "sRNAfrag"),aes(label=round(ratio,2), y = ratio + 0.015), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x = "Offset From Truth Position (bp)",
       y = "Proportion",
       title = "miRNA Predicted Peaks",
       subtitle = "Text representing sRNAfrag proportion") +
  facet_grid(~factor(type, levels = c("start","end"))) +
  ylim(0,1) +
  xlim(-5,5)

offset_graph <- joined_for_metrics %>%
  group_by(type, offset) %>%
  count() %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(ratio = n / sum(n)) %>%
  ggplot(aes(x = as.integer(offset), y = ratio)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=round(ratio,2)), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x = "Offset From Truth Position (bp)",
       y = "Proportion",
       title = "miRNA Predicted Peaks - sRNAfrag") +
  facet_grid(~factor(type, levels = c("start","end"))) +
  ylim(0,1)

offset_graph_flai <- joined_for_metrics_flai %>%
  group_by(type, offset) %>%
  count() %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(ratio = n / sum(n)) %>%
  ggplot(aes(x = as.integer(offset), y = ratio)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=round(ratio,2)), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x = "Offset From Truth Position (bp)",
       y = "Proportion",
       title = "miRNA Predicted Peaks - FlaiMapper") +
  facet_grid(~factor(type, levels = c("start","end"))) +
  ylim(0,1)

takedownnotes<-joined_for_metrics %>%
  group_by(type, offset) %>%
  count() %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(ratio = n / sum(n))

takedownnotes_flai<-joined_for_metrics_flai %>%
  group_by(type, offset) %>%
  count() %>%
  ungroup() %>%
  group_by(type) %>%
  mutate(ratio = n / sum(n))

# FIND TRUE # WE SHOULD EXPECT

expected_peaks <- sum((even_longer_df %>%
  filter(offset == 0) %>%
  semi_join((joined_for_metrics %>%
               select(source) %>%
               unique()), by = c("derived_id" = "source")) %>%
  group_by(derived_id) %>%
  count())$n)

num_mirnas <- even_longer_df %>%
  filter(offset == 0) %>%
  semi_join((joined_for_metrics %>%
               select(source) %>%
               unique()), by = c("derived_id" = "source")) %>%
  select(derived_id) %>%
  distinct()

plotting_df_no_filter <- labeled_hits %>%
  drop_na()

plotting_df_no_filter$offset <- as.integer(plotting_df_no_filter$offset)

plotting_df_no_filter$dist <- abs(plotting_df_no_filter$offset)

no_filter_dataset_for_offsetgraph <- plotting_df_no_filter %>%
  ggplot(aes(x = as.factor(dist), y = log(counts))) +
  geom_boxplot() +
  labs(x = "Absolute value of offset",
       y = "log2(Counts)",
       title = "Predicted Position Offset vs. Counts")

merged_counts <- read.csv("merged_counts.csv")
ref_table <- read.csv("ref_table.csv")

merged_filtered <- merged_counts %>%
  filter(!str_detect(outside_maps, ",")) %>%
  filter(!str_detect(outside_maps,";")) %>%
  filter(!str_detect(outside_maps, "''")) %>%
  filter(str_detect(outside_maps, "mi"))

merged_filtered$outside_maps <- str_replace_all(merged_filtered$outside_maps, '\\{', "")
merged_filtered$outside_maps <- str_replace_all(merged_filtered$outside_maps, '\\}', "")
merged_filtered$outside_maps <- str_replace_all(merged_filtered$outside_maps, "'", "")

merged_filtered <- merged_filtered %>%
  mutate(total = rowSums(.[2:10]))

merged_combined_filtered <- merged_filtered %>%
    select(outside_maps, total)

merged_combined_filtered <- merged_combined_filtered %>%
  arrange(desc(total))

duplicates <- duplicated(merged_combined_filtered$outside_maps, fromLast = FALSE)

merged_combined_filtered_dedup <- merged_combined_filtered[!duplicates, ]

joined_merged_counts <- merged_combined_filtered_dedup %>%
  inner_join(counts_df, by = c("outside_maps" = "ID"))

counts_metrics_merged <- joined_merged_counts %>%
  ggplot(aes(x = sums, y = total)) +
  geom_point() +stat_poly_line(formula = y~x, method = "lm") +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(use_label("R2"), label.y = 0.8) + 
  labs(x = "True Counts", y = "Pipeline Merged Counts",
       title = "Post-Merge vs. True Counts")

fig3 <- offset_graph / (counts_metrics + counts_metrics_merged) + 
  plot_annotation(tag_levels = 'A', tag_suffix= ")" , title = "Position & Count Calling")

ggsave("fig3.tiff", plot = fig3, height = 9, width = 9, dpi = 600)
ggsave("fig3.jpeg", plot = fig3, height = 9, width = 9, dpi = 600)
