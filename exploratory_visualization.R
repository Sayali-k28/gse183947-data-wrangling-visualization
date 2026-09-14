# script to visualize gene expression data (GSE183947)
# each section below makes ONE plot and saves it as a PNG file

library(tidyverse)

# load the tidy table made by data_wrangling.R
dat.long <- read.csv("data/dat_long.csv")


### 1. Density plot: is BRCA1 different between tumor and normal? ----------

brca1_data <- dat.long %>%
  filter(gene == "BRCA1")

density_plot <- ggplot(brca1_data, aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3) +
  labs(title = "BRCA1 expression: tumor vs normal")

ggsave("density_BRCA1.png", density_plot, width = 6.5, height = 4.5)


### 2. Violin plot: does BRCA1 differ by metastasis status? ----------------
# note: metastasis status is mostly only recorded for tumor samples, so
# this plot is closely related to the density plot above, not a fully
# separate finding.

violin_plot <- ggplot(brca1_data, aes(x = metastasis, y = FPKM, fill = metastasis)) +
  geom_violin() +
  labs(title = "BRCA1 expression by metastasis status")

ggsave("violin_BRCA1_metastasis.png", violin_plot, width = 5.5, height = 4.5)


### 3. Scatter plot: are BRCA1 and BRCA2 related? --------------------------

brca_wide <- dat.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM)

# correlation = a single number from -1 to 1 that says how closely two
# things move together. Near 0 means "not much relationship."
correlation <- cor(brca_wide$BRCA1, brca_wide$BRCA2)
correlation <- round(correlation, 2)
cat("BRCA1 vs BRCA2 correlation:", correlation, "\n")

scatter_plot <- ggplot(brca_wide, aes(x = BRCA1, y = BRCA2, colour = tissue)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "BRCA1 vs BRCA2 expression",
       subtitle = paste("correlation =", correlation))

ggsave("scatter_BRCA1_BRCA2.png", scatter_plot, width = 7, height = 5)


### 4. Heatmap: compare 5 genes at once ------------------------------------
# problem: some genes (like TP53) are naturally expressed much higher than
# others. If we colour by raw FPKM, TP53 drowns out everything else.
# fix: for each gene, convert its values to a "z-score" -- how far each
# sample is from THAT GENE's own average. This makes every gene fairly
# comparable on the same colour scale.

genes_of_interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")

heatmap_data <- dat.long %>%
  filter(gene %in% genes_of_interest) %>%
  group_by(gene) %>%
  mutate(z_score = (FPKM - mean(FPKM)) / sd(FPKM)) %>%
  ungroup() %>%
  arrange(tissue)

# keep samples grouped by tissue instead of a random order
heatmap_data$samples <- factor(heatmap_data$samples, levels = unique(heatmap_data$samples))

heatmap_plot <- ggplot(heatmap_data, aes(x = samples, y = gene, fill = z_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
  labs(title = "Gene expression heatmap (scaled per gene)") +
  theme(axis.text.x = element_blank())  # 60 sample names would be unreadable

ggsave("heatmap_scaled.png", heatmap_plot, width = 10, height = 4)

cat("Done! 4 plots saved.\n")
