#script to clean and reshape gene expression data (GSE183947) into a tidy table
#load libraries
library(dplyr)
library(tidyverse)
library (GEOquery)

#read in the data
# run fetch_data.R first -- it downloads this file straight from GEO's
# supplementary files, so this path works on any machine, not just the
# one that originally ran the analysis.
dat <- read.csv(file = "data/GSE183947_fpkm.csv")
dim(dat)

# get metadata

gse <- getGEO(GEO ='GSE183947' , GSEMatrix = TRUE)

gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

metadata.modified <- metadata %>%
  select(1,10,11,17)%>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue=trimws(gsub("tissue:","", tissue))) %>%
  mutate(metastasis=trimws(gsub("metastasis:","", metastasis)))

head(dat)

#reshaping data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = "samples", value = "FPKM", -gene)

#join dataframes = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# explore data

dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM),
           median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)

#save the tidy long-format table so data_visualization.R (or anything else)
#can pick it up without needing this script's R session in memory
write.csv(dat.long, "data/dat_long.csv", row.names = FALSE)
cat("Wrote data/dat_long.csv (", nrow(dat.long), "rows )\n")

  
