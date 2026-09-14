# Fetch the GSE183947 FPKM matrix directly from GEO's supplementary files.
# Replaces a hardcoded Desktop path with something that works for anyone
# who clones this repo.

library(GEOquery)

dir.create("data", showWarnings = FALSE)

files <- getGEOSuppFiles("GSE183947", baseDir = "data", filter_regex = "fpkm")
gz_path <- rownames(files)[1]

# unzip into data/GSE183947_fpkm.csv
R.utils::gunzip(gz_path, destname = "data/GSE183947_fpkm.csv", remove = FALSE, overwrite = TRUE)

# tidy up GEOquery's nested folder (data/GSE183947/GSE183947_fpkm.csv.gz)
unlink("data/GSE183947", recursive = TRUE)

cat("Wrote data/GSE183947_fpkm.csv\n")
