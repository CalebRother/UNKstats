# raw-data/prepare_datasets.R
# ---------------------------
# This script converts raw CSV data into preloaded datasets
library(usethis)

# Read in your dataset from CSV (adjust filenames if needed for other files)
afromontane <- read.csv("raw-data/TableS1_Species-PAM.csv")

# Save it as a compressed .rda object inside the package's /data folder
# so users can call data(afromontane)
use_data(afromontane, overwrite = TRUE)
