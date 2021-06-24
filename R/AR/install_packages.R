# Install BiocManager and required packages for ExomeDepth
install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("Biostrings", "IRanges", "Rsamtools", "GenomicRanges", "GenomicAlignments"))

# Install remaining (regular) packages
install.packages("optparse", repos="https://cloud.r-project.org"i)
install.packages("dplyr", repos="https://cloud.r-project.org")
install.packages("devtools", repos="https://cloud.r-project.org")
suppressPackageStartupMessages(library("devtools"))
devtools::install_github("Multiplicom/ExomeDepth", force = TRUE)
