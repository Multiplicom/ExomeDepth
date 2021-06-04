# Collection of helper functions to run ExomeDepth in Alissa Reporter

# Create a dataframe from the Agilent Regions or Covered bed-file
get_bed_frame <- function(bed_file, BIN_LENGTH_THRESHOLD = 50) {
  bed.frame <- read.delim(file = bed_file, header = FALSE, 
                          stringsAsFactors = FALSE)
  names(bed.frame)[1] <- "chromosome"
  names(bed.frame)[2] <- "start"
  names(bed.frame)[3] <- "end"
  bin.length <- bed.frame$end - bed.frame$start # assumption that the bed-file is 0 based. Correct?
  bed.frame <- bed.frame[which((bed.frame$end - bed.frame$start)>= BIN_LENGTH_THRESHOLD),] # TODO: keep track of bins that were excluded
  return(bed.frame)
}

# Function which creates the target counts dataframe
get_target_counts <- function(bam_dir, reference_fasta, bed_file){
  bed_frame = get_bed_frame(bed_file)
  bam_files = list.files(path=bam_dir, pattern="*.bam$", full.names=TRUE, recursive=FALSE)
  my_counts = ExomeDepth::getBamCounts(bed.frame=bed_frame, bam.files=bam_files, referenceFasta=reference_fasta)
  my_counts$names = paste(my_counts$chromosome, paste(my_counts$start, my_counts$end, sep = "-"), sep = ":")
  return(my_counts)
}

# Function to convert the target counts into text delimited files
get_coverage_files <- function(my_counts, file_dir){
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  coverage_columns <- names(my_counts)[!names(my_counts) %in% fixed_columns]
  for (cov_col in coverage_columns){
    sample_df <- my.counts[,c(fixed_columns, cov_col)]
    file_name <- paste(strsplit(cov_col, "\\.")[[1]][1], '_cov.txt', sep = "")
    out_file <- file.path(file_dir, file_name)
    write.table(sample_df, out_file, quote=FALSE, sep='\t', row.names = FALSE)
  }
}


