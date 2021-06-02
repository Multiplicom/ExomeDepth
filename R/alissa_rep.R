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
  bam_files = list.files(path=bam_dir, pattern="*.bam$", full.names=TRUE, recursive=False)
  my_counts = exomedepth.getBamCounts(bed_frame=bed_frame, bam_files=bam_files, referenceFasta=reference_fasta)
  return(my_counts)
}


