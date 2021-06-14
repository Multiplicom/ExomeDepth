# Collection of helper functions to run ExomeDepth in Alissa Reporter
# TODO: use roxygen to document each function

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

# Function which converts the target counts dataframe into a separate file
#dataframe_to_file <- function(input_df, file_name){
#  write.table(input_df,input_df,sep="\t",row.names=FALSE)
#}

# Function to convert the target counts into text delimited files
get_coverage_files <- function(my_counts, file_dir){
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  coverage_columns <- names(my_counts)[!names(my_counts) %in% fixed_columns]
  for (cov_col in coverage_columns){
    sample_df <- my_counts[,c(fixed_columns, cov_col)]
    file_name <- paste(strsplit(cov_col, "\\.")[[1]][1], '_cov.txt', sep = "")
    out_file <- file.path(file_dir, file_name)
    write.table(sample_df, out_file, quote=FALSE, sep='\t', row.names = FALSE)
  }
  my_counts_file <- file.path(file_dir, 'my_counts.txt')
  write.table(my_counts, my_counts_file, 
                                quote=FALSE, sep='\t', row.names = FALSE)
  return(my_counts_file)
}

# Function to select the reference samples for the target sample
perform_cnv_calling <- function(my_counts_file, target_sample, ref_samples, bias_correction = TRUE){
  my_counts <- read.table(my_counts_file, header = TRUE, sep = "\t", quote = "")
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  my.ref.samples <- names(my_counts)[names(my_counts) %in% ref_samples]
  my.test <- my_counts[,target_sample]
  my.reference.set <- as.matrix(my_counts[, my.ref.samples]) 
  my.choice <- select.reference.set(test.counts = my.test,
                                    reference.counts = my.reference.set,
                                    bin.length = (my_counts$end - my_counts$start)/10, 
                                    n.bins.reduced = 10000)
  cat("Selected reference samples: ", paste(my.choice$reference.choice, collapse = ", "), sep = "")
  # TODO: write choice dataframe to file and store this as fileobj?
  my.matrix <- as.matrix( my_counts[, my.choice$reference.choice, drop = FALSE]) # includes all selected samples
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  #TODO: also include here possibility to perform gc-bias correction
  if (bias_correction){
    model <- 'cbind(test, reference) ~ GC'
  }else{
    model <- 'cbind(test, reference) ~ 1'
  }
  all.exons <- new('ExomeDepth', test = my.test,
                   reference = my.reference.selected, 
                   formula = model)
  all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4,
                        chromosome = my_counts$chromosome, start = my_counts$start,
                        end = my_counts$end,
                        name = my_counts$names)
  return(all.exons)
  
}



