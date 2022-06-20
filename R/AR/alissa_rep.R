# Collection of helper functions to run ExomeDepth in Alissa Reporter
# TODO: use roxygen to document each function
library('dplyr')

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
  coverage_columns <- get_coverage_columns(fixed_columns, my_counts)
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
perform_cnv_calling <- function(my_counts_file, target_sample, ref_samples, file_dir, bias_correction = TRUE, transition.probability = 10^-4){
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
  #Perform gc-bias correction, if needed
  if (bias_correction){
    cat("GC content incorporated into the model")
    data <- data.frame(GC = my_counts$GC)
    model <- 'cbind(test, reference) ~ GC'
  }else{
    model <- 'cbind(test, reference) ~ 1'
    data <- NULL
  }
  all.exons <- new('ExomeDepth', data = data, test = my.test,
                   reference = my.reference.selected, 
                   formula = model)
  all.exons <- CallCNVs(x = all.exons, transition.probability = transition.probability,
                        chromosome = my_counts$chromosome, start = my_counts$start,
                        end = my_counts$end,
                        name = my_counts$names)
  # Write CNV calls to file
  file_name <- paste(strsplit(target_sample, "\\.")[[1]][1], '_cnv.txt', sep = "")
  cnv_calls_file <- file.path(file_dir, file_name)
  write.table(all.exons@CNV.calls, cnv_calls_file, 
              quote=FALSE, sep='\t', row.names = FALSE)
  # Write calculated dq-values to a file
  file_name <- paste(strsplit(target_sample, "\\.")[[1]][1], '_dq.txt', sep = "")
  dq_file <- file.path(file_dir, file_name)
  dq_df <- all.exons@annotations
  dq_df$test <- all.exons@test
  dq_df$reference <- all.exons@reference
  dq_df$ratio_observed <- all.exons@test/ (all.exons@reference + all.exons@test)
  dq_df$ratio_expected <- all.exons@expected
  pc <- 1e-16 #pseudocount to avoid zero division
  dq_df$dq <- (dq_df$ratio_observed + pc) / (all.exons@expected + pc) 
  write.table(dq_df, dq_file, 
              quote=FALSE, sep='\t', row.names = FALSE)
  return(list(cnv_calls_file, dq_file))
}

add_poisson_noise <- function(my_counts_file, lambda_vec){
  # Add poisson noise to the counts matrix
  my.counts <- read.table(my_counts_file, header = TRUE, sep = "\t", quote = "")
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  coverage_columns <- get_coverage_columns(fixed_columns, my.counts)
  my.counts.matrix <- my.counts[,coverage_columns]
  #lambda = c(10, 10, 50, 10, 10) # add different amounts of poisson noise to the samples
  n_els = nrow(my.counts.matrix)
  if (length(lambda_vec) == n_els){
    #Check whether the lambda vector has the same length as the number of coverage columns in the file
    my.counts.noise = matrix(, nrow = n_els, ncol = length(file_names), dimnames = list(c(NULL), as.vector(file_names)))
    for (i in 1:length(file_names)){
      my.counts.noise[,i] = my.counts.matrix[,i] +  rpois(n_els, lambda[i])
    }
    my.counts <- cbind(my.counts[,c("chromosome", "start", "end", "GC")],as.data.frame(my.counts.noise))
    rm(list = c("my.counts.matrix", "my.counts.noise"))
    #TODO: overwrite the existing coverage files or rename them?
    return(my.counts)
  }else{
    stop("The lambda vector should have the same length as the number of samples")
  }
}

add_cnv <- function(my_counts_file_in, my_counts_file_out, cnv_file){
  cnvs <- read.table(cnv_file, header = TRUE, sep = "\t", quote = "")
  my.counts <- read.table(my_counts_file_in, header = TRUE, sep = "\t", quote = "")
  #TODO: also that above! 
  names(my.counts) <- gsub("\\.bam", "", names(my.counts)) # remove the .bam extension in the dataframe 
  for (row in 1:nrow(cnvs)){
    cnv.copy_num  <- cnvs[row, "copy_number"]
    cnv_coordinates <- cnvs[row, "cnv_coordinates"]
    cnv.chrom <- strsplit(cnv_coordindates, ":")[[1]][1]
    cnv.start <- strsplit(strsplit(cnv_coordindates, ":")[[1]][2], "-")[[1]][1]
    cnv.end <- strsplit(strsplit(cnv_coordindates, ":")[[1]][2], "-")[[1]][2]
    cnv.sample <- cnvs[row, "sample"]
    if(!sample %in% names(my.counts)){
      stop("Given sample is not present in the counts file")
    }
    # Find the number of bins that fall within the cnv coordinates
    rowSelect <- my.counts$chromosome == cnv.chrom & my.counts$start >= cnv.start & my.counts$start <= cnv.end
    num_cnvs <- length(which(rowSelect))
    if(num_cnvs < 1){
      stop("Given cnv coordinates are not present in the counts file")
    }
    my.counts[rowSelect, cnv.sample] <- my.counts[rowSelect, cnv.sample]*cnv.copy_num/2
  }
  #Overwrite the existing coverage files or rename them?
  write.table(as.data.frame(my.counts), my_counts_file_out, 
              quote=FALSE, sep='\t', row.names = FALSE)
  return(my_counts_file_out)
}

get_coverage_columns <- function(fixed_columns, counts_df){
  coverage_columns <- names(counts_df)[!names(counts_df) %in% fixed_columns]
  return(coverage_columns)
}

