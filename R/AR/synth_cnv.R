suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ExomeDepth"))

option_list = list(
  make_option(c("-i", "--icounts"), type="character", default=NULL, 
              help="path to inputfile with target counts", metavar="character"),
  make_option(c("-o", "--ocounts"), type="character", default=NULL, 
              help="path to outputfile with target counts", metavar="character"),
  make_option(c("-c", "--cnv"), type="character", default=NULL,
              help="file with cnvs to be added to the counts", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


add_cnv <- function(my_counts_file_in, my_counts_file_out, cnv_file){
  cnvs <- read.table(cnv_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors=FALSE)
  my.counts <- read.table(my_counts_file_in, header = TRUE, sep = "\t", quote = "", stringsAsFactors=FALSE)
  #names(my.counts) <- gsub("\\.bam", "", names(my.counts)) # remove the .bam extension in the dataframe 
  for (row in 1:nrow(cnvs)){
    cnv.copy_num  <- cnvs[row, "copy_number"]
    cnv_coordinates <- cnvs[row, "cnv_coordinates"]
    cnv.chrom <- strsplit(cnv_coordinates, ":")[[1]][1]
    cnv.start <- strsplit(strsplit(cnv_coordinates, ":")[[1]][2], "-")[[1]][1]
    cnv.end <- strsplit(strsplit(cnv_coordinates, ":")[[1]][2], "-")[[1]][2]
    cnv.sample <- cnvs[row, "sample"]
    cat("Add CNV with coordinates ", cnv_coordinates, " to sample ", cnv.sample, "\n", sep = "")
    if(!c(cnv.sample) %in% names(my.counts)){
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

add_cnv(my_counts_file_in = opt$icounts, my_counts_file_out = opt$ocounts, cnv_file = opt$cnv)