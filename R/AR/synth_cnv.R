suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ExomeDepth"))

option_list = list(
  make_option(c("-i", "--icounts"), type="character", default=NULL, 
              help="path to inputfile with target counts", metavar="character"),
  make_option(c("-o", "--ocounts"), type="character", default=NULL, 
              help="path to outputfile with target counts", metavar="character"),
  make_option(c("-c", "--cnv"), type="character", default=NULL,
              help="file with cnvs to be added to the counts", metavar="character"),
  make_option(c("-n", "--noise"), type="character", default=NULL,
              help="file with poisson noise to be added to the counts", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$icounts)){
  print_help(opt_parser)
  stop("Argument with input counts file should be provided. \n", call.=FALSE)
}

if (is.null(opt$ocounts)){
  print_help(opt_parser)
  stop("Argument with output counts file should be provided. \n", call.=FALSE)
}

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
    if(!c(cnv.sample) %in% names(my.counts)){
      stop("Given sample is not present in the counts file")
    }
    cat("Add CNV with coordinates ", cnv_coordinates, " to sample ", cnv.sample, "\n", sep = "")
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

add_poisson_noise <- function(my_counts_file_in, my_counts_file_out, noise_file){
  # Add poisson noise to the counts matrix
  # Get the counts
  my.counts <- read.table(my_counts_file_in, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  # Get the noise data
  poisson_noise  <- read.table(noise_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  n_els = nrow(my.counts)
  for (row in 1:nrow(poisson_noise)){
    nsample <- poisson_noise[row, "sample"]
    noise <- poisson_noise[row, "noise"]
    cat("Add Poission noise ", noise, " to sample ", nsample, "\n", sep = "")
    if(!c(nsample) %in% names(my.counts)){
      stop("Given sample is not present in the counts file")
    }
    my.counts[,nsample] <- my.counts[,nsample] + rpois(n_els, noise)
  }
  # write the new counts matrix
  write.table(as.data.frame(my.counts), my_counts_file_out, 
              quote=FALSE, sep='\t', row.names = FALSE)
  return(my_counts_file_out)
}

if (!is.null(opt$noise) & !is.null(opt$cnv)){
  counts_file_out <- add_cnv(my_counts_file_in = opt$icounts, my_counts_file_out = opt$ocounts, cnv_file = opt$cnv)
  counts_file_out <- add_poisson_noise(my_counts_file_in = counts_file_out, my_counts_file_out = opt$ocounts, noise_file = opt$noise)
} else if (!is.null(opt$noise)){
  counts_file_out <- add_poisson_noise(my_counts_file_in = opt$icounts, my_counts_file_out = opt$ocounts, noise_file = opt$noise)
} else if (!is.null(opt$cnv)){
  counts_file_out <- add_cnv(my_counts_file_in = opt$icounts, my_counts_file_out = opt$ocounts, cnv_file = opt$cnv)
}
