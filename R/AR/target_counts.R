suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ExomeDepth"))


option_list = list(
  make_option(c("-b", "--bam"), type="character", default=NULL, 
              help="directory containing the bam files", metavar="character"),
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="path to reference fasta file", metavar="character"),
  make_option(c("-d", "--bed"), type="character", default=NULL,
              help="path to bed file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="directory containing the output files", metavar="character"),
  make_option(c("-m", "--mapq"), type="double", default=20,
              help="minimum mapping quality for reads to be included [default= %default]", metavar="double")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("Argument with bam-file directory should be provided. \n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Argument with output directory should be provided. \n", call.=FALSE)
}

if (is.null(opt$fasta)){
  print_help(opt_parser)
  stop("Argument with reference fasta file should be provided. \n", call.=FALSE)
}

if (is.null(opt$bed)){
  print_help(opt_parser)
  stop("Argument with bed file should be provided. \n", call.=FALSE)
}

# Create a dataframe from the Agilent Regions or Covered bed-file
get_bed_frame <- function(bed_file) {
  bed.frame <- read.delim(file = bed_file, header = FALSE, 
                          stringsAsFactors = FALSE)
  names(bed.frame)[1] <- "chromosome"
  names(bed.frame)[2] <- "start"
  names(bed.frame)[3] <- "end"
  return(bed.frame)
}

# Function which creates the target counts dataframe
get_target_counts <- function(bam_dir, reference_fasta, bed_file, min_mapq=20){
  bed_frame = get_bed_frame(bed_file)
  bam_files = list.files(path=bam_dir, pattern="*.bam$", full.names=TRUE, recursive=FALSE)
  my_counts = ExomeDepth::getBamCounts(bed.frame=bed_frame, bam.files=bam_files, referenceFasta=reference_fasta, min.mapq = min_mapq)
  my_counts$names = paste(my_counts$chromosome, paste(my_counts$start, my_counts$end, sep = "-"), sep = ":")
  names(my_counts) <- gsub("\\.bam", "", names(my_counts)) # remove the .bam extension in the dataframe 
  return(my_counts)
}

# Function to convert the target counts into text delimited files
get_coverage_files <- function(my_counts, file_dir){
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  coverage_columns <- get_coverage_columns(fixed_columns, my_counts)
  # add `#` as a first character to the column names - to comment the header (needed by the Bedtools intersect)
  coverage_columns[1] <- paste("#", coverage_columns[1], sep="")
  for (cov_col in coverage_columns){
    sample_df <- my_counts[,c(fixed_columns, cov_col)]
    file_name <- paste(strsplit(cov_col, "\\.")[[1]][1], '_cov.txt', sep = "") #removes .bam extension, but probably no longer needed as is been taken care of above
    out_file <- file.path(file_dir, file_name)
    cat("Write counts of ", cov_col, " to output file \n", sep = "")
    cat("Output file: ", out_file, "\n", sep = "")
    cat("Number of lines in sample_df: ", nrow(sample_df), "\n", sep = "")
    write.table(sample_df, out_file, quote=FALSE, sep='\t', row.names = FALSE)
  }
  my_counts_file <- file.path(file_dir, 'my_counts.txt')
  write.table(my_counts, my_counts_file, 
              quote=FALSE, sep='\t', row.names = FALSE)
  return(my_counts_file)
}

get_coverage_columns <- function(fixed_columns, counts_df){
  coverage_columns <- names(counts_df)[!names(counts_df) %in% fixed_columns]
  return(coverage_columns)
}

my_counts <- get_target_counts(bam_dir=opt$bam, reference_fasta=opt$fasta, bed_file=opt$bed, min_mapq=opt$mapq)
my_counts_file <- get_coverage_files(my_counts, file_dir=opt$out)
