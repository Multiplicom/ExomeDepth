suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ExomeDepth"))

option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="path to file with target counts", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="name of the target sample name for cnv calling", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=NULL,
              help="comma-separated list of reference samples", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="directory containing the output files", metavar="character"),
  make_option(c("-b", "--bias"), action="store_true", type="logical", default=FALSE,
              help="perform GC-bias correction [default= %default]", metavar="logical"),
  make_option(c("-p", "--prob"), type="double", default=1e-4,
              help="transition probability used for cnv calling [default= %default]", metavar="double")  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Argument with counts file should be provided. \n", call.=FALSE)
}

if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("Argument with target sample name should be provided. \n", call.=FALSE)
}

if (is.null(opt$ref)){
  print_help(opt_parser)
  stop("Argument with list of reference sample names should be provided. \n", call.=FALSE)
}

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Argument with output directory should be provided. \n", call.=FALSE)
}

# Function to select the reference samples for the target sample
perform_cnv_calling <- function(my_counts_file, target_sample, ref_samples, file_dir, bias_correction, transition.probability){
  my_counts <- read.table(my_counts_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  fixed_columns <- c('chromosome', 'start', 'end', 'GC', 'names')
  ref_samples_list <- strsplit(ref_samples, split = ",")[[1]]
  my.ref.samples <- names(my_counts)[names(my_counts) %in% ref_samples_list]
  my.test <- my_counts[,target_sample]
  my.reference.set <- as.matrix(my_counts[, my.ref.samples]) 
  my.choice <- ExomeDepth::select.reference.set(test.counts = my.test,
                                    reference.counts = my.reference.set,
                                    bin.length = (my_counts$end - my_counts$start)/10, 
                                    n.bins.reduced = 10000)
  cat("Selected reference samples: ", paste(my.choice$reference.choice, collapse = ", "), "\n", sep = "")
  # TODO: write choice dataframe to file and store this as fileobj?
  my.matrix <- as.matrix( my_counts[, my.choice$reference.choice, drop = FALSE]) # includes all selected samples
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  #Perform gc-bias correction, if needed
  if (bias_correction){
    cat("GC content incorporated into the model \n")
    data <- data.frame(GC = my_counts$GC)
    model <- 'cbind(test, reference) ~ GC'
  }else{
    model <- 'cbind(test, reference) ~ 1'
    data <- NULL
  }
  all.exons <- new('ExomeDepth', data = data, test = my.test,
                   reference = my.reference.selected, 
                   formula = model)
  all.exons <- ExomeDepth::CallCNVs(x = all.exons, transition.probability = transition.probability,
                        chromosome = my_counts$chromosome, start = my_counts$start,
                        end = my_counts$end,
                        name = my_counts$names)
  # Calculate dq-values
  dq_df <- all.exons@annotations
  dq_df$test <- all.exons@test
  dq_df$reference <- all.exons@reference
  dq_df$ratio_observed <- all.exons@test/ (all.exons@reference + all.exons@test)
  dq_df$ratio_expected <- all.exons@expected
  dq_df$dq <-  dq_df$ratio_observed/ all.exons@expected
  # TODO: create separate file which contains the information (i.e. dq) on the specific bins in the cnv
  #for (cnv in 1:nrow(all.exons@CNV.calls)){
  #  select_bins <- which(dq_df$chromosome == cnv$chromosome & dq_df$end >= cnv$start & dq_df$end <= cnv$end)
  #  dq_bin <- paste(round(dq_df[select_bins, "dq"], digits = 3), collapse = ",")
    #mean_dq <- mean(dq_df[select_bins, "dq"]) # this seems to be the same as the reads.ratio in the file, no need to recalculate this?
  #}
  # Write CNV calls to file
  file_name <- paste(strsplit(target_sample, "\\.")[[1]][1], '_cnv.txt', sep = "")
  cnv_calls_file <- file.path(file_dir, file_name)
  df_CNV_calls <- all.exons@CNV.calls
  if(dim(df_CNV_calls)[1] != 0) { # check whether there are cnv calls
    #reorder columns so that chr, start and end are the first three columns (needed by the Bedtools intersect)
    df_CNV_calls <-  df_CNV_calls[c("chromosome", "start", "end", "start.p", "end.p", "type", "nexons", "id", "BF",
                                    "reads.expected", "reads.observed", "reads.ratio")]
    # add `#` as a first character to the column names - to comment the header (also needed by the Bedtools intersect)
    colnames(df_CNV_calls)[1] <- paste("#", colnames(df_CNV_calls)[1], sep="")
  }
  write.table(df_CNV_calls, cnv_calls_file,
              quote=FALSE, sep='\t', row.names = FALSE)
  # Write calculated dq-values to a file
  file_name <- paste(strsplit(target_sample, "\\.")[[1]][1], '_dq.txt', sep = "")
  dq_file <- file.path(file_dir, file_name)

  # reorder columns so that chr, start and end are the first three columns (needed by the Bedtools intersect)
  dq_df <- dq_df[c("chromosome", "start", "end", "name", "test", "reference", "ratio_expected", "ratio_observed", "dq")]
  # add `#` as a first character to the column names - to comment the header (also needed by the Bedtools intersect)
  colnames(dq_df)[1] <- paste("#", colnames(dq_df)[1], sep="")

  write.table(dq_df, dq_file,
              quote=FALSE, sep='\t', row.names = FALSE)
  return(list(cnv_calls_file, dq_file))
}

perform_cnv_calling(my_counts_file=opt$counts, target_sample=opt$sample, 
                                ref_samples=opt$ref, file_dir=opt$out, 
                                bias_correction=opt$bias, transition.probability=opt$prob)







