#!/usr/bin/env Rscript

library("optparse")

# load the functions from the script
source('tools/DeepRibo/src/s_curve_cutoff_estimation.R')

# commandline parser
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="dataset file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# list the dataset and the path to which the png figure is stored
get_cutoff_values(path=opt$file, dest="figure")
#$min_RPKM
#  ....
#$min_coverage
#   	  ....

