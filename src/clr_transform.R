#!/usr/bin/Rscript

## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(tools))
suppressMessages(library(argparse))
suppressMessages(require(zCompositions)) # 0 imputation
suppressMessages(require(compositions)) # clr transformation

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--otus', help='OTU table after filtering', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)

    parser$add_argument('--ofile', help='Outputfile after normalization/scaling', required=TRUE)

    parser$add_argument('--src', help='Source directory', default="src")
    return(parser)
}

## ARGS
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## UTILS
source(sprintf('%s/utils.R', args$src))

## LOG
logfile <- sprintf('%s.log', tools::file_path_sans_ext(args$ofile))
logger <- add_logger(
    logger_file=logfile,
    logger_name='mylogger',
    logger_level='INFO'
)

# log session
log_session_info(logger)

# log args
log_args(args, logger)

# log file time stamps
ifiles <- c(args$otus, args$sample_meta)
flog.info(sprintf(
    "STATUS: Input files: variable name, path, timestamp\n%s",
    paste(ifiles, sapply(ifiles, fmtime_str), sep=': ', collapse='\n')
), name=logger)

## DATA
# OTUs
otus <- read_otus(otus_file=args$otus, flags=c('kept'))
rownames(otus) <- otus$OTU
flog.info(sprintf('OTUs: %d x %d', nrow(otus), ncol(otus)), name=logger)
# Sample meta info
meta <- read_meta(args$sample_meta, c('SampleID', 'Indication'))
meta <- meta[intersect(rownames(meta),colnames(otus)),]
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)

## COUNTS
counts <- t(otus[,rownames(meta)]) # samples x features
flog.info(sprintf('OTU counts: only samples from meta table and non-fitered OTUs: %d x %d', nrow(counts), ncol(counts)), name=logger)

## CLR transformation
counts_clr <- counts2clr(counts)
write.table(x=counts_clr, file=args$ofile, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
