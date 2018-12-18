#!/usr/bin/Rscript

## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
print(warnings())

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--multiqc_stats', help='MultiQC stat.s file (multiqc_general_stats.txt)', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)
    parser$add_argument('--sample_ids', help='Sample ID map', required=TRUE)

    parser$add_argument('--obname', help='Output basename, i.e. w/o extension', required=TRUE)

    parser$add_argument('--width', help='PDF width', default=6, type='integer')
    parser$add_argument('--height', help='PDF height', default=6, type='integer')
    parser$add_argument('--src', help='Source directory', default="src")
    return(parser)
}

## ARGS
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## UTILS
source(sprintf('%s/utils.R', args$src))

## LOG
logfile <- sprintf('%s.log', args$obname)
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
ifiles <- c(args$multiqc_stats, args$sample_meta, args$sample_ids)
flog.info(sprintf(
    "STATUS: Input files: variable name, path, timestamp\n%s",
    paste(ifiles, sapply(ifiles, fmtime_str), sep=': ', collapse='\n')
), name=logger)

## DATA
# MultiQC stat.s
multiqc <- read.csv(file=args$multiqc_stats, sep='\t', stringsAsFactors=FALSE, header=TRUE)
multiqc$SampleID <- sapply(multiqc$Sample, function(x){ unlist(strsplit(x,'_'))[1] })
multiqc$Amplicon <- sapply(multiqc$Sample, function(x){ unlist(strsplit(x,'_'))[2] })
multiqc$ReadID   <- sapply(multiqc$Sample, function(x){ unlist(strsplit(x,'_'))[3] })
# Sample meta info
meta <- read_meta(
    meta_file=args$sample_meta,
    cols=c('SampleID', 'Indication'),
    id_map_file=args$sample_ids
)
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)

## PLOTS
plots <- list()

## Numberof reads per sample/group
# data: w/o NTC and only R1; add Indication
readcount <- multiqc[multiqc$ReadID == 'R1' & multiqc$SampleID != "NTC", c('SampleID', 'Amplicon', 'FastQC_mqc.generalstats.fastqc.total_sequences')]
readcount$Indication <- meta[readcount$SampleID, 'Indication']
for(amplicon in sort(unique(readcount$Amplicon))){
    flog.info(sprintf('Median read count for %s: %.2f', amplicon, median(readcount[readcount$Amplicon==amplicon, 'FastQC_mqc.generalstats.fastqc.total_sequences'])), name=logger)
    for(indi in sort(unique(readcount$Indication))){
        flog.info(sprintf('Median read count for %s, %s: %.2f', amplicon, indi, median(readcount[readcount$Amplicon==amplicon & readcount$Indication==indi, 'FastQC_mqc.generalstats.fastqc.total_sequences'])), name=logger)
    }
}
# plot
set.seed(42)
plots[['readcount']] <-
    ggplot(readcount, aes(x=Indication, y=FastQC_mqc.generalstats.fastqc.total_sequences)) +
    geom_boxplot(fill=NA, outlier.shape = NA) + # no outliers
    geom_jitter(shape=21, size=2, fill=NA, width=0.25, colour="#666666") +
    facet_wrap(~Amplicon, ncol=2) + #, scales="free_y") +
    theme_classic() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background =element_rect(fill="#ffffff")
    ) +
    xlab('Indication') + ylab('Read count') + ggtitle('')

## PDF
pdf(sprintf('%s.pdf', args$obname), width=args$width, height=args$height)
for(pname in names(plots)){
    print(plots[[pname]])
}
dev.off()

## RDS
saveRDS(plots, file=sprintf('%s.rds', args$obname))
