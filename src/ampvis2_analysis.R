#!/usr/bin/Rscript

## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(require(ampvis2))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(require(ggsignif)) # to plot h. test p-value
suppressMessages(require(coin)) # for wilcox test with ties
print(warnings())

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--otus', help='OTU table after filtering', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)

    parser$add_argument('--obname', help='Output basename, i.e. w/o extension', required=TRUE)

    parser$add_argument('--width', help='PDF width', default=6, type='integer')
    parser$add_argument('--height', help='PDF height', default=6, type='integer')
    parser$add_argument('--src', help='Source directory', default="src")
    return(parser)
}

wilcox_wrapper <- function(x, y){
    test_result <- coin::wilcox_test(
        value ~ label,
        data.frame(
            value=c(x, y),
            label=c(rep('x', length(x)), rep('y', length(y)))
        ),
        alternative="two.sided",
        conf.level=0.95
    )
    return(list(
        p.value=pvalue(test_result)
    ))
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
ifiles <- c(args$otus,args$sample_meta)
flog.info(sprintf(
    "STATUS: Input files: variable name, path, timestamp\n%s",
    paste(ifiles, sapply(ifiles, fmtime_str), sep=': ', collapse='\n')
), name=logger)

## DATA
# OTUs
otus <- read_otus(otus_file=args$otus, flags=c('kept'))
flog.info(sprintf('OTUs: %d x %d', nrow(otus), ncol(otus)), name=logger)
# Sample meta info
meta <- read_meta(
    meta_file=args$sample_meta,
    cols=c('SampleID', 'Indication')
)
meta <- meta[intersect(rownames(meta),colnames(otus)),]
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)

## AMPVIS2
# OTUs in expected format
otus_           <- otus[,c('OTU', rownames(meta), ranks)]
colnames(otus_) <- c('OTU', rownames(meta), ranks_ampvis)
# create ampvis2 object
amp <- ampvis2::amp_load(otus_,meta)

## PLOTS
plots <- list()

## Rarecurves
plots[['amp_rarecurve']] <-
    amp_rarecurve(
        data=amp,
        color_by='Indication',
        stepsize=500
    ) +
    facet_wrap(~Indication) +
    guides(colour=FALSE) +
    ggtitle('Rarefaction curves')

## Alpha-div.
adiv <- amp_alphadiv(amp, measure = c("observed", "shannon", "simpson", "invsimpson"))
adiv_ <- reshape2::melt(data=adiv, id.vars=c('SampleID', 'Indication', 'Reads')) # for the plot
# save table
write.table(x=adiv, file=sprintf('%s_alphadiv.tsv', args$obname), row.names=FALSE)
# plot
set.seed(42)
plots[['amp_alphadiv']] <-
    ggplot(adiv_, aes(x=Indication, y=value)) +
    geom_boxplot(fill=NA, outlier.shape = NA) + # no outliers
    geom_jitter(shape=21, size=2, fill=NA, width=0.25, colour="#666666") +
    facet_wrap(~variable, ncol=2, scales="free_y") +
    theme_bw() +
    theme(
        axis.text.x=element_text(angle=90, hjust=1)
    ) +
    ggtitle('Alpha-diversity') +
    geom_signif(
        comparisons=get_comparisons_for_ggsignif(adiv_$Indication),
        test=wilcox_wrapper,
        map_signif_level=TRUE,
        step_increase=0.25,
        vjust=1.5
    )

## PDF
pdf(sprintf('%s.pdf', args$obname), width=args$width, height=args$height)
for(pname in names(plots)){
    print(plots[[pname]])
}
dev.off()

## RDS
saveRDS(plots, file=sprintf('%s.rds', args$obname))
