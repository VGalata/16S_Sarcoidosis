#!/usr/bin/Rscript

## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(library(ggplot2))
suppressMessages(require(ggtree))
suppressMessages(library(tools))
suppressMessages(library(gtools))

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--otu_counts', '-c', help='OTU.txt', required=TRUE)
    parser$add_argument('--otu_tax', '-t', help='hiera_RDP.txt', required=TRUE)
    parser$add_argument('--otu_blastn', '-b', help='Processed BLASTn hits', required=TRUE)
    parser$add_argument('--otu_tree', '-p', help='Tree.tre', required=TRUE)
    parser$add_argument('--sample_meta', '-m', help='Sample meta information', required=TRUE)
    parser$add_argument('--sample_ids', '-i', help='Sample ID map', required=TRUE)

    parser$add_argument('--otab', help='Output table file for OTUs', required=TRUE)
    parser$add_argument('--opdf', help='Output PDF file for plots', required=TRUE)
    parser$add_argument('--ords', help='Output RDS file for plots', required=TRUE)

    parser$add_argument('--contam_expr', help='Regular expression for contamination taxon', default='Eukaryota')
    parser$add_argument('--contam_rank', help='Taxonomy column to search for contamination taxa', default='taxon_superkingdom_name')
    parser$add_argument('--contam_case', help='Whether contamination expression is case INsensitive', action="store_true")

    parser$add_argument('--min_count', help='Minimal OTU count', default=1, type='integer')

    parser$add_argument('--ignor', help='Sample IDs to ignore', default=c(), nargs='*')
    parser$add_argument('--ctrls', help='Sample IDs of controls', default=c(), nargs='*')

    parser$add_argument('--width', help='PDF width', default=12, type='integer')
    parser$add_argument('--height', help='PDF height', default=6, type='integer')
    parser$add_argument('--src', help='Source directory', default="src")
    return(parser)
}

# Filter contaminated OTUs by taxonomy
get_otus_contam <- function(tax, contam_expr, contam_rank, contam_case){
    testit::assert(contam_rank %in% colnames(tax))
    ids <- as.character(
        rownames(tax)[  grepl(contam_expr, tax[,contam_rank], ignore.case=contam_case) ]
    )
    return(ids)
}

## ARGS
args <- get_argparser()$parse_args(commandArgs(trailingOnly=TRUE))

## UTILS
source(sprintf('%s/utils.R', args$src))

## LOG
logfile <- sprintf('%s.log', tools::file_path_sans_ext(args$otab))
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
ifiles <- c(args$otu_counts, args$otu_tax, args$otu_blastn, args$otu_tree, args$sample_meta, args$sample_ids)
flog.info(sprintf(
    "STATUS: Input files: variable name, path, timestamp\n%s",
    paste(ifiles, sapply(ifiles, fmtime_str), sep=': ', collapse='\n')
), name=logger)

## DATA
# OTUs
counts <- read.csv(file=args$otu_counts, sep='\t', header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
flog.info(sprintf('Counts: %d x %d', nrow(counts), ncol(counts)), name=logger)

ltaxa  <- read.csv(file=args$otu_tax, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
rownames(ltaxa) <- ltaxa$OTU
ltaxa  <- ltaxa[,setdiff(colnames(ltaxa), 'OTU')]
flog.info(sprintf('LotuS taxa: %d x %d', nrow(ltaxa), ncol(ltaxa)), name=logger)

btaxa  <- read.csv(file=args$otu_blastn, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE); rownames(btaxa) <- btaxa$qseqid
flog.info(sprintf('BLAST taxa: %d x %d', nrow(btaxa), ncol(btaxa)), name=logger)

otu_ids <- rownames(counts)

otree  <- ape::read.tree(file=args$otu_tree)

# Sample meta info
meta <- read_meta(
    meta_file=args$sample_meta,
    cols=c('SampleID', 'Indication'),
    id_map_file=args$sample_ids
)
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)
samples <- setdiff(colnames(counts),c(args$ctrls, args$ignor))

## FILTERING
# By taxonomy
otus_contam <- get_otus_contam(tax=btaxa, contam_expr=args$contam_expr, contam_rank=args$contam_rank, contam_case=args$contam_case)
flog.info(sprintf('Filtering OTUs: by BLAST taxonomy: %d', length(otus_contam)), name=logger)

# By controls
otus_inctrl <- c()
if(length(args$ctrls) > 0){
    otus_inctrl <- rownames(counts)[ apply(counts[,args$ctrls,drop=FALSE], 1, function(x){ any(x > 0) }) ]
}
flog.info(sprintf('Filtering OTUs: by presence in controls: %d', length(otus_inctrl)), name=logger)

# By min. count
otus_mincount <- rownames(counts)[ apply(counts[,samples,drop=FALSE], 1, function(x){ all(x < args$min_count) }) ]
flog.info(sprintf('Filtering OTUs: by minimal count: %d', length(otus_mincount)), name=logger)

# Filtering flag
otus_flag <- sapply(otu_ids, function(x){
    flags <- c()
    if(x %in% otus_contam){
        flags <- c(flags, 'contamination')
    }
    if(x %in% otus_inctrl){
        flags <- c(flags, 'controls')
    }
    if(x %in% otus_mincount){
        flags <- c(flags, 'count')
    }
    if(length(flags) == 0){
        return('kept')
    } else {
        return(paste(sort(flags), collapse=';'))
    }
}, USE.NAMES=TRUE)

## RESULT
otus <- data.frame(
    OTU=otu_ids,
    ltaxa[otu_ids,],
    BLASTN_NT=btaxa[otu_ids, 'lineage'],
    FLAG=otus_flag[otu_ids],
    counts[otu_ids,samples,drop=FALSE],
    row.names=otu_ids, check.names=FALSE, stringsAsFactors=FALSE
)

## TAB
write.table(x=otus, file=args$otab, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

## PLOTS
plots <- list()

# Perc. of OTU counts per sample affected by filtering
df <- NULL
for(x in samples){
    tmp <- aggregate(counts[,x], by=list(FLAG=otus[rownames(counts),'FLAG']), FUN=sum)
    tmp$x <- 100 * tmp$x / sum(counts[,x])
    tmp$s <- x
    colnames(tmp) <- c('Flag', 'Pct', 'Sample')
    if(is.null(df)){
        df <- tmp
    } else {
        df <- rbind(df, tmp)
    }
}
df$SampleIDNew <- meta[df$Sample,'SampleIDNew'] # add new sample IDs
df$SampleIDNew <- factor(x=df$SampleIDNew, levels=gtools::mixedsort(unique(df$SampleIDNew)), ordered=TRUE)
plots[['otu_filtering_pct_count_per_sample']] <-
    ggplot(data=df, aes(x=SampleIDNew, y=Pct, fill=Flag)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=flag_colors, name='Filtering flag') +
    xlab("") + ylab("Percentage of OTU counts [%]") + ggtitle("") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, size=10),
        axis.ticks.x=element_blank()
    )

# OTU tree with flag highlighting
df <- data.frame(
    OTU=rownames(otus),
    Flag=otus$FLAG,
    row.names=rownames(otus),
    stringsAsFactors=FALSE
)
plots[['otu_filtering_otu_tree']] <-
    ggtree(otree, layout="rectangular") %<+% df +
    theme_tree2() +
    geom_tippoint(aes(colour=Flag)) +
    scale_colour_manual(values=flag_colors, name='Filtering flag') +
    ggtitle("LotuS tree") +
    theme(legend.position="right")

## PDF
pdf(args$opdf, width=args$width, height=args$height)
for(pname in names(plots)){
    print(plots[[pname]])
}
dev.off()

## RDS
saveRDS(plots, file=args$ords)
