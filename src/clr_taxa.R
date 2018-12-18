#!/usr/bin/Rscript

## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(require(zCompositions)) # 0 imputation
suppressMessages(require(compositions)) # clr transformation
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(require(grid))
suppressMessages(require(pheatmap)) # heatmaps
suppressMessages(require(viridis)) # colors

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--otus', help='OTU table after filtering', required=TRUE)
    parser$add_argument('--trans', help='Transformed OTU data', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)
    parser$add_argument('--sample_ids', help='Sample ID map', required=TRUE)
    parser$add_argument('--taxa', help='Complete, correct and unique taxon names', required=TRUE, nargs='+')

    parser$add_argument('--obname', help='Output basename, i.e. w/o extension', required=TRUE)

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

# Transformed values
trans  <- read.csv(file=args$trans,  sep='\t', header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
flog.info(sprintf('trans: %d x %d', nrow(trans), ncol(trans)), name=logger)

# Sample meta info
meta <- read_meta(args$sample_meta, id_map_file=args$sample_ids)
meta <- meta[intersect(rownames(meta),colnames(otus)),]
rownames(trans) <- meta[rownames(trans),'SampleIDNew']
rownames(meta)  <- meta$SampleIDNew
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)

## TAXA
taxa_otus <- NULL
for(taxon in args$taxa){
    taxon_rank <- NULL
    taxon_otus <- NULL
    for(rank in ranks){
        if(any(otus[,rank] == taxon)){
            testit::assert(is.null(taxon_rank))
            taxon_rank <- rank
            taxon_otus <- rownames(otus)[otus[,rank] == taxon]
        }
    }
    flog.info(sprintf('Taxon \"%s\": rank \"%s\", %d OTUs: %s', taxon, taxon_rank, length(taxon_otus), paste(taxon_otus, collapse=', ')), name=logger)
    taxa_otus <- rbind(taxa_otus, data.frame(Taxon=taxon, OTU=taxon_otus, check.names=FALSE, stringsAsFactors=FALSE))
}
rownames(taxa_otus) <- taxa_otus$OTU

## PLOTS
plots <- list()

# Heatmap
plots[['taxa_clr_pheatmap']] <- pheatmap(
    mat=t(trans[,taxa_otus$OTU]),
    color=viridis(30),
    border_color=NA,
    scale="none",
    show_rownames=TRUE,
    show_colnames=TRUE,
    fontsize_col=7,
    clustering_distance_rows='euclidean',
    clustering_distance_cols='euclidean',
    clustering_method='complete',
    cluster_rows=FALSE,
    cluster_cols=TRUE,
    # treeheight_row=0, # too many OTUs to see tree structure
    legend=TRUE,
    annotation_row=taxa_otus[, 'Taxon', drop=FALSE],
    annotation_col=meta[rownames(trans), 'Indication', drop=FALSE],
    annotation_colors=list(
        Indication=meta_colors$Indication
    ),
    main='CLR transformed values',
    silent=TRUE
)
for(taxon in args$taxa){
    plots[[sprintf('taxa_clr_pheatmap_%s', taxon)]] <- pheatmap(
        mat=t(trans[,taxa_otus$OTU[taxa_otus$Taxon == taxon], drop=FALSE]),
        color=viridis(30),
        border_color=NA,
        scale="none",
        show_rownames=TRUE,
        show_colnames=TRUE,
        fontsize_col=7,
        clustering_distance_rows='euclidean',
        clustering_distance_cols='euclidean',
        clustering_method='complete',
        cluster_rows=FALSE,
        cluster_cols=TRUE,
        # treeheight_row=0, # too many OTUs to see tree structure
        legend=TRUE,
        annotation_row=taxa_otus[, 'Taxon', drop=FALSE],
        annotation_col=meta[rownames(trans), 'Indication', drop=FALSE],
        annotation_colors=list(
            Indication=meta_colors$Indication
        ),
        main='CLR transformed values',
        silent=TRUE
    )
}

## PDF
pdf(sprintf('%s.pdf', args$obname), width=9, height=6)
for(pname in names(plots)){
    if(grepl('pheatmap', pname)){
        grid.newpage()
        grid.draw(plots[[pname]]$gtable)
    } else {
        print(plots[[pname]])
    }
}
dev.off()

## RDS
saveRDS(plots, file=sprintf('%s.rds', args$obname))
