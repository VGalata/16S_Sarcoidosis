## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(require(grid))
suppressMessages(require(ggplot2))
suppressMessages(require(viridis)) # colors
suppressMessages(require(ggtree)) # plot trees
suppressMessages(require(pheatmap)) # heatmaps
suppressMessages(require(ggfortify)) # PCA biplot
print(warnings())

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--trans', help='Transformed OTU data', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)
    parser$add_argument('--sample_annot', help='Additional sample annotation information', required=TRUE)
    parser$add_argument('--sample_ids', help='Sample ID map', required=TRUE)

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
# ifiles <- c(args$otus, args$otu_tree, args$trans, args$sample_meta, args$sample_annot, args$sample_ids)
ifiles <- c(args$trans, args$sample_meta, args$sample_annot, args$sample_ids)
flog.info(sprintf(
    "STATUS: Input files: variable name, path, timestamp\n%s",
    paste(ifiles, sapply(ifiles, fmtime_str), sep=': ', collapse='\n')
), name=logger)

# Transformed values
trans  <- read.csv(file=args$trans,  sep='\t', header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
flog.info(sprintf('trans: %d x %d', nrow(trans), ncol(trans)), name=logger)

# Sample meta info
meta <- read_meta(
    meta_file=args$sample_meta,
    # cols=c('F-Index', 'R-Index', 'SampleID', 'Indication', 'FID'),
    cols=c('F-Index', 'R-Index', 'SampleID', 'Indication'),
    id_map_file=args$sample_ids
)
# meta <- meta[intersect(rownames(meta),colnames(otus)),]
meta <- meta[intersect(rownames(meta),rownames(trans)),]
rownames(trans) <- meta[rownames(trans),'SampleIDNew']
rownames(meta)  <- meta$SampleIDNew
flog.info(sprintf('Meta: %d x %d', nrow(meta), ncol(meta)), name=logger)
# Add. sample annotation info
annot <- read_annot(annot_file=args$sample_annot)
# Merge both
meta2 <- merge_meta_annot(meta, annot)
testit::assert(all(meta2$Indication == meta2$Group))

## Plots
plots <- list()
# plot_cols <- c('Indication', 'Zugangsweg', 'Smoker', 'CD4_CD8', 'Sex', 'Comorbidity', 'Comedication')
plot_cols <- c('Indication', 'Smoker')

## PCA
# Reference: https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
proj_pca <- prcomp(x=trans, center=TRUE, scale.=FALSE)
pdf(sprintf('%s_trans_pca.pdf', args$obname), width=7, height=5)
for(col in plot_cols){
    plots[[sprintf('trans_pca_%s', col)]] <-
        autoplot(proj_pca, data=meta2, colour=col) +
        ggtitle('PCA') +
        theme_bw()
    plots[[sprintf('trans_pca_labels_%s', col)]] <-
        autoplot(proj_pca, data=meta2, colour=col, shape=FALSE, label.size=2) +
        ggtitle('PCA') +
        theme_bw()
    plots[[sprintf('trans_pca_biplot_%s', col)]] <-
        autoplot(proj_pca, data=meta2, colour=col, loadings=TRUE) +
        ggtitle('PCA biplot') +
        theme_bw()
    print(plots[[sprintf('trans_pca_%s', col)]])
    print(plots[[sprintf('trans_pca_labels_%s', col)]])
    print(plots[[sprintf('trans_pca_biplot_%s', col)]])
}
dev.off()

## Heatmap + hclust
# Reference: https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf
# How to plot: https://support.rstudio.com/hc/en-us/community/posts/239529128-Notebooks-grid-newpage
pdf(sprintf('%s_trans_pheatmap.pdf', args$obname), width=7, height=10, onefile=FALSE)
plots[['trans_pheatmap']] <- pheatmap(
    mat=t(trans),
    color=viridis(30),
    scale="none",
    show_rownames=FALSE,
    show_colnames=TRUE,
    fontsize_col=7,
    clustering_distance_rows='euclidean',
    clustering_distance_cols='euclidean',
    clustering_method='complete',
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    treeheight_row=0, # too many OTUs to see tree structure
    legend=TRUE,
    annotation_col=meta2[rownames(trans),plot_cols,drop=FALSE],
    annotation_colors=list(
        Indication=meta_colors$Indication,
        Smoker=meta_colors$Smoker
    ),
    main='Transformed values',
    silent=TRUE
)
grid.newpage()
grid.draw(plots[['trans_pheatmap']]$gtable)
dev.off()

## RDS
saveRDS(plots, file=sprintf('%s.rds', args$obname))
