## IMPORTS
suppressMessages(library(futile.logger))
suppressMessages(library(argparse))
suppressMessages(require(ALDEx2)) # diff. abundance analysis

suppressMessages(require(ggplot2))
suppressMessages(require(patchwork)) # combine plots
suppressMessages(require(scales))
print(warnings())

## FUNCTIONS
# Argument parser
get_argparser <- function(){
    parser <- ArgumentParser()
    parser$add_argument('--otus', help='OTU table after filtering', required=TRUE)
    parser$add_argument('--sample_meta', help='Sample meta information', required=TRUE)

    parser$add_argument('--obname', help='Output basename, i.e. w/o extension', required=TRUE)

    parser$add_argument('--pv_col', help='Relevant p-value column', default="wi.ep", choices=c('we.ep', 'we.eBH', 'wi.ep', 'wi.eBH'))
    parser$add_argument('--pv_cutoff', help='P-value cutoff', default=0.05, type='double')
    parser$add_argument('--eff_cutoff', help='Effect size cutoff', default=0.5, type='double')

    parser$add_argument('--src', help='Source directory', default="src")
    return(parser)
}

# Plot ALDEx2 analysis results overview
plot_aldex_results <- function(tab, effect_cutoff=0.5, pvalue_cutoff=0.05, pvalue_col='wi.ep', title=''){
    # effect size flag
    tab$color_effect <- abs(tab$effect) > effect_cutoff
    # p-value flags
    tab$color_pvalue <- tab[,pvalue_col] < pvalue_cutoff
    # plotting params
    point_size <- 2.5
    point_shape<- 21
    point_alpha<- 0.7
    color_effect_values <- c('FALSE'='#333333', 'TRUE'='#FF3333')
    shape_pvalue_values <- c('FALSE'=21, 'TRUE'=24)
    # effec plot: within vs. between
    p1 <-
        ggplot(data=tab, aes(x=diff.win, y=diff.btw, fill=color_effect, shape=color_pvalue)) +
        geom_point(size=point_size, alpha=point_alpha, colour='white') +
        scale_shape_manual(values=shape_pvalue_values, guide=FALSE)+
        scale_fill_manual(values=color_effect_values, guide=FALSE) +
        geom_abline(slope=1,  intercept=0, linetype="dashed", colour="#333333") +
        geom_abline(slope=-1, intercept=0, linetype="dashed", colour="#333333") +
        ggplot2::xlim(min(0, min(tab$diff.win)), max(0, max(tab$diff.win))) +
        ggplot2::ylim(min(0, min(tab$diff.btw)), max(0, max(tab$diff.btw))) +
        xlab('Within group difference') +
        ylab('Between group difference') +
        labs(title=title, subtitle='Effect size plot') +
        # ggtitle('Effect size plot') +
        theme_bw()
    # bland-altman plot
    p2 <-
        ggplot(data=tab, aes(x=rab.all, y=diff.btw, fill=color_effect, shape=color_pvalue)) +
        geom_point(size=point_size, alpha=point_alpha, colour='white') +
        scale_shape_manual(values=shape_pvalue_values, guide=FALSE)+
        scale_fill_manual(values=color_effect_values, guide=FALSE) +
        xlab('Abundance (CLR)') +
        ylab('Between group difference') +
        labs(title='', subtitle='Bland-Altman plot') +
        # ggtitle('Bland-Altman plot') +
        theme_bw()
    # difference vs. p-value plot
    p3 <-
        ggplot(data=tab, aes_string(x='diff.btw', y=pvalue_col, fill='color_effect', shape='color_pvalue')) +
        geom_point(size=point_size, alpha=point_alpha, colour='white') +
        scale_shape_manual(values=shape_pvalue_values, guide=FALSE)+
        scale_fill_manual(values=color_effect_values, guide=FALSE) +
        scale_y_continuous(
            trans = log10_trans(),
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        geom_hline(yintercept=pvalue_cutoff, linetype="dotted") +
        xlab('Between group difference') +
        ylab('Log(p-value)') +
        labs(title='', subtitle='Difference vs. expected p-value') +
        # ggtitle('Difference vs. p-value') +
        theme_bw()
    # effect vs. p-value plot
    p4 <-
        ggplot(data=tab, aes_string(x='effect', y=pvalue_col, fill='color_effect', shape='color_pvalue')) +
        geom_point(size=point_size, alpha=point_alpha, colour='white') +
        scale_shape_manual(values=shape_pvalue_values, guide=FALSE)+
        scale_fill_manual(values=color_effect_values, guide=FALSE) +
        scale_y_continuous(
            trans = log10_trans(),
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        ) +
        geom_hline(yintercept=pvalue_cutoff, linetype="dotted") +
        geom_vline(xintercept=-effect_cutoff, linetype="dotted") +
        geom_vline(xintercept= effect_cutoff, linetype="dotted") +
        xlab('Effect') +
        ylab('Log(p-value)') +
        labs(title='', subtitle='Effect vs. expected p-value') +
        # ggtitle('Effect vs. p-value') +
        theme_bw()
    # all together
    p <- p1 + p2 + p3 + p4 + plot_layout(ncol=2)
    return(p)
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
ifiles <- c(args$otus, args$otu_tree, args$clr, args$sample_meta)
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

## DIFF
# CLR-transformation
flog.info('ALDEx2: CLR-transformation', name=logger)
set.seed(42)
otu_clr <- aldex.clr(
    reads=otus[,rownames(meta)], # OTUs x samples (!!!)
    conds=meta$Indication,
    denom="iqlr"
)

# Calculate effect sizes and differences between conditions
flog.info('ALDEx2: Calculate effect sizes', name=logger)
set.seed(23)
otu_eff <- aldex.effect(
    clr=otu_clr,
    conditions=meta$Indication,
    include.sample.summary=TRUE
)

# Calculate Welch's t-test and Wilcoxon test statistics
flog.info('ALDEx2: Calculate test statistics', name=logger)
set.seed(21)
otu_ttest <- aldex.ttest(
    clr=otu_clr,
    conditions=meta$Indication,
)

# Put together into one table
otu_diff <- data.frame(otu_eff, otu_ttest)
# change column order
otu_diff <- otu_diff[,c(
    "rab.win.Control", "rab.win.Sarcoidosis",
    "diff.btw", "diff.win", "effect", "overlap", "we.ep", "we.eBH", "wi.ep", "wi.eBH",
    # rest
    "rab.all",
    colnames(otu_diff)[grepl('^rab\\.sample\\.', colnames(otu_diff))]
)]
# sort w.r.t. p-value (ascending) and effect size (descending)
otu_diff <- otu_diff[with(otu_diff, order(otu_diff[,args$pv_col], -otu_diff[,'effect'])), ]
# OTU as column
otu_diff <- cbind(OTU=rownames(otu_diff), otu_diff)

# Save
write.table(x=otu_diff, file=sprintf('%s_aldex2.tsv', args$obname), sep='\t', row.names=FALSE)

## Plots
plots <- list()

# Overview plot
plots[['aldex2_ov']] <- plot_aldex_results(
    tab=otu_diff,
    effect_cutoff=args$eff_cutoff,
    pvalue_cutoff=args$pv_cutoff,
    pvalue_col=args$pv_col,
    title=''
)
pdf(sprintf('%s_aldex2_ov.pdf', args$obname), width=7, height=5)
print(plots[['aldex2_ov']])
dev.off()

## RDS
saveRDS(plots, file=sprintf('%s.rds', args$obname))
