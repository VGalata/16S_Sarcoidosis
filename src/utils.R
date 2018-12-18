#!/usr/bin/Rscript

##################################################
# LOGGING
##################################################
add_logger <- function(logger_file, logger_name=NULL, logger_level='INFO'){
    # Levels: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
    if(is.null(logger_name)){
        logger_name <- sprintf('LOGGER_%s', logger_level)
    }
    switch(logger_level,
        'INFO'  = { flog.logger(name=logger_name, INFO,  appender=appender.file(logger_file)) },
        'TRACE' = { flog.logger(name=logger_name, TRACE, appender=appender.file(logger_file)) },
        'DEBUG' = { flog.logger(name=logger_name, DEBUG, appender=appender.file(logger_file)) },
        'WARN'  = { flog.logger(name=logger_name, WARN,  appender=appender.file(logger_file)) },
        'ERROR' = { flog.logger(name=logger_name, ERROR, appender=appender.file(logger_file)) },
        'FATAL' = { flog.logger(name=logger_name, FATAL, appender=appender.file(logger_file)) },
        { flog.logger(name=logger_name, INFO, appender=appender.file(logger_file)) }
    )
    return(logger_name)
}

log_session_info <- function(logname){
    # log session information
    flog.info(
        sprintf(
            'SESSION INFO START\n%s\nSESSION INFO END',
            paste(capture.output(sessionInfo()), collapse='\n')
        ),
        name=logname
    )
}

log_args <- function(args, logname){
    # log given arguments
    flog.info(
        sprintf(
            'ARGS\n%s',
            paste(sapply(sort(names(args)), function(x){ sprintf('%s: %s', x, as.character(args[x])) }), collapse='\n')
        ),
        name=logname
    )
}

log_tab <- function(tab, logname, tab_name="Table", nrows=5, ncols=5){
    ncols <- min(ncols, ncol(tab))
    nrows <- min(nrows, nrow(tab))
    flog.info(
        sprintf(
            'Table "%s": %d x %d:\n%s',
            tab_name, nrow(tab), ncol(tab),
            capture_print(head(tab[,1:ncols,drop=FALSE], nrows))
        ),
        name=logname
    )
}

capture_print <- function(obj){
    # capture print output (for logging)
    return(paste(capture.output(print(obj)), collapse='\n'))
}

##################################################
# DATA
##################################################
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
ranks_ampvis <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
annot_map <- list(
    'Group'=c('1'='Sarcoidosis', '2'='Control'),
    'Sex'=c('1'='female', '2'='male'),
    'Smoker'=c('1'='current', '2'='former', '3'='never')
)

#' File modification timestamp as string
fmtime_str <- function(fname){
    sprintf('%s', file.mtime(fname))
}

read_meta <- function(meta_file, cols=NULL, id_map_file=NULL){
    if(is.null(cols) | length(cols)==0){
        cols <- c('SampleID', 'Indication')
    }
    df <- read.csv(file=meta_file, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    df <- df[,cols,drop=FALSE]
    df <- df[!base::duplicated(df),,drop=FALSE]
    if(!is.null(id_map_file)){
        id_map <- read.csv(file=id_map_file, sep='\t', header=FALSE, stringsAsFactors=FALSE)
        colnames(id_map) <- c('OLD', 'NEW')
        testit::assert( length(intersect(df$SampleID, id_map$OLD)) == length(unique(df$SampleID)) )
        df$SampleIDNew <- sapply(df$SampleID, function(x){ id_map[id_map$OLD == x, 'NEW'] })
    }
    rownames(df) <- df$SampleID
    return(df)
}

read_otus <- function(otus_file, flags=c('kept')){
    df <- read.csv(file=otus_file, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
    # filter by flag
    if(!is.null(flags) & length(flags) > 0 & 'FLAG' %in% colnames(df)){
        df <- df[sapply(df$FLAG, function(x){ x %in% flags }),]
    }
    return(df)
}

rm_tips_from_tree <- function(phy, tip_names, logger=''){
    return(
        ape::drop.tip(
            phy=phy,
            tip=tip_names,
            trim.internal=TRUE,
            subtree=FALSE,
            collapse.singles=TRUE,
            interactive=FALSE
        )
    )
}

read_annot <- function(annot_file, cols=NULL, convert_cols=TRUE){
    annot <- read.csv(file=annot_file, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA', '?'))
    # columns sub-set
    if(!is.null(cols)){
        annot <- annot[,cols,drop=FALSE]
    }
    # convert values in some columns
    if(convert_cols){
        for(col in names(annot_map)){
            if(col %in% colnames(annot)){
                annot[col] <- factor(
                    sapply(annot[col], function(x){ ifelse(is.na(x), NA, annot_map[[col]][x]) }),
                    levels=annot_map[[col]]
                )
            }
        }
    }
    return(annot)
}

merge_meta_annot <- function(meta, annot){
    testit::assert('F-Index' %in% colnames(meta))
    testit::assert('R-Index' %in% colnames(meta))
    testit::assert('F Index' %in% colnames(annot))
    testit::assert('R Index' %in% colnames(annot))
    row_map <- sapply(1:nrow(meta), function(x){
        x_match <- which( annot['F Index'] == meta[x,'F-Index'] &  annot['R Index'] == meta[x,'R-Index'])
        testit::assert(sprintf('No annot match for %s', paste(meta[x,], collapse=';')), length(x_match) == 1)
        return(x_match)
    })
    testit::assert(length(intersect(colnames(meta), colnames(annot))) == 0)
    meta <- cbind(meta, annot[row_map,setdiff(colnames(annot), c('F Index', 'R Index')),drop=FALSE])
    return(meta)
}

##################################################
# Transformations
##################################################
#' Preprocessing and CLR transformation
#' @param tab: Table counts (samples x features)
#' @return CLR transformed values (samples x features)
# Reference: https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot
# 1) impute zeros: cmultRepl(X,  label=0, method="CZM")
#       - ‘output="prop"’ by default
#       - X is samples x features
# 2) apply to each feature the CLR-transformation: apply(X, 2, function(x){log(x) - mean(log(x))})
#       - X is samples x features
#       - log(x) - mean(log(x)) = compositions::clr(x)
counts2clr <- function(tab){
    # rm zero features
    tab <- tab[,apply(tab, 2, function(x){ any(x > 0) })]
    flog.info(sprintf('RM zero features: %d x %d', nrow(tab), ncol(tab)), name=logger)
    log_tab(tab, logger)
    # impute zeros and calculate proportions/relative abundance
    tab <- zCompositions::cmultRepl(tab, method="CZM", label=0, output="prop")
    flog.info(sprintf('Impute zeros & calc. prop.: %d x %d', nrow(tab), ncol(tab)), name=logger)
    log_tab(tab, logger)
    # clr transformation
    tab <- t(apply(tab, 1, compositions::clr)) # samples x features
    flog.info(sprintf('CLR transformation: %d x %d', nrow(tab), ncol(tab)), name=logger)
    log_tab(tab, logger)
    return(tab)
}

##################################################
# PLOTS
##################################################
flag_colors <- c(
    'contamination'='#f75f55', # red
    'controls'='#619CFF', # blue
    'counts'='#ffd11a', # yellow
    'contamination;controls'='#d9b3ff', # purple
    'contamination;count'='#ff9933', # orange
    'controls;count'='#1aff66', # green
    'contamination;controls;count'='white',
    'kept'='#666666'
)

meta_colors <- list(
    Indication=c(Control='#91b6d4', Sarcoidosis='#ff8080'),
    Smoker=c(never='#80ffaa', former='#ffd633', current='#ffa64d')
)

#' Create list of comparisons for ggsignif::geom_signif()
#' @param v: Vector of group labels (str/factors; (un)sorted, (non)unique)
#' @return list with vectors of length 2 (all pairs of given unique labels)
get_comparisons_for_ggsignif <- function(v){
    # values -> str -> unique -> sort
    v <- sort(unique(as.character(v)))
    # all pairwise combinations as a list
    comparisons <- combn(v, 2, simplify=FALSE)
    return(comparisons)
}
