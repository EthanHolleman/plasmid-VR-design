library(stringr)
library(ggplot2)
library(ggpubr)

skew <- function(a_count, b_count interval){

    (b_count - a_count) / (a_count + b_count)

}

content <- function(a_count, b_count, seq_len){

    (a_count + b_count) / seq_len

}

get_nucleotide_counts <- function(seq){

    counts <- list()
    for (nuc in c('A', 'T', 'G', 'C')){
        counts[[nuc]] <- str_count(seq, nuc)
    }
    counts

}

seq_skew_content_sliding_window <- function(seq, window_size=25){

    windows <- list()

    for (i in seq(1, length(seq)-window_size, window_size){

        window_start = i
        window_end = i + window_size
        window_seq = substr(seq, window_start, window_end)

        nuc_counts <- get_nucleotide_counts(window_seq)
        
        gc_skew <- skew(nuc_counts['G'], nuc_counts['C'])
        at_skew <- skew(nuc_counts['A'], nuc_counts['T'])

        gc_content <- content(nuc_counts['G'], nuc_counts['C'], length(window_seq))
        at_content <- content(nuc_counts['A'], nuc_counts['T'], length(window_seq))

        windows[[i]] <- c(gc_skew, at_skew, gc_content, at_content)

    }

    df <- as.data.frame(do.call(rbind, windows))
    colnames(df) <- 'gc_skew', 'at_skew', 'gc_content', 'at_content'
    df


}