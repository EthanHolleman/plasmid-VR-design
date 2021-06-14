library(stringr)
library(ggplot2)
library(ggpubr)

skew <- function(a_count, b_count, interval){

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

seq_skew_content_sliding_window <- function(seq, window_size=10){

    windows <- list()
    k <- 1
    for (i in 1:(nchar(seq)-window_size))

    {

        window_start = i
        window_end = i + window_size
        window_seq = substr(seq, window_start, window_end)

        nuc_counts <- get_nucleotide_counts(window_seq)
        gc_skew <- skew(nuc_counts[['G']], nuc_counts[['C']])
        at_skew <- skew(nuc_counts[['A']], nuc_counts[['T']])

        gc_content <- content(nuc_counts[['G']], nuc_counts[['C']], nchar(window_seq))
        at_content <- content(nuc_counts[['A']], nuc_counts[['T']], nchar(window_seq))

        print(i)
        windows[[k]] <- c(i, gc_skew, 'GC_skew')
        windows[[k+1]] <- c(i, gc_content, 'GC_content')
        windows[[k+2]] <- c(i, gc_content, 'GC_content')
        windows[[k+3]] <- c(i, at_skew, 'AT_skew')
        windows[[k+4]] <- c(i, at_content, 'AT_content')
        k = k + 5

    }

    df <- as.data.frame(do.call(rbind, windows))
    colnames(df) <- c('window_number', 'value', 'metric')
    df


}


plot_skew_content_windows <- function(windows.df){

    plot <- ggplot(windows.df, aes(x=as.numeric(window_number), y=as.numeric(value), color=metric)) +
            geom_point() + geom_line() + theme_pubr() + scale_color_brewer(palette='Dark2')
    ggsave('test.png', plot)
        

}

seq='ATTTGTGTACCACAGTGTGTACACATGTGTGACACATGTACAACTGGTGTGAACCAACACAGTGTGACACGTGTGAC'
df <- seq_skew_content_sliding_window(seq)
plot_skew_content_windows(df)