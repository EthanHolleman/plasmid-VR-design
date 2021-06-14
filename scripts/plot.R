library(stringr)
library(ggplot2)
library(ggpubr)


read_variable_region_tsv <- function(file.path){

    as.data.frame(read.table(file.path, sep='\t', header=T))

}



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

## Calculate skew and content metrics for a given character vector
## using a specific window size.
seq_skew_content_sliding_window <- function(seq, window_size=30){
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


plot_skew_content_windows <- function(windows.df, region_name){

    plot <- ggplot(windows.df, aes(x=as.numeric(window_number), y=as.numeric(value), color=metric)) +
            geom_point() + geom_line() + theme_pubr() + scale_color_brewer(palette='Dark2') +
            labs(title=region_name)
    plot
        

}

## generate skew and content plots from a variable region dataframe and
## the row number of a variable region.
skew_content_plot_from_df_row <- function(df, row_num, window_size=30){

    row <- df[row_num, ]
    seq <- row$Sequence
    windows.df <- seq_skew_content_sliding_window(seq, window_size)
    plot_skew_content_windows(windows.df, row$name)

}

variable_region_table <- function(df, row_num){

    text_table <- df[row_num, c('name', 'GC_content', 'GC_skew', 'AT_content', 'AT_skew')]
    ggtexttable(text_table, rows=NULL)

}


calculated_windowed_skew_content_vs_params <- function(df, row_num, window_size=30){
    
    row <- df[row_num, ]
    seq <- row$Sequence
    windows.df <- seq_skew_content_sliding_window(seq, window_size)
    
    diff_vals <- list()
    for (i in 1:nrow(windows.df)){
        diff_vals[[i]] <- as.numeric(windows.df[i, ]$value) - as.numeric(df[row_num, windows.df[i, ]$metric])
    }

    windows.df$diff_vals <- unlist(diff_vals)

    ggplot(windows.df, aes(x=as.factor(as.numeric(window_number)), y=diff_vals), fill=metric) +
           geom_bar(stat='identity') +
           scale_fill_brewer(palette='Dark2') + theme_pubr() + facet_wrap(~metric) +
           theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
           scale_x_discrete(breaks = seq(1, max(as.numeric(windows.df$window_number)), by = 10))

}


main <- function(){

    input.path <- as.character(snakemake@input)
    output.path <- as.character(snakemake@output)
    #save.image('plot.RData')
    df <- read_variable_region_tsv(input.path)
    
    pdf(output.path, width=10, height=10)
    for (i in 1:nrow(df)){
        skew_content_plot <- skew_content_plot_from_df_row(df, i)
        vr_table <- variable_region_table(df, i)
        diffs_barplot <- calculated_windowed_skew_content_vs_params(df, i)
        
        aranged <- ggarrange(skew_content_plot, vr_table, diffs_barplot, nrow=3, ncol=1)
        print(aranged)
    }
    dev.off()

}


if (!interactive()){

    main()

}