library(stringr)
library(ggplot2)
library(ggpubr)
library(ggridges)

WINDOW_SIZE <- 30


read_variable_region_tsv <- function(file.path){

    as.data.frame(read.table(file.path, sep='\t', header=T))

}

## Cacluate skew of nucleotide B over A
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
            labs(title=region_name) + labs(x='Window Number', y='Metric value')
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

    text_table_reg <- df[row_num, c('GC_content', 'GC_skew', 'AT_content', 'AT_skew')]
    text_table_clust <- df[row_num, c('Cluster_length', 'Clustered_nucleotide', 'Clustering.method')]
    
    row <- df[row_num, ]
    seq <- row$Sequence
    nuc_counts <- get_nucleotide_counts(seq)

    gc_skew <- skew(nuc_counts[['G']], nuc_counts[['C']])
    at_skew <- skew(nuc_counts[['A']], nuc_counts[['T']])
    gc_content <- content(nuc_counts[['G']], nuc_counts[['C']], nchar(seq))
    at_content <- content(nuc_counts[['A']], nuc_counts[['T']], nchar(seq))

    calc_metrics_table <- data.frame(
        gc_skew=gc_skew, at_skew=at_skew, gc_content=gc_content, at_content=at_content
    )

    t1 <- ggtexttable(text_table_reg, rows=NULL, theme = ttheme("light"))
    t1 <- tab_add_footnote(t1, 'User defined parameters.')
    t2 <- ggtexttable(text_table_clust, rows=NULL, theme = ttheme("light"))
    t3 <- ggtexttable(calc_metrics_table, rows=NULL, theme = ttheme("light"))
    t3 <- tab_add_footnote(t3, 'Values calculated from sequence.')
    ggarrange(t1, t2, t3, nrow=3, ncol=1)

}

distance_to_next_same_nucleotide <- function(seq, nuc_index){

    nuc <- substr(seq, nuc_index, nuc_index)
    dist <- nchar(seq)
    k <- 1
    right_index <- nuc_index
    left_index <- nuc_index
    while (right_index <= nchar(seq) | left_index >= 1){
        k <- k + 1
        left_index <- left_index - 1
        right_index <- right_index + 1
        for (index in c(left_index, right_index)){
            if (index > 0 & index < nchar(seq)){
                index_nuc <- substr(seq, index, index)
                if (index_nuc == nuc){
                    dist <- abs(abs(index) - nuc_index)
                    return(dist)
                }
            }
        }

        if (k > nchar(seq)){
            
            break  # emergency stop
        }
    }
    dist

}

plot_distance_to_next_same_nucleotide <- function(df, row_num){

    row <- df[row_num, ]
    seq <- row$Sequence

    dists <- list()
    for (nuc_index in 1:nchar(seq)){
        nuc <- substr(seq, nuc_index, nuc_index)
        dists[[nuc_index]] <- c(
            distance_to_next_same_nucleotide(seq, nuc_index),
            nuc
            )
    }
    dist.df <- as.data.frame(do.call(rbind, dists))
    colnames(dist.df) <- c('Distance', 'Nucleotide')
    plot <- ggplot(dist.df, aes(y=Nucleotide, x=as.numeric(Distance), fill=Nucleotide)) + 
           geom_density_ridges(alpha=0.7) + theme_pubr() + 
            scale_fill_brewer(palette='Dark2') + theme(legend.position = "none") +
            labs(x='Distance to nucleotide of same species')
    plot


}


nucleotide_proportions <- function(df, row_num){

    row <- df[row_num, ]
    seq <- row$Sequence
    nuc_counts <- get_nucleotide_counts(seq)
    df <- data.frame(nucleotide=names(nuc_counts), count=unlist(nuc_counts))
    ggplot(df, aes(x=nucleotide, y=count, fill=nucleotide)) +
           geom_bar(stat='identity', color='black') + theme_pubr() + 
           scale_color_brewer(palette='Dark2') +  theme(legend.position = "none") +
           labs(x='Nucleotide', y='Count')
    
}

## Calculate difference between sequence metric in each sliding window
## and the value of the metric that was originally specified
calculated_windowed_skew_content_vs_params <- function(df, row_num, window_size=30){
    
    row <- df[row_num, ]
    seq <- row$Sequence
    windows.df <- seq_skew_content_sliding_window(seq, window_size)
    
    diff_vals <- list()
    for (i in 1:nrow(windows.df)){
        diff_vals[[i]] <- as.numeric(windows.df[i, ]$value) - as.numeric(df[row_num, windows.df[i, ]$metric])
    }

    windows.df$diff_vals <- unlist(diff_vals)

    barplot <- ggplot(windows.df, aes(x=as.factor(as.numeric(window_number)), y=diff_vals), fill=metric) +
           geom_bar(stat='identity') +
           scale_fill_brewer(palette='Dark2') + theme_pubr() + facet_wrap(~metric) +
           theme(legend.position = "none") +
           scale_x_discrete(breaks = seq(1, max(as.numeric(windows.df$window_number)), by = 25)) +
           labs(x='Window Number', y='Difference from parameter')

    boxplot <- ggplot(windows.df, aes(y=diff_vals, x=metric, fill=metric)) +
                geom_boxplot() + theme_pubr() + theme(legend.position = "none") +
                scale_fill_brewer(palette='Dark2') +
                theme(axis.text.x = element_text(angle = 45, hjust=1)) +
                labs(x='Metric', y='Difference from parameter')
    
    ggarrange(barplot, boxplot, nrow=1, ncol=2, widths = c(1, 0.5))


}

## Format nucleotide sequence for printing to
## plot
format_sequence_string <- function(df, i){

    wrap = 40
    row <- df[i, ]
    name <- paste('>', row$name, sep='')
    seq_str <- ''
    for (i in seq(1, nchar(row$Sequence)-wrap, wrap)){
        seq_str <- paste(seq_str, substr(row$Sequence, i, i+wrap), sep='\n')
    }

    formated <- seq_str
    formated

}

## Plot the actual nucleotide sequence as a ggparagraph
plot_seq <- function(df, i){

    formated_seq <- format_sequence_string(df, i)
    ggparagraph(formated_seq, face = "bold", size = 14, color = "black")
}


# plot_text <- function(df, i, window_size){

#     row <- df[i, ]
#     text <- paste(
#         'Sequence metrics for variable region',
#         row$name,
#         '.',
#         'A: Lineplot showing GC AT skew and content over',
#         window_size,
#         'nulceotide intervals.'

#     )
#     ggparagraph(text = text, face = "italic", size = 11, color = "black")

# }

extract_vr_name_from_rlooper_filepath <- function(file.path){
    split <- unlist(strsplit(as.character(file.path), '/'))
    split[length(split)-1]

}


merge_rlooper_calculations <- function(df, rlooper.filepaths){

    # use the variable region names which act as IDs to get average local
    # energy calculations in the same order as the
    calc.filepaths <- list()
    for (i in 1:length(rlooper.filepaths)){
        name <- extract_vr_name_from_rlooper_filepath(rlooper.filepaths[[i]])
        calc.filepaths[[i]] <- c(name, rlooper.filepaths[[i]])
    }
    calc.df <- as.data.frame(do.call(rbind, calc.filepaths))
    colnames(calc.df) <- c('name', 'rlooper_filepath')
    # merge rlooper filepaths into the dataframe
    df.merge <- merge(df, calc.df, by='name')
    
    df.merge

}


read_average_local_energy_file <- function(rlooper.filepath){

    print('Reacing File')
    print(rlooper.filepath)
    ale.df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(ale.df) <- c('ale')
    ale.df$position <- 1:nrow(ale.df)
    ale.df$ale <- as.numeric(ale.df$ale)

    ale.df
}


plot_average_local_energy <- function(df, i){
    
    file.path <- df[i, ]$rlooper_filepath
    ale.df <- read_average_local_energy_file(file.path)
    ggplot(ale.df, aes(x=position, y=ale)) + geom_point(color='#e32b17') + geom_line(color='#e32b17') +
            theme_pubr() + scale_x_discrete(breaks = seq(1, max(ale.df$position), by = 10)) +
            labs(x='Nucleotide position', y='Average G rlooper')

}



main <- function(){

    input.path <- as.character(snakemake@input['variable_regions'])
    output.path <- as.character(snakemake@output)
    #save.image('plot.RData')
    df <- read_variable_region_tsv(input.path)
    df <- merge_rlooper_calculations(df, 
    snakemake@input['rlooper_calculations']$rlooper_calculations)
    
    pdf(output.path, width=18, height=14)
    for (i in 1:nrow(df)){
        skew_content_plot <- skew_content_plot_from_df_row(df, i, WINDOW_SIZE)
        vr_table <- variable_region_table(df, i)
        diffs_barplot <- calculated_windowed_skew_content_vs_params(df, i, WINDOW_SIZE)
        nuc_props <- nucleotide_proportions(df, i)
        sequence <- plot_seq(df, i)
        clustering <- plot_distance_to_next_same_nucleotide(df, i)
        rlooper <- plot_average_local_energy(df, i)

        skew_and_rlooper <- ggarrange(
            skew_content_plot, rlooper, nrow=2, ncol=1, heights=c(1, 0.75)
        )

        seq_and_table <- ggarrange(
            vr_table,
            nuc_props,
            ggarrange(clustering, sequence, nrow=2, ncol=1, heights=c(1, 0.5)),
            nrow=1, ncol=3
        )

        
        #figure_text <- plot_text(df, i, WINDOW_SIZE)

        aranged <- ggarrange(
            skew_and_rlooper, seq_and_table, diffs_barplot,
            nrow=3, ncol=1)
        print(aranged)
    }
    dev.off()

}


if (!interactive()){

    main()

}