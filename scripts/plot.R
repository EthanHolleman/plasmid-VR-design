library(stringr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(RColorBrewer)

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


read_name_and_id_num_from_parsed_RNA <- function(parsed.RNA.filepath){

    split <- unlist(strsplit(as.character(parsed.RNA.filepath), '/'))  
    id.num <- split[length(split)-2]
    name <-  split[length(split)-3]

    c(name, id.num)
}



merge_SPOT_RNA_predictions <- function(parsed.RNA.filepaths, df){

    parsed.rna.list <- list()
    for (i in 1:length(parsed.RNA.filepaths)){
        filepath <-  parsed.RNA.filepaths[[i]]
        name_id_num <- read_name_and_id_num_from_parsed_RNA(filepath)
        file.contents <- as.data.frame(read.table(filepath), sep='\t')
        colnames(file.contents) <- c(
            'structure_file', 'prop_unpaired', 'length', 
            'num_hairpins', 'prop_hairpin'
        )
        row <- c(
            name_id_num, file.contents[, 'prop_unpaired'], 
            file.contents[, 'prop_hairpin']
            )
        names(row) <- c('name', 'id_num', 'prop_unpaired', 'prop_hairpin')
        parsed.rna.list[[i]] <- row
    }

    parsed.rna.df <- as.data.frame(do.call(rbind, parsed.rna.list))
    df.merge <- merge(df, parsed.rna.df, by=c('name', 'id_num'))

    df.merge

}


read_expectation_files <- function(expectations.paths){

    expect.list <- list()
    for (i in 1:length(expectations.paths)){

        df <- as.data.frame(read.table(expectations.paths[[i]], sep='\t', header=T))
        expect.list[[i]] <- df

    }

    df <- as.data.frame(do.call(rbind, expect.list))
    df

}


distrabution_from_exectation_parameters <- function(expectation.df, nsamples=1000){

    print(expectation.df)
    dists.list <- list()
    for (i in 1:nrow(expectation.df)){
        expect.mean <- expectation.df[i, ]$mean
        expect.sd <- expectation.df[i, ]$sd
        samples <- rnorm(nsamples, m=as.numeric(expect.mean), sd=as.numeric(expect.sd))  # for now assume no neg values
        df <- data.frame(
            value=samples,
            mean=expect.mean,
            sd=expect.sd,
            length=expectation.df[i, ]$length,
            attribute=expectation.df[i, ]$attribute
        )
        dists.list[[i]] <- df
    }
    all.expect.df <- do.call(rbind, dists.list)

}


plot_deviation_from_expectation_metrics <- function(expectations.filepaths, df, i){

    expect.df <- read_expectation_files(expectations.filepaths)
    expect.dists <- distrabution_from_exectation_parameters(expect.df)
    expect.dists <- subset(expect.dists, length=df[i, 'length'])
    plots <- list()
    metrics <- unique(expect.df$attribute)
    colors <- brewer.pal(length(metrics), 'Dark2')
    for (j in 1:length(metrics)){
        metric <- metrics[[j]]
        var.region.metric_value <- df[i, metric]
        expect.dists.metric <- subset(expect.dists, attribute==metric)
        
        plot.metric <- ggplot(expect.dists.metric, aes(x=as.numeric(value))) + 
                        geom_density(color=colors[[j]], alpha=0.7) +
                        theme_pubr() + 
                        geom_vline(
                            xintercept=as.numeric(var.region.metric_value), 
                            linetype='dashed'
                        ) +
                         geom_vline(
                            xintercept=as.numeric(expect.dists.metric[, 'mean']),
                            linetype='dashed',
                            color=colors[[j]]
                        ) +
                        labs(x=metric) +
                        theme(
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank()
                        ) +
                        + theme(
                            axis.text.x = element_text(
                                angle = 90, vjust = 0.5, hjust=1
                                )
                            )
        plots[[j]] <- plot.metric
    }

    ggarrange(plotlist=plots)

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
    split[length(split)-2]

}

extract_vr_id_num_from_rlooper_filepath <- function(file.path){
    split <- unlist(strsplit(as.character(file.path), '/'))
    split[length(split)-1]

}


merge_rlooper_calculations <- function(df, rlooper.filepaths, attribute){

    # use the variable region names which act as IDs to get average local
    # energy calculations in the same order as the
    calc.filepaths <- list()
    for (i in 1:length(rlooper.filepaths)){
        name <- extract_vr_name_from_rlooper_filepath(rlooper.filepaths[[i]])
        id <- extract_vr_id_num_from_rlooper_filepath(rlooper.filepaths[[i]])
        mean.attribute <- mean(read_rlooper_wig_file(rlooper.filepaths[[i]])$value)
        calc.filepaths[[i]] <- c(name, as.numeric(id), rlooper.filepaths[[i]], mean.attribute)
        
    }

    calc.df <- as.data.frame(do.call(rbind, calc.filepaths))
    colnames(calc.df) <- c('name', 'id_num', paste(attribute, 'path', sep='_'), attribute)
    # merge rlooper filepaths into the dataframe
    df.merge <- merge(df, calc.df, by=c('name', 'id_num'))
    
    df.merge

}


read_rlooper_wig_file <- function(rlooper.filepath){

    df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(df) <- c('value')
    df$position <- 1:nrow(df)
    df$value <- as.numeric(df$value)

    df
}


plot_rlooper_calcs <- function(df, i, rlooper_attributes){
    
    plots <- list()
    colors <- brewer.pal(length(rlooper_attributes), 'Dark2')
    for (i in 1:length(rlooper_attributes)){
        attribute <- rlooper_attributes[[i]]
        message(attribute)
        cn <- paste(attribute, 'path', sep='_')
        message(cn)
        file.path <- df[i, c(cn)]
        message(file.path)

        df.seq <- read_rlooper_wig_file(file.path)
        seq <- df[i, ]$Sequence
        p <- ggplot(df.seq, aes(x=position, y=value)) + 
                geom_point(color=colors[[i]]) + 
                geom_line(color=colors[[i]]) +
                theme_pubr() + 
                labs(x='Nucleotide position', y=attribute)
        plots[[i]] <- p
    }

    ggarrange(plotlist=plots, nrow=length(plots), ncol=1)
    
}



main <- function(){
    save.image('plot.image.RData')
    print('================================')
    print('Reading / parsing inputs')
    print('================================')
    input.path <- as.character(snakemake@input['variable_regions'])
    output.path <- as.character(snakemake@output)
    #save.image('plot.RData')
    df <- read_variable_region_tsv(input.path)
    df <- merge_rlooper_calculations(
        df, 
        snakemake@input['rlooper_lae']$rlooper_lae,
        'local_average_energy'
    )
    df <- merge_rlooper_calculations(
        df, 
        snakemake@input['rlooper_bprob']$rlooper_bprob,
        'bp_prob'
    )
    df <- merge_SPOT_RNA_predictions(snakemake@input['parsed_RNA']$parsed_RNA, df)

    expectation.filepaths <- unique(snakemake@input['expectation_files']$expectation_files)

    print('================================')
    print('Making plots')
    print('================================')
    pdf(output.path, width=18, height=14)
    for (i in 1:nrow(df)){
        skew_content_plot <- skew_content_plot_from_df_row(df, i, WINDOW_SIZE)
        vr_table <- variable_region_table(df, i)
        diffs_barplot <- calculated_windowed_skew_content_vs_params(df, i, WINDOW_SIZE)
        nuc_props <- nucleotide_proportions(df, i)
        sequence <- plot_seq(df, i)
        clustering <- plot_distance_to_next_same_nucleotide(df, i)
        rlooper <- plot_rlooper_calcs(df, i, c('local_average_energy', 'bp_prob'))
        expectations <- plot_deviation_from_expectation_metrics(expectation.filepaths, df, i)

        skew_and_rlooper <- ggarrange(
            skew_content_plot, rlooper, nrow=2, ncol=1, heights=c(1, 0.75)
        )

        seq_and_table <- ggarrange(
            vr_table,
            expectations,
            nuc_props,
            ggarrange(clustering, sequence, nrow=2, ncol=1, heights=c(1, 0.5)),
            nrow=1, ncol=4
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