library(GGally)
library(viridis)
library(ggpubr)
library(RColorBrewer)
# make a parallel coordinates plot (spagetti plot) that shows all candidate
# sequences for a particular set of variable region parameters.


read_seq_rankings <- function(file.path){

    as.data.frame(read.table(file.path, sep='\t', header=TRUE))

}

get_metric_rankings <- function(rank.df){
    ranks <- c()
    for (each_col in colnames(rank.df)){
        is_rank <- grepl('rank', each_col, fixed=TRUE)
        if (is_rank){
            ranks <- c(ranks, each_col)
        }


    }

    ranks

}


make_spagetti <- function(rank.df){

    ranks <- get_metric_rankings(rank.df)
    rank.only.df <- rank.df[, ranks]
    # rank.only.df$id_num <- rank.df$id_num
    ggparcoord(
        rank.only.df,
        columns=c(1:length(ranks)),
        groupColumn=5,
        alphaLines=0.4,
    ) +
    scale_color_viridis() +
    theme_pubr() + labs(y='Sequence rank', x='Metric', title='Sequence rankings')

}

metrics <- function(rank.df){
    # get the fieldnames of the actual metric values not their rankings
    metrics.ranks <- get_metric_rankings(rank.df)
    ranks <- c()
    for (i in 1:length(metrics.ranks)){
        ranks <- c(ranks, gsub('_ranks', '', metrics.ranks[[i]]))  # remove _ranks prefix
    }
    ranks

}


metric_distribution_plot <- function(rank.df, metric, color){

    # get top metric
    top_rank = subset(rank.df, overall_rank==0)
    top_rank.metric.value = top_rank[, c(metric)]
    ggplot(rank.df, aes_string(x=metric)) + geom_density(fill=color, alpha=0.7) +
    geom_vline(xintercept=top_rank.metric.value,  linetype="dashed") + theme_pubr() + labs(x=metric)

}

all_metric_distribution_plots <- function(rank.df){
    metrics.names <- metrics(rank.df)
    print(metrics.names)
    colors <- brewer.pal(n=length(metrics.names), "Dark2")
    plots <- list()
    for (i in 1:length(metrics.names)){
        plots[[i]] <- metric_distribution_plot(
            rank.df, metrics.names[[i]], colors[[i]])
    }
    ggarrange(plotlist=plots)


}



main <- function(){

    input.ranks <- as.character(snakemake@input)
    output.plot <- as.character(snakemake@output)
    rank.df <- read_seq_rankings(input.ranks)
    plot <- make_spagetti(rank.df)
    dists <- all_metric_distribution_plots(rank.df)
    main.plot <- ggarrange(plot, dists)
    ggsave(output.plot, main.plot, dpi=300, width=24, height=12)


}

if (!interactive()){
    main()
}