library(GGally)
library(viridis)
library(ggpubr)

# make a parallel coordinates plot (spagetti plot) that shows all candidate
# sequences for a particular set of variable region parameters.


read_seq_rankings <- function(file.path){

    as.data.frame(read.table(file.path, sep='\t'))

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
    rank.only.df <- ranks.df[, ranks]
    rank.only.df$id_num <- rank.df$id_num
    ggparcoord(
        rank.only.df
        columns=1:length(ranks),
        groupColumn=nrow(rank.only.df),
        order='anyClass',
        showPoints=TRUE,
        alphaLines=0.4,
    ) +
    scale_color_viridis(discrete=TRUE) +
    theme_pubr() + labs(y='Sequence rank', x='Metric', title='Sequence rankings')

}

main <- function(){

    input.ranks <- as.character(snakemake@input)
    output.plot <- as.character(snakemake@output)

    rank.df <- read_seq_rankings(input.ranks)
    plot <- make_spagetti(rank.df)

    ggsave(output.plot, plot, dpi=300)


}

if (!interactive()){
    main()
}