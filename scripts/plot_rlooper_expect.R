# read all files at once?
library(ggplot2)
library(ggpubr)
library(ggridges)
library(RColorBrewer)

read_rlooper_wig <- function(rlooper.filepath){

    print(rlooper.filepath)
    df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(df) <- c('value')
    c(mean(df$value), nrow(df))

}


read_all_wig_means <- function(wig.paths){

    means <- list()
    for (i in 1:length(wig.paths)){

        means[[i]] <- read_rlooper_wig(wig.paths[[i]])

    }
    means.df <- as.data.frame(do.call(rbind, means))
    colnames(means.df) <- c('mean', 'length')
    means.df

}


plot_means <- function(means.df, stat_name, palette='Blues'){

    num_seqs <- nrow(means.df)
    ggplot(means.df, aes(x=mean, y=as.factor(length), fill=length)) +
        theme_pubr() + geom_density_ridges(alpha=0.7) + 
        labs(title=paste('Mean', stat_name), x=stat_name, y='Sequence length') +
        theme(legend.position = "none")
    
}

main <- function(){

    ale.files <- snakemake@input$ale
    bpprob.files <- snakemake@input$bpprob
    output.path <- as.character(snakemake@output)

    ale.means.df <- read_all_wig_means(ale.files)
    plot.ale <- plot_means(
        ale.means.df, 'Average local energy (G)'
        )

    bpprob.means.df <- read_all_wig_means(bpprob.files)
    plot.bprpob <- plot_means(
        bpprob.means.df, 'R-loop probability',
        palette='Greens'
        )

    all.plots <- ggarrange(plot.ale, plot.bprpob, nrow=2, ncol=1)
    ggsave(output.path, all.plots, dpi=500)

}

if (! interactive()){
    main()
}