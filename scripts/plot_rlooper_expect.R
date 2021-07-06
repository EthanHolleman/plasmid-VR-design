# read all files at once?
library(ggplot2)
library(ggpubr)
library(ggridges)
library(RColorBrewer)

read_rlooper_wig <- function(rlooper.filepath){


    df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(df) <- c('value')
    c(mean(df$value), sd(df$value), nrow(df))

}


read_all_wig_means <- function(wig.paths){

    means <- list()
    for (i in 1:length(wig.paths)){

        means[[i]] <- read_rlooper_wig(wig.paths[[i]])

    }
    means.df <- as.data.frame(do.call(rbind, means))
    colnames(means.df) <- c('mean', 'sd', 'length')
    means.df$se <- means.df$sd / (sqrt(means.df$length))
    means.df

}


plot_means <- function(means.df, stat_name, palette='Blues'){

    num_seqs <- nrow(means.df)
    dists <- ggplot(means.df, aes(x=mean, y=as.factor(length), fill=length)) +
        theme_pubr() + geom_density_ridges(alpha=0.7) + 
        labs(title=paste('Mean', stat_name), x=paste('Mean', stat_name), y='Sequence length') +
        theme(legend.position = "none")
    

    means <- ggplot(means.df, aes(x=as.numeric(length), y=mean)) +
             geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
             geom_line(color='dodgerblue') + geom_point(color='black') +
             theme_pubr() + labs(x='Mean local average energy', y='Sequence length')
    
    dists
    
}

agg_variable_mean_sd <- function(big.df, var_name, group_name){

    means <- aggregate(big.df[, 'mean'], list(big.df[, group_name]), mean)
    sds <- aggregate(big.df[, 'mean'], list(big.df[, group_name]), sd)
    merge.df <- merge(means, sds, by='Group.1')
    merge.df$attribute <- var_name
    colnames(merge.df) <- c(group_name, 'mean', 'sd', 'attribute')
    merge.df

}


write_expectation_table <- function(means.prob, means.energy, output.path){

    means.energy <- agg_variable_mean_sd(means.energy, 'local_average_energy', 'length')
    means.prob <- agg_variable_mean_sd(means.prob, 'bp_prob', 'length')
    all.means <- rbind(means.energy, means.prob)
    write.table(all.means, file=output.path, sep='\t', quote=FALSE, row.names=FALSE)

}


main <- function(){

    ale.files <- snakemake@input$ale
    bpprob.files <- snakemake@input$bpprob
    output.path.plot <- as.character(snakemake@output$plot)
    output.path.expect <- as.character(snakemake@output$expect)
    
    ale.means.df <- read_all_wig_means(ale.files)
    plot.ale <- plot_means(
        ale.means.df, 'local average energy'
        )

    bpprob.means.df <- read_all_wig_means(bpprob.files)
    plot.bprpob <- plot_means(
        bpprob.means.df, 'R-loop probability',
        palette='Greens'
        )

    all.plots <- ggarrange(plot.ale, plot.bprpob, nrow=2, ncol=1)
    write_expectation_table(bpprob.means.df, ale.means.df, output.path.expect)
    ggsave(output.path.plot, all.plots, dpi=500)

}

if (! interactive()){
    main()
}