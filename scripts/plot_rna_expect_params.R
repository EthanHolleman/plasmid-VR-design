library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)


read_all_agg_metric_files <- function(file.path.list){

    dfs <- list()
    for (i in 1:length(file.path.list)){

        file.path <- as.character(file.path.list[[i]])
        dfs[[i]] <- as.data.frame(read.table(file.path, sep='\t'))

    }

    do.call(rbind, dfs)

}


grid_plots <- function(big.df){

    # at this point need df with only GC skew content and RNA metrics
    big.df.essential <- big.df[, c('GC_skew', 'GC_content', 'prop_unpaired', 'prop_hairpin')]

    big.df.agg <- aggregate(.~GC_skew+GC_content, big.df.essential, mean)
    colnames(big.df.agg) <-  c('GC_skew', 'GC_content', 'mean_prop_unpaired', 'mean_prop_hairpin')

    prop_hairpin <- ggplot(big.df, aes(GC_skew, GC_content, fill=mean_prop_hairpin)) +
            geom_tile() + theme_pubr() + scale_fill_viridis(discrete=FALSE) +
            labs(x='GC skew', y='GC content') +
             theme(legend.key.size = unit(2, 'cm')) +
            theme(text = element_text(size=24))

    prop_unpaired <- ggplot(big.df, aes(GC_skew, GC_content, fill=mean_prop_unpaired)) +
            geom_tile() + theme_pubr() + scale_fill_viridis(discrete=FALSE) +
            labs(x='GC skew', y='GC content') +
             theme(legend.key.size = unit(2, 'cm')) +
            theme(text = element_text(size=24))
    
    ggarrange(prop_hairpin, prop_unpaired, nrow=1, ncol=2)


}

main() <- function(){


    agg_metrics.files <- snakemake@input
    big.df <- read_all_agg_metric_files(agg_metrics.files)
    plots.grid <- grid_plots(big.df)

    output.path <- as.character(snakemake@output)
    ggsave(output.path, plots.grid)


}

if (! interactive()){

    main()

}



