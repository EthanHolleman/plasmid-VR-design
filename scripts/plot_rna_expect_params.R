library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)


read_all_agg_metric_files <- function(file.path.list){

    dfs <- list()
    for (i in 1:length(file.path.list)){

        file.path <- as.character(file.path.list[[i]])
        dfs[[i]] <- as.data.frame(read.table(file.path, sep='\t', header=TRUE))

    }

    do.call(rbind, dfs)

}



grid_plots <- function(big.df){

    # at this point need df with only GC skew content and RNA metrics
    big.df.essential <- big.df[, c('GC_skew', 'GC_content', 'prop_unpaired', 'prop_hairpin')]

    big.df.agg <- aggregate(.~GC_skew+GC_content, big.df.essential, mean)
    colnames(big.df.agg) <-  c('GC_skew', 'GC_content', 'mean_prop_unpaired', 'mean_prop_hairpin')

    prop_hairpin <- ggplot(big.df.agg, aes(GC_skew, GC_content, fill=mean_prop_hairpin)) +
            geom_tile() + theme_pubr() + scale_fill_viridis(discrete=FALSE) +
            labs(x='GC skew', y='GC content', fill='Mean proportion bases in hairpin') +
             theme(legend.key.size = unit(2, 'cm')) +
            theme(text = element_text(size=24))

    prop_unpaired <- ggplot(big.df.agg, aes(GC_skew, GC_content, fill=mean_prop_unpaired)) +
            geom_tile() + theme_pubr() + scale_fill_viridis(discrete=FALSE) +
            labs(x='GC skew', y='GC content', fill='Mean proportion unpaired bases') +
             theme(legend.key.size = unit(2, 'cm')) +
            theme(text = element_text(size=24))
    
    ggarrange(prop_hairpin, prop_unpaired, nrow=1, ncol=2, labels=c('A', 'B'))


}

main <- function(){


    agg_metrics.files <- snakemake@input$metrics
    tsv.files <- snakemake@input$tsv_files
    save.image('debug.plot_rna_expect_params.RData')

    metrics.df <- read_all_agg_metric_files(agg_metrics.files)
    metrics.df$name <- metrics.df$p_name
    tsv.df <- read_all_agg_metric_files(tsv.files)

    big.df <- merge(metrics.df, tsv.df, by=c('name', 'id_num'))

    plots.grid <- grid_plots(big.df)

    output.path <- as.character(snakemake@output)
    ggsave(output.path, plots.grid, width=18, height=14)


}

if (! interactive()){

    main()

}



