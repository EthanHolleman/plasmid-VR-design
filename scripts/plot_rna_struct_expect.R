library(ggplot2)
library(ggpubr)
library(ggridges)

# read output produced by parsing structure file
read_rna_struct_tsv <- function(file.path){

    df <- as.data.frame(read.table(file.path, header=F))
    colnames(df) <- c('filepath', 'prop_unpaired', 'length', 
                      'num_hairpins', 'prop_hairpin'
                      )
    df

}

read_all_rna_struct_tsvs <- function(file.path.list){

    dfs <- list()
    for (i in 1:length(file.path.list)){

        file.path <- as.character(file.path.list[[i]])
        dfs[[i]] <- read_rna_struct_tsv(file.path)

    }

    do.call(rbind, dfs)

}


plot_rna_struct <- function(big.df){

    prop_unpaired_vs_prop_hairpin <- ggplot(big.df, aes(x=prop_unpaired, y=prop_hairpin), color=num_hairpins) +
                                     geom_point(alpha=0.5) + theme_pubr() + 
                                     labs(
                                         x='Proportion sequence unpaired',
                                         y='Proportion sequence in hairpin'
                                    )
    dist_hairpin <- ggplot(big.df, aes(x=prop_hairpin, y=as.factor(length))) +
                    geom_density_ridges(alpha=0.7) + theme_pubr() +
                    labs(x='Proportion sequence in hairpin', y='Sequence length')
    dist_unpaired <- ggplot(big.df, aes(x=prop_unpaired, y=as.factor(length))) +
                    geom_density_ridges(alpha=0.7) + theme_pubr() +
                    labs(x='Proportion sequence unpaired', y='Sequence length')
    
    ggarrange(prop_unpaired_vs_prop_hairpin,
            ggarrange(dist_hairpin, dist_unpaired, nrow=1, ncol=2),
            nrow=2, ncol=1
    )
}


agg_variable_mean_sd <- function(big.df, var_name, group_name){

    means <- aggregate(big.df[, var_name], list(big.df[, group_name]), mean)
    sds <- aggregate(big.df[, var_name], list(big.df[, group_name]), sd)
    merge.df <- merge(means, sds, by='Group.1')
    merge.df$attribute <- var_name
    colnames(merge.df) <- c(group_name, 'mean', 'sd', 'attribute')
    merge.df

}

write_expectation_table <- function(big.df, output.path){

    means.hairpin <- agg_variable_mean_sd(big.df, 'prop_hairpin', 'length')
    means.prop_unpaired <- agg_variable_mean_sd(big.df, 'prop_unpaired', 'length')
    all.means <- rbind(means.hairpin, means.prop_unpaired)
    write.table(all.means, file=output.path, sep='\t', quote=FALSE, row.names=FALSE)

}

main <- function(){

    input.files.list <- snakemake@input
    output.path.plot <- as.character(snakemake@output$plot)
    output.path.expect <- as.character(snakemake@output$expect)
    big.df <- read_all_rna_struct_tsvs(input.files.list)
    plot <- plot_rna_struct(big.df)
    expect_table <- write_expectation_table(big.df, output.path.expect)
    ggsave(output.path.plot, plot, dpi=500, width=10, height=10)


}

if (!interactive()){
    main()
}