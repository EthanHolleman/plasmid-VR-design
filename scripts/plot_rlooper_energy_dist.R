# read all files at once?
library(ggplot2)
library(ggpubr)

read_average_local_energy_file <- function(rlooper.filepath){

    ale.df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(ale.df) <- c('ale')
    mean(ale.df$ale)

}


read_all_local_average_energy_means <- function(all.ale.filepaths){

    means <- list()
    for (i in 1:length(all.ale.filepaths)){

        means[[i]] <- read_average_local_energy_file(all.ale.filepaths[[i]])

    }
    ale.means.df <- as.data.frame(do.call(rbind, means))
    colnames(ale.means.df) <- c('mean')
    ale.means.df

}


plot_ale_means <- function(ale.means.df){

    num_seqs <- nrow(ale.means.df)
    ggplot(ale.means.df, aes(x=mean)) + geom_density(alpha=0.7, fill='firebrick') +
        theme_pubr() + 
        labs(title=paste('Mean LAE of', num_seqs, 'of length 200'), x='LAE')
    
}

main <- function(){

    all.input.files <- snakemake@input
    output.path <- as.character(snakemake@output)
    ale.means.df <- read_all_local_average_energy_means(all.input.files)
    print(ale.means.df)
    plot <- plot_ale_means(ale.means.df)
    ggsave(output.path, plot, dpi=500)

}

if (! interactive()){
    main()
}