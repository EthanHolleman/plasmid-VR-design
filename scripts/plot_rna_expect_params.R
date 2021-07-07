library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)


read_param_seq_tsv_files <- function(tsv.files){

    tsv.rows <- list()
    for (i in 1:length(tsv.files)){
        cur_file = as.character(tsv.files[[i]])
        tsv.rows[[i]] <- as.data.frame(
            read.table(cur_file, sep='\t', header=T)
        )
    
    }
    do.call(rbind, tsv.rows)


}

