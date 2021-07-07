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

read_rlooper_wig <- function(rlooper.filepath, p_name, id_num, measurement){

    df <- as.data.frame(read.table(as.character(rlooper.filepath), skip=4, header=FALSE))
    colnames(df) <- c('value')
    c(mean(df$value), sd(df$value), nrow(df))

    mean_measurement=paste('mean', measurement, sep='_')
    sd_measurement=paste('sd', measurement, sep='_')
    l <- c(
       mean(df$value), sd(df$value), nrow(df), p_name, id_num
    )
    names(l) <- c(mean_measurement, sd_measurement, 'length', 'name', 'id_num')
    l

}

read_list_rlooper_wig_files <- function(file.list, p_names, id_nums, measurement){

    rows <- list()
    for (i in 1:length(file.list)){
        rows[[i]] <- read_rlooper_wig(
            file.list[[i]], p_names[[i]], id_nums[[i]], measurement
        )

    }

    as.data.frame(do.call(rbind, rows))

}

get_p_names_id_nums_from_rlooper_filelist <- function(file.paths){

    l <- list()
    for (i in 1:length(file.paths)){
        split <- unlist(strsplit(as.character(file.paths[[i]]), '/'))
        # split <- unlist(strsplit(as.character(file.path), '/'))
        # print(split)
        name = split[length(split)]
        p_name =  unlist(strsplit(as.character(name), '[.]'))[1]
        suffix = unlist(strsplit(as.character(name), '[.]'))[2]
        id_num <- as.numeric(unlist(strsplit(as.character(suffix), '_'))[1])
        l[[i]] <- list('name'=p_name, 'id_num'=id_num)
    }

    as.data.frame(do.call(rbind, l))

}


merge_tsv_files_and_param_rlooper_measurements <- function(tsv.df, param.df){

    merge(tsv.df, param.df, on=c('name', 'id_num'))

}


plot_densities <- function(params.merge.df, random.merge.df){

# horizontal line where the random mean for each metric is at
    params.merge.df$GC_skew <- as.numeric(as.character(params.merge.df$GC_skew))
    params.merge.df$GC_content <- as.numeric(as.character(params.merge.df$GC_content))
    params.merge.df$mean_bpprob <- as.numeric(as.character(params.merge.df$mean_bpprob))
    params.merge.df$mean_ale <- as.numeric(as.character(params.merge.df$mean_ale))
    

    # just get the data we need
    params.merge.df.cats <- params.merge.df[, 
                                            c('GC_skew', 'GC_content', 'mean_ale',
                                                'mean_bpprob')
                                            ]
    # label it as parameter informed (not 100% random)
    params.merge.df.cats.type <- params.merge.df.cats
    params.merge.df.cats.type$type <- 'param'

    # add random data
    print('Merging random data!')
    print(colnames(random.merge.df))
    random.merge.df$mean_bpprob <- as.numeric(as.character(random.merge.df$mean_bpprob))
    random.merge.df$mean_ale <- as.numeric(as.character( random.merge.df$mean_ale))
    random.merge.df$type <- 'random'

    gc_skew_levels <- unique(params.merge.df.cats$GC_skew)
    gc_content_levels <- unique(params.merge.df.cats$GC_content)

    extended.random.dfs <- list()
    counter <- 1
    for (each_skew in gc_skew_levels){

        for (each_content in gc_content_levels){

            df <- random.merge.df
            df$GC_skew <- each_skew
            df$GC_content <- each_content
            extended.random.dfs[[counter]] <- df
            counter <- counter + 1
        }

    }
    random.merge.df <- do.call(rbind, extended.random.dfs)

    params.merge.df.cats.all <- as.data.frame(rbind(params.merge.df.cats.type, random.merge.df))

    bpprob <- ggplot(params.merge.df.cats.all, aes(y=mean_bpprob, x=type, fill=type)) +
              geom_boxplot() + stat_compare_means(method='t.test') +
            facet_grid(rows=vars(GC_skew), cols=vars(GC_content)) +
            theme_pubr() + 
            labs(x='Mean per base R-loop probability') +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ale <- ggplot(params.merge.df.cats.all, aes(y=mean_ale, x=type, fill=type)) +
              geom_boxplot() + stat_compare_means(method='t.test') +
            facet_grid(rows=vars(GC_skew), cols=vars(GC_content)) +
            theme_pubr() + 
            labs(x='Mean per base average local energy') +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


    params.merge.df.cats.agg.mean <- aggregate(.~GC_skew+GC_content, params.merge.df.cats, mean)
    colnames(params.merge.df.cats.agg.mean) <- c('GC_skew', 'GC_content', 'average_mean_ale', 'average_mean_bpprob')

    params.merge.df.cats.agg.sd <- aggregate(.~GC_skew+GC_content, params.merge.df.cats, sd)
    colnames(params.merge.df.cats.agg.sd) <-  c('GC_skew', 'GC_content', 'sd_mean_ale', 'sd_mean_bpprob')

    params.merge.df.cats.merge <- merge(
        params.merge.df.cats.agg.mean,  params.merge.df.cats.agg.sd, by=c('GC_skew', 'GC_content')
        )
    bpprob_scatter <- ggplot(params.merge.df.cats.merge, aes(x=GC_skew, y=GC_content, fill=average_mean_bpprob)) +
                        geom_tile() + theme_pubr() +  scale_fill_viridis(discrete=FALSE)
    
    lae_scatter <- ggplot(params.merge.df.cats.merge, aes(x=GC_skew, y=GC_content, fill=average_mean_ale)) +
                        geom_tile() + theme_pubr() +  scale_fill_viridis(discrete=FALSE)
    

    ggarrange(bpprob, ale, bpprob_scatter, lae_scatter, nrow=2, ncol=2)

}




main <- function(){

    random.ale.files <- snakemake@input$random_ale
    random.bpprob.files <- snakemake@input$random_bpprob

    param.ale.files <-  snakemake@input$ale
    params.bprob.files <- snakemake@input$bpprob

    tsv.files <- snakemake@input$tsv_files

    p_names <- snakemake@params$p_names
    id_nums <- snakemake@params$id_nums

    save.image('param.plot.Rdata')

    # read in tsv files for each generated sequence
    param.seqs.df <- read_param_seq_tsv_files(tsv.files)

    random.ale.df <- read_list_rlooper_wig_files(
        random.ale.files, p_names=rep(NA, length(random.ale.files)),
        id_nums=rep(NA, length(random.ale.files)), measurement='ale'
    )
    random.bpprob.df <- read_list_rlooper_wig_files(
        random.bpprob.files, p_names=rep(NA, length(random.bpprob.files)),
        id_nums=rep(NA, length(random.ale.files)), measurement='bpprob'
    )

    random.bpprob.df.merge <- as.data.frame(
        cbind(random.bpprob.df$mean_bpprob,  random.ale.df$mean_ale)
    )
    colnames(random.bpprob.df.merge) <- c('mean_bpprob', 'mean_ale')

    param.ale.file.ids <- get_p_names_id_nums_from_rlooper_filelist(param.ale.files)
    param.ale.df <- read_list_rlooper_wig_files(
          param.ale.files, p_names=param.ale.file.ids$name, 
          id_nums=param.ale.file.ids$id_num,
          measurement='ale'
    )
   

    param.bpprob.file.ids <- get_p_names_id_nums_from_rlooper_filelist(params.bprob.files)
    param.bpprob.df <- read_list_rlooper_wig_files(
           params.bprob.files, p_names= param.bpprob.file.ids$name, 
          id_nums= param.bpprob.file.ids$id_num,
          measurement='bpprob'
    )
    

    # add ale and bprob to param tsv df
    param.seqs.df.ale <- merge_tsv_files_and_param_rlooper_measurements(
         param.seqs.df,  param.ale.df
    )
    stopifnot(nrow(param.seqs.df.ale)==nrow(param.ale.df))


    param.seqs.df.bpprob <- merge_tsv_files_and_param_rlooper_measurements(
         param.seqs.df, param.bpprob.df
    )
    stopifnot(nrow(param.bpprob.df) == nrow(param.seqs.df.bpprob))


    param.seqs.df.ale.bpprob <- cbind(param.seqs.df.ale, param.seqs.df.bpprob[, c('mean_bpprob', 'sd_bpprob')])
    plot_dense <- plot_densities(param.seqs.df.ale.bpprob, random.bpprob.df.merge)
    ggsave(as.character(snakemake@output$plot), plot_dense, height=24, width=24)

}

if (! interactive()){
    main()
}