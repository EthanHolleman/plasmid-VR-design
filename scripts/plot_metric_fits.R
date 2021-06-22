library(ggplot2)
library(tidyr)

read_agg_distances <- function(dists.filepath){

    df <- as.data.frame(read.table(dists.filepath, sep='\t', header=T))

    # gather over distances to metrics
    df.gather <- gather(df, key='metric_type', value='dist', 
                        distance_bp_prob, distance_local_average_energy, 
                        distance_prop_hairpin, distance_prop_unpaired
                        )
    df.gather

}

