

#Date: 28may2025
# This code is used to check what edge degree threshold we need to make grande break
# at ~291 cells

library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
# library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggridges)
library(ghibli)
library(ggpubr)
library(tidyverse)
library(ggnewscale)
library(cowplot)
library(png)
library(grid)

library(slider)


theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

setwd("/Users/pipe/Desktop/edge_degree_selection_28may2025")

# Loading the network simulation data ####

# frag_df=data.frame()
# 
# for(j in seq(5, 20)){
#   temp_df=read.csv(paste("grande_", j, "e_100g_random/fragmentation_inf.csv", sep=""), header=TRUE)
#   temp_df$strain='Grande'
#   temp_df$edge_degree_threshold=j
#   frag_df=rbind(frag_df, temp_df)
# }
# 
# write.csv(frag_df, "~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig3_edge_degree_selection/edge_selection_grande_frag_28may2026.csv",
#           row.names = FALSE)

frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig3_edge_degree_selection/edge_selection_grande_frag_28may2026.csv",
                 header=TRUE)
summary(frag_df)

ggplot(frag_df, aes(x=edge_degree_threshold, y=cluster_size, fill=factor(edge_degree_threshold)))+
  geom_violin()+
  guides(fill='none')+
  scale_y_log10()+
  geom_hline(yintercept=291)+
  NULL

summ_frag_df=frag_df %>%
  group_by(strain, edge_degree_threshold) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size))

ggplot(summ_frag_df, aes(x=edge_degree_threshold, y=mean_size, col=strain))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(x=edge_degree_threshold, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=strain), 
              alpha=0.3, linetype='blank')+
  guides(fill='none', col='none')+
  # scale_y_log10()+
  geom_hline(yintercept=291)+
  NULL

# Fracture size at edge degree of 14: 224.34376
# Fracture size at edge degree of 15: 339.82627
# Fracture size at edge degree of 16: 451.26935
# I need to do the test with edge degree of 15, because the results of the paper are already with 14
# Note: even if we think of using the 95% percentile, that is 394 cells, so we can still argue that
# 15 is a better value.


# Cluster size distribution ####
#format of the csv files
#(Date)_(strain [gob8/gob21])-(replicate [1-5])_(image number)_measurements_rois_border.csv

load_one_csv_cluster_size <- function(file) {
  df <- read.csv(file, row.names = 1) #reading the csv file
  measure_vars <- colnames(df) #getting the names of the columns
  filename <- basename(file) #saving only the name of the file, so not saving the path to the file
  
  #removing the last part of the name of the file "_measurements_rois.csv" to only keep the important information
  
  metadata <- strsplit(filename, split = "_", fixed = TRUE)[[1]] #splitting the metadata into the parts of the information
  temp_date=metadata[1]
  temp_strain_replicate=metadata[2]
  if(strsplit(temp_strain_replicate, split='-', fixed=TRUE)[[1]][1]=='gob21'){
    temp_strain='petite'
  } else {
    temp_strain='grande'
  }
  temp_image_num=metadata[3]
  temp_replicate=strsplit(temp_strain_replicate, split='-', fixed=TRUE)[[1]][2]
  
  df$date=temp_date
  df$strain=temp_strain
  df$image_num=temp_image_num
  df$replicate=temp_replicate
  
  #creating the final data frame
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}

# Loading cluster size data ####
#Directory where the images are located


# in_dir="~/work_dir/observed_synchrony/data/Microscopy/all_ancestor_measurements/cluster_measurements"
# 
# data_clust <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
#                     .fun = load_one_csv_cluster_size)
# data_clust$volume=(4/3)*pi*((data_clust$Major/2)^3)

data_clust=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/ancestor_cluster_measurements_23may2024.csv",
                    header=TRUE)
data_clust$strain=factor(data_clust$strain, levels=c('grande', 'petite'))
summary(data_clust)



# Number of cells per cluster ####


# Format of the csv files:
# (Date)_(strain [gob8])-(replicate [1-5])_(image number)_measurements_by_group.csv

load_one_csv_cell_distribution <- function(file) {
  df <- read.csv(file, row.names = 1) # reading the csv file
  filename <- basename(file) # saving only the name of the file
  
  # Parse metadata from filename
  metadata <- strsplit(filename, split = "_", fixed = TRUE)[[1]]
  temp_date <- metadata[1]
  temp_strain_replicate <- metadata[2]
  temp_image_num <- metadata[3]
  
  # Extract strain information
  strain_part <- strsplit(temp_strain_replicate, split = '-', fixed = TRUE)[[1]][1]
  if(strain_part == 'gob21') {
    temp_strain <- 'petite'
  } else {
    temp_strain <- 'grande'
  }
  
  # Extract replicate
  temp_replicate <- strsplit(temp_strain_replicate, split = '-', fixed = TRUE)[[1]][2]
  
  # Filter for rows where group is NA (these are the group segmentations)
  group_rows <- df[is.na(df$group), ]
  
  # Create group_id from row names
  group_rows$group_id <- rownames(group_rows)
  
  # Calculate num_cells for each group
  # For each group_id, count how many individual cells belong to that group
  group_rows$num_cells <- sapply(group_rows$group_id, function(group_id) {
    # Count rows where group column equals this group_id
    sum(df$group == as.numeric(group_id), na.rm = TRUE)
  })
  
  # Add metadata columns
  group_rows$date <- temp_date
  group_rows$strain <- temp_strain
  group_rows$replicate <- temp_replicate
  group_rows$image_num <- temp_image_num
  
  # Remove the group column since it's all NA for these rows
  group_rows$group <- NULL
  
  # Filter out groups with no cells (segmentation errors)
  group_rows <- group_rows[group_rows$num_cells > 0, ]
  
  return(group_rows)
}


#Directory where the images are located
# in_dir="~/work_dir/observed_synchrony/data/Microscopy/2025may19_ancestors_measurements/cell_distribution/cell_counts_results"
# 
# cell_dist <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
#                    .fun = load_one_csv_cell_distribution)

cell_dist=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/ancestor_cluster_cell_counts_23may2024.csv",
                   header=TRUE)
cell_dist$strain=factor(cell_dist$strain, levels=c('grande', 'petite'))
summary(cell_dist)


ggplot(cell_dist, aes(x=as.factor(replicate), y=num_cells))+
  geom_violin(fill="#F8766D")+
  theme_classic()+
  facet_wrap(~strain)+
  ylab("# of cells per cluster")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_dist, aes(x=as.factor(replicate), y=num_cells))+
  geom_violin(fill="#F8766D")+
  theme_classic()+
  facet_wrap(~strain)+
  ylab("# of cells per cluster")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  scale_y_continuous(trans='log10')+
  NULL



ggplot(cell_dist, 
       aes(x=num_cells, group=interaction(replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("# of cells per cluster")+
  scale_x_continuous(trans='log10')+
  facet_wrap(~strain, ncol=1)+
  geom_vline(xintercept = 25)+
  geom_vline(xintercept = 2.739821, linetype='dashed')+
  geom_vline(xintercept = 274.3413, linetype='dashed')+
  NULL

mean(cell_dist$num_cells)
# 79.53882
median(cell_dist$num_cells)
# 3

mean(cell_dist[cell_dist$num_cells>25,]$num_cells)
# 274.3413
median(cell_dist[cell_dist$num_cells>25,]$num_cells)
# 236.5


# Number of clusters imaged in total
table(cell_dist$replicate)
# 1   2   3   4   5 
# 483 550 702 548 525 

table(cell_dist[cell_dist$num_cells>25,]$replicate)
# 1   2   3   4   5 
# 144 145 165 180 160


ggplot(cell_dist[cell_dist$num_cells>25,], aes(x=as.factor(replicate), y=num_cells))+
  geom_violin(fill="#F8766D")+
  theme_classic()+
  facet_wrap(~strain)+
  ylab("# of cells per cluster")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  scale_y_continuous(trans='log10')+
  NULL




# Is there a correlation between cluster area an number of cells in cluster?
ggplot(cell_dist, aes(x=Area, y=num_cells, col=replicate))+
  geom_point(alpha=0.5)+
  theme_classic()+
  ylab("# of cells per cluster")+
  guides(col='none')+
  NULL


ggplot(cell_dist, 
       aes(x=num_cells, group=interaction(replicate)))+
  geom_histogram(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("# of cells per cluster")+
  ylim(c(0,50))+
  #scale_x_continuous(trans='log10')+
  facet_wrap(~replicate, ncol=1)+
  # geom_vline(xintercept = 25)+
  # geom_vline(xintercept = 2.739821, linetype='dashed')+
  # geom_vline(xintercept = 274.3413, linetype='dashed')+
  NULL

ggplot(data_clust, 
       aes(x=(Major+Minor)/4, group=interaction(replicate)))+
  geom_histogram(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  # xlab("# of cells per cluster")+
  # ylim(c(0,25))+
  #scale_x_continuous(trans='log10')+
  facet_wrap(~replicate, ncol=1)+
  NULL

quantile((data_clust$Major+data_clust$Minor)/4, c(0.97))


ggplot(data_clust, 
       aes(x=(sqrt(Area/pi)), group=interaction(replicate)))+
  geom_histogram(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  # xlab("# of cells per cluster")+
  # ylim(c(0,25))+
  #scale_x_continuous(trans='log10')+
  facet_wrap(~replicate, ncol=1)+
  NULL


# Number of cells:
quantile(cell_dist$num_cells, c(0.90, 0.95, 0.96, 0.97))

# Cluster radius:
quantile(sqrt(data_clust[data_clust$strain=='grande',]$Area/pi), c(0.9, 0.95, 0.96, 0.97))


# paper figure ####

quant_radius=quantile(sqrt(data_clust[data_clust$strain=='grande',]$Area/pi), c(0.9))

cluster_radius=ggplot(data_clust[data_clust$strain=='grande',], 
       aes(x=(sqrt(Area/pi)), fill=strain))+
  geom_histogram(col='black', fill='#AE93BEFF')+
  theme_classic()+
  guides(alpha = "none")+
  ylab("# of Clusters")+
  xlab('Cluster Radius')+
  xlim(c(0,50))+
  guides(fill='none')+
  geom_vline(xintercept=quant_radius, linetype='dashed', size=0.5)+
  theme_classic(base_size = 10)+
  NULL
cluster_radius


quant_cell_dist=quantile(cell_dist$num_cells, c(0.90))

cells_per_cluster=ggplot(cell_dist, aes(x=num_cells, fill=strain))+
  geom_histogram(col='black', fill='#AE93BEFF')+
  theme_classic()+
  guides(fill = "none")+
  xlab("# of Cells per Cluster")+
  ylab('Count')+
  # ylim(c(0, 250))+
  scale_x_log10()+
  geom_vline(xintercept=quant_cell_dist, linetype='dashed', size=0.5)+
  theme_classic(base_size = 10)+
  NULL
cells_per_cluster


edge_degree=ggplot(summ_frag_df, aes(x=edge_degree_threshold, y=mean_size, col=strain))+
  geom_line(col='#AE93BEFF')+
  geom_point(col='#AE93BEFF')+
  geom_ribbon(aes(x=edge_degree_threshold, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=strain), 
              alpha=0.3, linetype='blank', fill='#AE93BEFF')+
  guides(fill='none', col='none')+
  # scale_y_log10()+
  geom_hline(yintercept=quant_cell_dist, linetype='dashed', size=0.5)+
  labs(x="Edge Degree Threshold", y='Fracture Size')+
  theme_classic(base_size = 10)+
  NULL
edge_degree


supp_fig=plot_grid(cluster_radius, cells_per_cluster, edge_degree,
                                   labels=c('A', 'B', 'C'), ncol=1, label_size=11)
supp_fig

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_3_edge_degree_selection.png',
       plot=supp_fig, dpi='retina', width=3.5, height=6.3, bg='white')
