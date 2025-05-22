
# Date: 15 May 2025
# Code to make figure about the cluster size distribution in the grande and petite ancestors


library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggridges)
library(ghibli)
# library(ggpubr)
# library(tidyverse)
library(ggnewscale)
library(cowplot)
library(png)
library(grid)



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

in_dir="~/work_dir/observed_synchrony/data/Microscopy/2025may19_ancestors_measurements/cluster_measurements"

data_clust <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                    .fun = load_one_csv_cluster_size)
data_clust$volume=(4/3)*pi*((data_clust$Major/2)^3)
data_clust$strain=factor(data_clust$strain, levels=c('petite', 'grande'))
summary(data_clust)


# Grande and petite volume measurements ####

summary(data_clust)


median_volume=data_clust %>%
  group_by(strain, replicate) %>%
  summarise(mean_vol=mean(volume),
            median_vol=median(volume),
            percentile_25=quantile(volume, c(.25)),
            percentile_75=quantile(volume, c(.75)))



# Cluster Volume
grande_mean_volume=mean(data_clust[data_clust$strain=='grande',]$volume)
grande_mean_volume
# 106836.2 um^3
petite_mean_volume=mean(data_clust[data_clust$strain=='petite',]$volume)
petite_mean_volume
# 54170.07 um^3

# Cluster Radius
mean(data_clust[data_clust$strain=='grande',]$Major)/2
# 27.59116
mean(data_clust[data_clust$strain=='petite',]$Major)/2
# 21.89868


ggplot(data_clust, aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("log10 Cluster Volume (Microns^3)")+
  geom_vline(xintercept=grande_mean_volume, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col='#00BFC4', linetype='dashed')+
  scale_x_continuous(trans='log10')+
  NULL



ggplot(data_clust, aes(x=as.factor(replicate), y=volume, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_wrap(~strain)+
  scale_y_continuous(trans='log10')+
  ylab("log10 Cluster Volume (Microns^3)")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL


ggplot(data_clust, aes(x=volume, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("log10 Cluster Volume (Microns^3)")+
  geom_vline(xintercept=petite_mean_volume, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=grande_mean_volume, col='#00BFC4', linetype='dashed')+
  scale_x_continuous(trans='log10')+
  NULL


# Loading Cell Properties Data ####


#format of the csv files
#(Date)_(strain [gob8/gob21])-(replicate [1-5])_(image number)_measurements_rois_border.csv

load_one_csv_cell_measurements <- function(file) {
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

#Directory where the images are located
in_dir="~/work_dir/observed_synchrony/data/Microscopy/2025may19_ancestors_measurements/cell_properties"

cell_data <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                   .fun = load_one_csv_cell_measurements)

cell_data$strain=factor(cell_data$strain, levels=c('petite', 'grande'))
summary(cell_data)

# Number of cells characterized in the second day
table(cell_data$replicate, cell_data$strain)
#   petite grande
# 1   1043   1285
# 2   1540   1462
# 3    892   1533
# 4   1058   1517
# 5   1334   1770

#### Aspect Ratio ####

median_ar=cell_data %>%
  group_by(strain, replicate) %>%
  summarise(mean_ar=mean(AR),
            median_ar=median(AR),
            percentile_25=quantile(AR, c(.25)),
            percentile_75=quantile(AR, c(.75)))



grande_mean_ar=mean(cell_data[cell_data$strain=='grande',]$AR)
grande_mean_ar
# mean: 1.19564
# median: 1.18784
petite_mean_ar=mean(cell_data[cell_data$strain=='petite',]$AR)
petite_mean_ar
# mean: 1.222915
# median: 1.21388


ggplot(cell_data, aes(x=as.factor(replicate), y=AR, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_wrap(~strain)+
  ylab("Cell Aspect Ratio")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data, 
       aes(x=AR, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("Cell Aspect Ratio")+
  geom_vline(xintercept=petite_mean_ar, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=grande_mean_ar, col='#00BFC4', linetype='dashed')+
  facet_wrap(~strain, ncol=1)+
  NULL


#### Cell Diameter ####

summ_diam=cell_data %>%
  group_by(strain, replicate) %>%
  summarise(mean_diam=mean(Major),
            median_diam=median(Major),
            percentile_25=quantile(Major, c(.25)),
            percentile_75=quantile(Major, c(.75)))


grande_mean_diam=mean(cell_data[cell_data$strain=='grande',]$Major)
grande_mean_diam
# mean: 5.170423 um
# median: 5.14785 um
petite_mean_diam=mean(cell_data[cell_data$strain=='petite',]$Major)
petite_mean_diam
# mean: 4.844126 um
# median: 4.82025 um


ggplot(cell_data, aes(x=as.factor(replicate), y=Major, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_wrap(~strain)+
  ylab("Cell Diameter")+
  xlab("Replicate")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data, 
       aes(x=Major, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("Cell Diameter")+
  geom_vline(xintercept=petite_mean_diam, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=grande_mean_diam, col='#00BFC4', linetype='dashed')+
  # scale_x_continuous(trans='log10')+
  facet_wrap(~strain, ncol=1)+
  NULL



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
in_dir="~/work_dir/observed_synchrony/data/Microscopy/2025may19_ancestors_measurements/cell_distribution/cell_counts_results"

cell_dist <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                   .fun = load_one_csv_cell_distribution)

cell_dist$strain=factor(cell_dist$strain, levels=c('petite', 'grande'))
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


# Paper Figure ####


#### version 2 ####

p_cluster_vol_log=ggplot(data_clust, aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.75))+
  theme_classic(base_size = 14)+
  guides(alpha = "none")+
  ylab("Density")+
  xlab(expression(log[10]~"Cluster Volume"~(mu*m^3)))+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF"))+
  geom_vline(xintercept=grande_mean_volume, col="#F0D77BFF", linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col="#AE93BEFF", linetype='dashed')+
  scale_x_continuous(trans='log10')+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.line.y = element_blank()) +  # Remove y-axis line
  labs(fill = "Strain")+
  NULL
p_cluster_vol_log

p_cluster_vol=ggplot(data_clust, aes(x=volume, fill=strain))+
  # geom_density(aes(fill=strain, alpha=0.5))+
  geom_histogram(position="identity", alpha=0.75, bins=30) +
  theme_classic(base_size = 14)+
  guides(alpha = "none")+
  ylab("# of clusters")+
  xlab(expression("Cluster Volume"~(µ*m^3)))+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF")) +
  geom_vline(xintercept=grande_mean_volume, col="#AE93BEFF", linetype='dashed') +
  geom_vline(xintercept=petite_mean_volume, col="#F0D77BFF", linetype='dashed') +
  xlim(c(0, 250000))+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(fill = "Strain")+
  NULL
p_cluster_vol

p_ar=ggplot(cell_data, 
            aes(x=AR, fill=strain))+
  geom_histogram(alpha=0.75, position="identity", bins=30) +
  theme_classic(base_size = 14)+
  guides(alpha = "none", fill="none")+
  ylab("# of cells")+
  xlab("Aspect Ratio")+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF")) +
  geom_vline(xintercept=grande_mean_ar, col="#AE93BEFF", linetype='dashed')+
  geom_vline(xintercept=petite_mean_ar, col="#F0D77BFF", linetype='dashed')+
  NULL
p_ar


p_diam=ggplot(cell_data, 
              aes(x=Major, fill=strain))+
  geom_histogram(alpha=0.75, position="identity", bins=30) +
  theme_classic(base_size = 14)+
  guides(alpha = "none", fill="none")+
  ylab("# of cells")+
  # xlab(expression("Cell Diameter"~(mu*m)))+
  xlab("Cell Diameter (µm)")+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF"))+
  geom_vline(xintercept=grande_mean_diam, col="#AE93BEFF", linetype='dashed')+
  geom_vline(xintercept=petite_mean_diam, col="#F0D77BFF", linetype='dashed')+
  NULL
p_diam


img_petite_clust <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/gob21-2_day2_001_620x620_100um-scale.png")
img_plot_petite_clust <- rasterGrob(img_petite_clust, interpolate = TRUE)

img_grande_clust <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/gob8-1_day2_001_620x620_100um-scale.png")
img_plot_grande_clust <- rasterGrob(img_grande_clust, interpolate = TRUE)


img_petite_ar <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/gob21-1_day2_003_no_masks_10um-scale.png")
img_plot_petite_ar <- rasterGrob(img_petite_ar, interpolate = TRUE)

img_grande_ar <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/gob8-1_day2_006_no_masks_10um-scale.png")
img_plot_grande_ar <- rasterGrob(img_grande_ar, interpolate = TRUE)



col_clust_vol=plot_grid(img_plot_petite_clust, img_plot_grande_clust, 
                        ncol = 1, 
                        rel_heights = c(1, 1))

col_ar=plot_grid(img_plot_petite_ar, img_plot_grande_ar, 
                 ncol = 1, 
                 rel_heights = c(1, 1))

plot_c_ar_diam=plot_grid(p_ar, p_diam, nrow=1, labels=c('C', 'D'), label_size = 16)
plot_c_ar_diam

fig_1_v2=plot_grid(p_cluster_vol, col_clust_vol, plot_c_ar_diam, col_ar, 
                   labels=c('A', 'B','','E'), nrow=2, label_size=16, rel_widths=c(1, 0.35), rel_heights=c(1,1))
fig_1_v2

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_1_clust_size_ancestors_22may2025_v2.svg',
       plot=fig_1_v2, dpi='retina', width=8, height=6, bg='white')
# note: image is saved as an svg to later add the strain labels for the images
