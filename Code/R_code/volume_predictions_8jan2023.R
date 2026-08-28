
#Date: 8 January 2024
# Code used to analyze the results of snowflake_volume_overlap_approximation.m

library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(igraph)

setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/table1_volume_predictions/')


# PA ancestors measurements ####
# This code is from ancestor_clust_size_diff_22may2025.R

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

in_dir="~/work_dir/observed_synchrony/data/Microscopy/all_ancestor_measurements/cluster_measurements"

data_clust <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                    .fun = load_one_csv_cluster_size)
data_clust$volume=(4/3)*pi*((data_clust$Major/2)^3)
data_clust$strain=factor(data_clust$strain, levels=c('grande', 'petite'))
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
# 105010.9 um^3

petite_mean_volume=mean(data_clust[data_clust$strain=='petite',]$volume)
petite_mean_volume
# 58179.76 um^3

grande_mean_volume-petite_mean_volume
# 46831.11

#### All data ####

# Cluster volume
# Grande mean: 105010.9 um^3
# Petite mean: 58179.76 um^3

# Aspect Ratio
# Grande mean: 1.197828
# Petite mean: 1.221764

# Cell Diameter
# Grande mean: 5.085207
# Petite mean: 4.844297



# Aspect ratio ####
#### Threshold prediction ####

aspr_volume_df=read.csv("exponential_clust_volume_1.197aspr.csv", header=TRUE)

summary(aspr_volume_df)

mean_volume=aspr_volume_df %>%
  group_by(overlap_thresh) %>%
  summarize(mean_volume=mean(volume_clust),
            mean_num_cells=mean(num_cells))

ggplot(aspr_volume_df, aes(x=overlap_thresh, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_volume, aes(x=overlap_thresh, y=mean_volume))+
  theme_bw()+
  geom_hline(yintercept=grande_mean_volume, col='red', linetype='dashed')+
  # scale_x_continuous(trans='log10')+
  NULL


ggplot(aspr_volume_df, aes(x=overlap_thresh, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_volume, aes(x=overlap_thresh, y=mean_num_cells))+
  theme_bw()+
  # scale_x_continuous(trans='log10')+
  NULL

ggplot(aspr_volume_df, aes(x=overlap_thresh, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='green')+
  theme_bw()+
  # scale_x_continuous(trans='log10')+
  geom_hline(yintercept=grande_mean_volume, col='red', linetype='dashed')+
  ylab('Cluster Volume')+
  xlab('Overlap Threshold')+
  NULL

ggplot(aspr_volume_df, aes(x=overlap_thresh, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  theme_bw()+
  # scale_x_continuous(trans='log10')+
  NULL

fit_volume=glm(aspr_volume_df$volume_clust ~ aspr_volume_df$overlap_thresh, family = gaussian)
fit_volume

b0 <- fit_volume$coefficients[1]
b1 <- fit_volume$coefficients[2]

(grande_mean_volume-b0)/b1
#Overlap needed to simulate clusters of the correct volume
# 63.45398 #making it check the overlap after each cell added



#### Size difference AR ####

aspr_df=read.csv("constant_threshold_changing_aspr_2025june12.csv", header=TRUE)

summary(aspr_df)

mean_aspr=aspr_df %>%
  group_by(aspect_ratio) %>%
  summarise(mean_vol=mean(volume_clust),
            mean_num_cells=mean(num_cells))

ggplot(aspr_df, aes(x=aspect_ratio, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_aspr, aes(x=aspect_ratio, y=mean_vol))+
  theme_bw()+
  # scale_y_continuous(trans='log10')+
  NULL

ggplot(aspr_df, aes(x=aspect_ratio, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  # geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='red')+
  theme_bw()+
  ylab('Cluster Volume')+xlab('Cellular Aspect Ratio')+
  NULL
#Better approximated by a square function

ggplot(aspr_df, aes(x=aspect_ratio, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  # geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='red', se=TRUE) +
  theme_bw()+
  ylab('Number of cells')+xlab('Cellular Aspect Ratio')+
  NULL


# Expected difference in volume
volume_model=glm(aspr_df$volume_clust ~ poly(aspr_df$aspect_ratio, 2, raw=TRUE), family=gaussian)
vol_coeffs=volume_model$coefficients

# Aspect Ratio
# Grande mean: 1.197828
# Petite mean: 1.221764
petite_asp_r=1.221764
grande_asp_r=1.197828

predic_petite_aspr=vol_coeffs[1]+(vol_coeffs[2]*petite_asp_r)+(vol_coeffs[3]*(petite_asp_r^2))
# 69875.07
predic_grande_aspr=vol_coeffs[1]+(vol_coeffs[2]*grande_asp_r)+(vol_coeffs[3]*(grande_asp_r^2))
# 63982.21
predic_grande_aspr-predic_petite_aspr # -5892.858
# Petite strain has 5892.858 um^3 more

# Expected difference in number of cells
cells_model=glm(aspr_df$num_cells ~ poly(aspr_df$aspect_ratio, 2, raw=TRUE), family=gaussian)
cell_coeffs=cells_model$coefficients


cell_coeffs[1]+(cell_coeffs[2]*petite_asp_r)+(cell_coeffs[3]*(petite_asp_r^2))
# 331.0742
cell_coeffs[1]+(cell_coeffs[2]*grande_asp_r)+(cell_coeffs[3]*(grande_asp_r^2))
# 319.7258
319.7258-331.0742 # -11.3484
#Petite strain has 11 more cells



# Cell Diameter ####
#### Threshold prediction ####

diam_volume_df=read.csv("exponential_clust_volume_5.085cell_diam.csv", header=TRUE)

summary(diam_volume_df)

mean_volume=diam_volume_df %>%
  group_by(overlap_thresh) %>%
  summarize(mean_volume=mean(volume_clust),
            mean_num_cells=mean(num_cells))

ggplot(diam_volume_df, aes(x=overlap_thresh, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_volume, aes(x=overlap_thresh, y=mean_volume))+
  theme_bw()+
  # scale_x_continuous(trans='log10')+
  geom_hline(yintercept = grande_mean_volume)+
  NULL


ggplot(diam_volume_df, aes(x=overlap_thresh, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_volume, aes(x=overlap_thresh, y=mean_num_cells))+
  theme_bw()+
  scale_x_continuous(trans='log10')+
  NULL

ggplot(diam_volume_df, aes(x=overlap_thresh, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  # scale_x_continuous(trans='log10')+
  geom_hline(yintercept = grande_mean_volume)+
  NULL

ggplot(diam_volume_df, aes(x=overlap_thresh, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  theme_bw()+
  scale_x_continuous(trans='log10')+
  NULL


fit_volume_diam=glm(diam_volume_df$volume_clust ~ diam_volume_df$overlap_thresh, family = gaussian)
fit_volume_diam

b0_diam <- fit_volume_diam$coefficients[1]
b1_diam <- fit_volume_diam$coefficients[2]

(grande_mean_volume-b0_diam)/b1_diam
# Threshold needed: 67.06724



#### Size Difference DIAM ####

diam_df=read.csv("constant_threshold_changing_cell_diam_2025june12.csv", header=TRUE)

summary(diam_df)


# Cell Diameter
# Grande mean: 5.085207
# Petite mean: 4.844297

petite_cell_diam=4.844297
grande_cell_diam=5.085207

mean_diam=diam_df %>%
  group_by(cell_diam) %>%
  summarise(mean_vol=mean(volume_clust),
            mean_num_cells=mean(num_cells))

ggplot(diam_df, aes(x=cell_diam, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_diam, aes(x=cell_diam, y=mean_vol))+
  theme_bw()+
  # scale_y_continuous(trans='log10')+
  NULL

ggplot(diam_df, aes(x=cell_diam, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="red", fill="red", se=TRUE)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  geom_vline(xintercept=grande_cell_diam, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_cell_diam, col='#00BFC4', linetype='dashed')+
  xlab('Cell Diameter (microns)')+ylab('Cluster Volume')+
  NULL


ggplot(diam_df, aes(x=cell_diam, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  NULL


# Expected difference in volume
volume_model=glm(diam_df$volume_clust ~ diam_df$cell_diam, family=gaussian)
vol_coeffs=volume_model$coefficients


petite_diam_pred=vol_coeffs[1]+(vol_coeffs[2]*petite_cell_diam)
petite_diam_pred
# 63225.55 
grande_diam_pred=vol_coeffs[1]+(vol_coeffs[2]*grande_cell_diam)
grande_diam_pred
# 74104.42
grande_diam_pred-petite_diam_pred
# Grande strain  has 10878.87 um^3 more

# Note: for some reason the obtained volume is much lower from the 77 thousand which was


# Expected difference in number of cells
cells_model=glm(diam_df$num_cells ~ diam_df$cell_diam, family=gaussian)
cell_coeffs=cells_model$coefficients

petite_diam_cell_pred=cell_coeffs[1]+(cell_coeffs[2]*petite_cell_diam)
petite_diam_cell_pred
# 27.97786
grande_diam_cell_pred=cell_coeffs[1]+(cell_coeffs[2]*grande_cell_diam)
grande_diam_cell_pred
# 19.79501
#Petite strain has 8 cells more, which makes sense as bigger cells can overlap more easily and have
#a higher volume to overlap



# Cell synchrony ####


sync_volume_df=read.csv('grande_clust_volume_sync_aprox_2025june12.csv', header=TRUE)

summary(sync_volume_df)

mean_network_volume=sync_volume_df %>%
  group_by(overlap_thresh, file_number) %>%
  summarize(mean_volume=mean(volume_clust),
            mean_cells=mean(num_cells))



ggplot(mean_network_volume, aes(x=overlap_thresh, y=mean_volume))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  scale_x_continuous(trans='log10')+
  geom_hline(yintercept = grande_mean_volume)+
  NULL

ggplot(mean_network_volume, aes(x=overlap_thresh, y=mean_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  scale_x_continuous(trans='log10')+
  # geom_line(data=temp_df, aes(x=temp_seq, y=temp_y), col='blue')+
  NULL

fit_synchrony=glm(mean_network_volume$mean_volume ~ mean_network_volume$overlap_thresh, family = gaussian)
fit_synchrony

b0_sync <- fit_synchrony$coefficients[1]
b1_sync <- fit_synchrony$coefficients[2]

(grande_mean_volume-b0_sync)/b1_sync
#233.1009




#### Size Difference SYNC ####

# sync_grande_diff=read.csv('grande_volume_pred_synchrony_2025june12.csv', header=TRUE)
# sync_petite_diff=read.csv('petite_volume_pred_synchrony_2025june12.csv', header=TRUE)
# sync_df=rbind(sync_grande_diff, sync_petite_diff)
# 
# write.csv(sync_df, file='synchrony_volume_pred_2025june12.csv', row.names=FALSE)

sync_df=read.csv('synchrony_volume_pred_2025june12.csv', header=TRUE)
summary(sync_df)

ggplot(sync_df, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  xlab('Strain')+ylab('Cluster Volume')+
  NULL

ggplot(sync_df, aes(x=strain, y=num_cells, fill=strain))+
  geom_violin()+
  theme_bw()+
  NULL

aggregate(volume_clust ~ strain, data = sync_df, FUN = mean)
# grande  105066.66
# petite  89587.77
#Difference: 15478.89


aggregate(num_cells ~ strain, data = sync_df, FUN = mean)
# grande  652.9434
# petite  613.2746
#Difference: 39.6688



# Combined size difference prediction ####


sync_grande_diff=read.csv('grande_volume_pred_synchrony_2025june12.csv', header=TRUE)
sync_petite_diff=read.csv('petite_combined_effect_volume_pred_synchrony_2025june12.csv', header=TRUE)
combined_effect_df=rbind(sync_grande_diff, sync_petite_diff)


ggplot(combined_effect_df, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  xlab('Strain')+ylab('Cluster Volume')+
  NULL


aggregate(volume_clust ~ strain, data = combined_effect_df, FUN = mean)
# grande    105066.66
# petite     85934.39
# Difference: 19132.27

