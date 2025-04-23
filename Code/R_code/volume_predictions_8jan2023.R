
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

setwd('~/work_dir/observed_synchrony/physics_model/overlap_volume_prediction_15dec2023/')


# PA ancestors measurements ####
# This code is from test_ypd_conc_effect_30jun2023.R

load_one_csv <- function(file) {
  df <- read.csv(file, row.names = 1) #reading the csv file
  measure_vars <- colnames(df) #getting the names of the columns
  filename <- basename(file) #saving only the name of the file, so not saving the path to the file
  
  #removing the last part of the name of the file "_measurements_rois.csv" to only keep the important information
  
  metadata <- strsplit(filename, split = "_", fixed = TRUE)[[1]] #splitting the metadata into the parts of the information
  temp_date=metadata[1] 
  if(metadata[2]=='gob21'){
    temp_strain='petite'
  } else {
    temp_strain='grande'
  }
  temp_treated=metadata[3]
  temp_media_conc=metadata[4]
  temp_dilution=metadata[5]
  temp_image_num=metadata[6]
  
  df$date=temp_date
  df$strain=temp_strain
  df$treatment=paste(temp_treated, temp_media_conc, collapse="_")
  df$dilution=temp_dilution
  df$image_num=temp_image_num
  
  df$full_id=paste(c(metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6]), collapse = "_") 
  
  #creating the final data frame
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}


in_dir="~/work_dir/observed_synchrony/data/Microscopy/media_concentration_27-28jun2023"

data_clust <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                    .fun = load_one_csv)
data_clust$volume=(4/3)*pi*((data_clust$Major/2)^3)
summary(data_clust)

data_clust$treatment=factor(data_clust$treatment, levels = c("before 1x","after 0.5x", "after 1x",
                                                             "after 1.5x", "after 2x"))

table(data_clust$treatment)


ggplot(data_clust, aes(x=strain, y=Major, fill=treatment))+
  geom_violin()+
  ylab("Diameter")+
  NULL

ggplot(data_clust[data_clust$treatment=='after 1x',],
       aes(x=strain, y=Major, fill=strain))+
  geom_violin()+
  theme_bw()+
  stat_summary(fun = "mean",
               geom = "crossbar")+
  NULL

ggplot(data_clust[data_clust$treatment=='after 1x',], aes(x=Major))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("Cluster Diameter (Microns)")+
  geom_vline(xintercept=52.9, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=43.51, col='#00BFC4', linetype='dashed')+
  NULL

aggregate(Major/2 ~ strain, data = data_clust[data_clust$treatment=='after 1x',], FUN = mean)


grande_mean_volume=mean(data_clust[data_clust$treatment=='after 1x' 
                              & data_clust$strain=='grande',]$volume)
petite_mean_volume=mean(data_clust[data_clust$treatment=='after 1x' 
                                   & data_clust$strain=='petite',]$volume)

ggplot(data_clust[data_clust$treatment=='after 1x',], aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("log10 Cluster Volume (Microns^3)")+
  # geom_vline(xintercept=77951.81, col='#F8766D', linetype='dashed')+
  # geom_vline(xintercept=43039.50111, col='#00BFC4', linetype='dashed')+
  geom_vline(xintercept=grande_mean_volume, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col='#00BFC4', linetype='dashed')+
  scale_x_continuous(trans='log10')+
  NULL

aggregate(volume ~ strain, data = data_clust[data_clust$treatment=='after 1x',], FUN = mean)

grande_mean_volume-petite_mean_volume
# 40291.29

# Aspect ratio ####
#### Threshold prediction ####

# aspr_volume_df=read.csv("exponential_clust_volume_1.256aspr_2024jan8.csv", header=TRUE)
aspr_volume_df=read.csv("exponential_clust_volume_1.256aspr_more_overlaps_v2.csv", header=TRUE)

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
# 60.98398 #making it check the overlap after each cell added



#### Size difference AR ####

aspr_df=read.csv("constant_threshold_changing_aspr_2024jan10.csv", header=TRUE)

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

petite_asp_r=1.259
grande_asp_r=1.256

vol_coeffs[1]+(vol_coeffs[2]*petite_asp_r)+(vol_coeffs[3]*(petite_asp_r^2))
# 58296.85
vol_coeffs[1]+(vol_coeffs[2]*grande_asp_r)+(vol_coeffs[3]*(grande_asp_r^2))
# 57649.16
# Petite strain has 647.69 um^3 more

# Expected difference in number of cells
cells_model=glm(aspr_df$num_cells ~ poly(aspr_df$aspect_ratio, 2, raw=TRUE), family=gaussian)
cell_coeffs=cells_model$coefficients

petite_asp_r=1.259
grande_asp_r=1.256

cell_coeffs[1]+(cell_coeffs[2]*petite_asp_r)+(cell_coeffs[3]*(petite_asp_r^2))
# 351.9224
cell_coeffs[1]+(cell_coeffs[2]*grande_asp_r)+(cell_coeffs[3]*(grande_asp_r^2))
# 350.5404
#Petite strain has 1 cell more


#### Plot close to petite and grande values ####

# table(aspr_df$aspect_ratio)
# #1.25 petite
# #1.2 grande
# 
# temp_grande_ar=aspr_df[aspr_df$aspect_ratio==1.2,]
# temp_grande_ar$strain='grande'
# temp_petite_ar=aspr_df[aspr_df$aspect_ratio==1.25,]
# temp_petite_ar$strain='petite'
# ar_strains=rbind(temp_grande_ar, temp_petite_ar)
# 
# ar_strains[, 1:4] <- lapply(ar_strains[, 1:4], as.numeric)
# summary(ar_strains)
# 
# ggplot(ar_strains, aes(x=strain, y=volume_clust, fill=strain))+
#   geom_violin()+
#   theme_bw()+
#   NULL

temp_petite_ar=read.csv('petite_volume_prediction_1.259ar.csv', header=TRUE)
temp_petite_ar$strain='petite'
temp_grande_ar=read.csv('grande_volume_prediction_1.256ar.csv', header=TRUE)
temp_grande_ar$strain='grande'
ar_strains=rbind(temp_grande_ar, temp_petite_ar)

ar_strains[, 1:4] <- lapply(ar_strains[, 1:4], as.numeric)
summary(ar_strains)

ggplot(ar_strains, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  NULL


#### Recreating aspect ratio volume plot ####
#Note: These plots are using the mean radius, which is not the same variable which
#is used in the article, as in the article the weighted mean radius is used instead


pa_data=read.csv('PA_cluster_radius_and_volume.csv', header=TRUE)

summary(pa_data)

ancestors_data=data.frame(matrix(c('petite', 1.259, 43039.50111,
                                   'grande', 1.256, 77951.81), ncol=3, byrow=TRUE))
colnames(ancestors_data)=c('strain', 'mean_aspect_ratio', 'mean_volume')
ancestors_data$mean_aspect_ratio=as.numeric(ancestors_data$mean_aspect_ratio)
ancestors_data$mean_volume=as.numeric(ancestors_data$mean_volume)
summary(ancestors_data)

ggplot(aspr_df, aes(x=aspect_ratio, y=volume_clust))+
  geom_point(alpha=0.1, color='blue', shape='.')+
  # geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown', se=TRUE, fill='yellow') +
  theme_bw()+
  geom_point(data=pa_data, aes(x=mean_aspect_ratio, y=mean_volume, col=strain))+
  # geom_point(data=pa_data[pa_data$time_point<=300,], 
  #            aes(x=mean_aspect_ratio, y=mean_volume, col=strain))+
  geom_point(data=ancestors_data, aes(x=mean_aspect_ratio, y=mean_volume))+
  ylab('Cluster Volume')+xlab('Cellular Aspect Ratio')+
  scale_y_continuous(trans='log10')+ylab('Log10(Cluster Volume)')+
  NULL

ggplot(aspr_df[aspr_df$aspect_ratio<2.5 & aspr_df$aspect_ratio>1,], 
       aes(x=aspect_ratio, y=volume_clust))+
  geom_point(alpha=0.1, color='blue', shape='.')+
  # geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown', se=TRUE, fill='yellow') +
  theme_bw()+
  geom_point(data=pa_data[pa_data$time_point<=250,], 
             aes(x=mean_aspect_ratio, y=mean_volume, col=strain))+
  geom_point(data=ancestors_data, aes(x=mean_aspect_ratio, y=mean_volume, shape=strain))+
  ylab('Cluster Volume')+xlab('Cellular Aspect Ratio')+
  NULL

# Using radius as in the original plot

aspr_df$radius=((3*aspr_df$volume_clust)/(4*pi))^(1/3)
ancestors_data$radius=((3*ancestors_data$mean_volume)/(4*pi))^(1/3)

ggplot(aspr_df[aspr_df$aspect_ratio<2.5 & aspr_df$aspect_ratio>1,], 
       aes(x=aspect_ratio, y=radius))+
  geom_point(alpha=0.1, color='blue', shape='.')+
  geom_smooth(method=lm , color="brown", fill="brown", se=TRUE)+
  # stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown', se=TRUE, fill='yellow') +
  theme_bw()+
  # geom_point(data=pa_data, aes(x=mean_aspect_ratio, y=mean_volume, col=strain))+
  geom_point(data=pa_data[pa_data$time_point<=300,],
             aes(x=mean_aspect_ratio, y=mean_radius, col=strain))+
  geom_point(data=pa_data[pa_data$mean_aspect_ratio<2.4,],
             aes(x=mean_aspect_ratio, y=mean_radius, col=strain))+
  geom_point(data=ancestors_data, aes(x=mean_aspect_ratio, y=radius))+
  ylab('Cluster Radius')+xlab('Cellular Aspect Ratio')+
  # scale_y_continuous(trans='log10')+ylab('Log10(Cluster Volume)')+
  ylim(c(15,80))+
  NULL


#### Using weighted mean radius ####

weighted_pa_df=read.csv("fig2d_raw_data.csv", header=TRUE)

summary(weighted_pa_df)

aspr_df$radius=((3*aspr_df$volume_clust)/(4*pi))^(1/3)
ancestors_data$radius=((3*ancestors_data$mean_volume)/(4*pi))^(1/3)

ggplot(aspr_df[aspr_df$aspect_ratio<2.5 & aspr_df$aspect_ratio>1,], 
       aes(x=aspect_ratio, y=radius))+
  geom_point(alpha=0.1, color='blue', shape='.')+
  # geom_smooth(method=lm , color="black", fill="red", se=TRUE)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown', se=TRUE, fill='yellow') +
  theme_bw()+
  geom_point(data=weighted_pa_df[weighted_pa_df$Aspect_ratio<=2.4,],
             aes(x=Aspect_ratio, y=Cluster_radius_.µm., col=Population))+
  geom_point(data=ancestors_data, aes(x=mean_aspect_ratio, y=radius))+
  ylab('Cluster Radius')+xlab('Cellular Aspect Ratio')+
  # scale_y_continuous(trans='log10')+ylab('Log10(Cluster Volume)')+
  ylim(c(15,80))+
  NULL

weighted_pa_df$Cluster_radius_.µm.

# Cell Diameter ####
#### Threshold prediction ####

diam_volume_df=read.csv("exponential_clust_volume_13.59cell_diam_2024jan8.csv", header=TRUE)

summary(diam_volume_df)

mean_volume=diam_volume_df %>%
  group_by(overlap_thresh) %>%
  summarize(mean_volume=mean(volume_clust),
            mean_num_cells=mean(num_cells))

ggplot(diam_volume_df, aes(x=overlap_thresh, y=volume_clust))+
  geom_point(alpha=0.1, color='blue')+
  geom_line(data=mean_volume, aes(x=overlap_thresh, y=mean_volume))+
  theme_bw()+
  scale_x_continuous(trans='log10')+
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
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  scale_x_continuous(trans='log10')+
  geom_hline(yintercept = grande_mean_volume)+
  NULL

ggplot(diam_volume_df, aes(x=overlap_thresh, y=num_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  theme_bw()+
  scale_x_continuous(trans='log10')+
  NULL

coeffs=glm(diam_volume_df$volume_clust ~ poly(diam_volume_df$overlap_thresh, 2, raw=TRUE), family=gaussian)
coeffs

#The overlap that I need to use is of 1.17975 by finding the root in geogebra

table(diam_volume_df$overlap_thresh)
mean(diam_volume_df[diam_volume_df$overlap_thresh==1.122018,]$volume_clust)


#### Size Difference DIAM ####

diam_df=read.csv("constant_threshold_changing_cell_diam_2024jan10.csv", header=TRUE)

summary(diam_df)

petite_cell_diam=7.38
grande_cell_diam=13.59

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
# 9905.454 
grande_diam_pred=vol_coeffs[1]+(vol_coeffs[2]*grande_cell_diam)
grande_diam_pred
# 34381.15
grande_diam_pred-petite_diam_pred
# Grande strain  has 24475.7 um^3 more

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

#### Plot close to petite and grande values ####

table(diam_df$cell_diam)
#7.4 petite
#13.6 grande

temp_grande_diam=diam_df[diam_df$cell_diam==13.6,]
temp_grande_diam$strain='grande'
temp_petite_diam=diam_df[diam_df$cell_diam==7.4,]
temp_petite_diam$strain='petite'
diam_strains=rbind(temp_grande_diam, temp_petite_diam)

diam_strains[, 1:4] <- lapply(diam_strains[, 1:4], as.numeric)
summary(diam_strains)

ggplot(diam_strains, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  NULL

# Cell synchrony ####


sync_volume_df=read.csv('grande_clust_volume_2024jan8.csv', header=TRUE)

summary(sync_volume_df)

mean_network_volume=sync_volume_df %>%
  group_by(overlap_thresh, file_number) %>%
  summarize(mean_volume=mean(volume_clust),
            mean_cells=mean(num_cells))

temp_seq=seq(10,1000, 1)
temp_y=-0.1405 * (temp_seq)^2 + 256.1017 * temp_seq + 22955.5576
temp_df=data.frame(temp_seq, temp_y)
summary(temp_df)

ggplot(mean_network_volume, aes(x=overlap_thresh, y=mean_volume))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  scale_x_continuous(trans='log10')+
  # geom_line(data=temp_df, aes(x=temp_seq, y=temp_y), col='blue')+
  geom_hline(yintercept = grande_mean_volume)+
  NULL

ggplot(mean_network_volume, aes(x=overlap_thresh, y=mean_cells))+
  geom_point(alpha=0.1, color='blue')+
  geom_smooth(method=lm , color="black", fill="red", se=TRUE, linewidth=0.2)+
  stat_smooth(method='lm', formula = y ~ poly(x,2), linewidth = 1, color='brown') +
  theme_bw()+
  scale_x_continuous(trans='log10')+
  # geom_line(data=temp_df, aes(x=temp_seq, y=temp_y), col='blue')+
  NULL

fit_synchrony=glm(mean_network_volume$mean_volume ~ mean_network_volume$overlap_thresh, family = gaussian)
fit_synchrony

b0_sync <- fit_synchrony$coefficients[1]
b1_sync <- fit_synchrony$coefficients[2]

(grande_mean_volume-b0_sync)/b1_sync
#468.9796

glm(mean_network_volume$mean_volume ~ poly(mean_network_volume$overlap_thresh, 2, raw=TRUE), family=gaussian)
#333.0398


#### Size Difference SYNC ####

sync_df=read.csv('synchrony_volume_pred_2024jan10.csv', header=TRUE)
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
#grande: 93131.86
#petite: 79931.94
#Difference: 13199.92


aggregate(num_cells ~ strain, data = sync_df, FUN = mean)
# grande  834.1804
# petite  788.2379
#Difference: 45.9425


# Plot 3 variables together ####

p1=ggplot(ar_strains, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  xlab('Strain')+ylab('Cluster Volume')+
  stat_summary(fun = "mean",
               geom = "crossbar")+
  ggtitle('Aspect Ratio')+
  NULL
p1

p2=ggplot(diam_strains, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  xlab('Strain')+ylab('Cluster Volume')+
  theme(axis.title.y=element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar")+
  ggtitle('Cell Diameter')+
  NULL
p2

p3=ggplot(sync_df, aes(x=strain, y=volume_clust, fill=strain))+
  geom_violin()+
  theme_bw()+
  xlab('Strain')+ylab('Cluster Volume')+
  theme(axis.title.y=element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar")+
  ggtitle('Cell Synchrony')+
  NULL
p3

ggarrange(p1,p2,p3, legend='none', ncol=3)

