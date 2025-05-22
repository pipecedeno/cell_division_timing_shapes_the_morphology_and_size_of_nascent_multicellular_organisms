
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
#(Date)_(strain [gob8/gob21])-(replicate [1-2])_day(1/2)_(image number)_measurements_rois_border.csv

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
  temp_day=as.integer(sub(".*?(\\d+)$", "\\1", metadata[3]))
  temp_image_num=metadata[4]
  temp_replicate=strsplit(temp_strain_replicate, split='-', fixed=TRUE)[[1]][2]
  
  df$date=temp_date
  df$strain=temp_strain
  df$day=temp_day
  df$image_num=temp_image_num
  df$replicate=temp_replicate
  
  # df$full_id=paste(c(metadata[1],metadata[2],metadata[3],metadata[4],metadata[5]), collapse = "_") 
  
  #creating the final data frame
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}

# Loading cluster size data ####
#Directory where the images are located

in_dir="~/work_dir/observed_synchrony/data/Microscopy/ancestor_measurements_22may2024/cluster_area_images/csv_results"

data_clust <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                    .fun = load_one_csv_cluster_size)
data_clust$volume=(4/3)*pi*((data_clust$Major/2)^3)
data_clust$strain=factor(data_clust$strain, levels=c('petite', 'grande'))
summary(data_clust)


# Grande and petite volume measurements ####

summary(data_clust)


median_volume=data_clust %>%
  group_by(strain, day, replicate) %>%
  summarise(mean_vol=mean(volume),
            median_vol=median(volume),
            percentile_25=quantile(volume, c(.25)),
            percentile_75=quantile(volume, c(.75)))

ggplot(median_volume, aes(x=day, y=median_vol, col=strain, group=interaction(strain, replicate)))+
  geom_line()+
  geom_ribbon(data=median_volume, aes(x=day, ymin=percentile_25, ymax=percentile_75, 
                                      fill=strain, group=interaction(strain, replicate)), alpha=0.2, linetype='blank')+
  theme_classic()+
  ylab("Median Cluster Volume (Microns^3)")+
  # scale_y_continuous(trans='log10')+ylab("log10 Cluster Volume (Microns^3)")+
  xlab("Days")+
  NULL


# Cluster Volume
grande_mean_volume=mean(data_clust[data_clust$strain=='grande' & data_clust$day==2,]$volume)
# 103600.1 um^3
petite_mean_volume=mean(data_clust[data_clust$strain=='petite' & data_clust$day==2,]$volume)
# 61110.92 um^3

# Cluster Radius
mean(data_clust[data_clust$strain=='grande' & data_clust$day==2,]$Major)/2
# 27.26306
mean(data_clust[data_clust$strain=='petite' & data_clust$day==2,]$Major)/2
# 23.14429


#this plot has the error that the mean is repeated to be the same each day
ggplot(data_clust, aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("log10 Cluster Volume (Microns^3)")+
  geom_vline(xintercept=grande_mean_volume, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col='#00BFC4', linetype='dashed')+
  scale_x_continuous(trans='log10')+
  facet_wrap(~day, ncol=1)+
  NULL

ggplot(data_clust, aes(x=as.factor(day), y=volume, fill=strain))+
  geom_violin()+
  theme_classic()+
  scale_y_continuous(trans='log10')+
  ylab("log10 Cluster Volume (Microns^3)")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(data_clust, aes(x=as.factor(replicate), y=volume, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_grid(strain~day)+
  scale_y_continuous(trans='log10')+
  ylab("log10 Cluster Volume (Microns^3)")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL


ggplot(data_clust[data_clust$day==2,], aes(x=volume, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("log10 Cluster Volume (Microns^3)")+
  geom_vline(xintercept=grande_mean_volume, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col='#00BFC4', linetype='dashed')+
  scale_x_continuous(trans='log10')+
  facet_wrap(~day, ncol=1)+
  NULL


# Loading Cell Properties Data ####


#format of the csv files
#(Date)_(strain [gob8/gob21])-(replicate [1-2])_day(1/2)_(image number)_measurements_rois_border.csv

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
  temp_day=as.integer(sub(".*?(\\d+)$", "\\1", metadata[3]))
  temp_image_num=metadata[4]
  temp_replicate=strsplit(temp_strain_replicate, split='-', fixed=TRUE)[[1]][2]
  
  df$date=temp_date
  df$strain=temp_strain
  df$day=temp_day
  df$image_num=temp_image_num
  df$replicate=temp_replicate
  
  # df$full_id=paste(c(metadata[1],metadata[2],metadata[3],metadata[4],metadata[5]), collapse = "_") 
  
  #creating the final data frame
  id_vars <- setdiff(colnames(df), measure_vars)
  df <- dplyr::select(df, all_of(c(id_vars, measure_vars)))
  return(df)
}

#Directory where the images are located
in_dir="~/work_dir/observed_synchrony/data/Microscopy/ancestor_measurements_22may2024/calcofluor_images/csv_results"

cell_data <- ldply(.data = list.files(path = in_dir, pattern = "*.csv", full.names = TRUE),
                   .fun = load_one_csv_cell_measurements)

cell_data$strain=factor(cell_data$strain, levels=c('petite', 'grande'))
summary(cell_data)

# Number of cells characterized in the second day
table(cell_data[cell_data$day==2,]$replicate, cell_data[cell_data$day==2,]$strain)
#     grande petite
# 1   3113   1052
# 2   2409   1100
# 3   2466   1428
# 4   2899   1392
# 5   2668   2107

#### Aspect Ratio ####

median_ar=cell_data %>%
  group_by(strain, day, replicate) %>%
  summarise(mean_ar=mean(AR),
            median_ar=median(AR),
            percentile_25=quantile(AR, c(.25)),
            percentile_75=quantile(AR, c(.75)))

ggplot(median_ar, aes(x=day, y=median_ar, col=strain, group=interaction(strain, replicate)))+
  geom_line()+
  geom_ribbon(data=median_ar, aes(x=day, ymin=percentile_25, ymax=percentile_75, 
                                  fill=strain, group=interaction(strain, replicate)), alpha=0.2, linetype='blank')+
  theme_classic()+
  ylab("Median Cell Aspect Ratio")+
  # scale_y_continuous(trans='log10')+ylab("log10 Cluster Volume (Microns^3)")+
  xlab("Days")+
  NULL


grande_mean_ar=mean(cell_data[cell_data$strain=='grande' & cell_data$day==2,]$AR)
grande_mean_ar
# mean: 1.199049
# median: 1.19058
petite_mean_ar=mean(cell_data[cell_data$strain=='petite' & cell_data$day==2,]$AR)
petite_mean_ar
# mean: 1.22081
# median: 1.20949

ggplot(cell_data, aes(x=as.factor(day), y=AR, fill=strain))+
  geom_violin()+
  theme_classic()+
  ylab("Cell Aspect Ratio")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data, aes(x=as.factor(replicate), y=AR, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_grid(strain~day)+
  ylab("Cell Aspect Ratio")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data[cell_data$day==2,], 
       aes(x=AR, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("Cell Aspect Ratio")+
  geom_vline(xintercept=grande_mean_ar, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_ar, col='#00BFC4', linetype='dashed')+
  # scale_x_continuous(trans='log10')+
  # facet_wrap(~day, ncol=1)+
  facet_wrap(~strain, ncol=1)+
  NULL


#### Cell Diameter ####

summ_diam=cell_data %>%
  group_by(strain, day, replicate) %>%
  summarise(mean_diam=mean(Major),
            median_diam=median(Major),
            percentile_25=quantile(Major, c(.25)),
            percentile_75=quantile(Major, c(.75)))

ggplot(summ_diam, aes(x=day, y=median_diam, col=strain, group=interaction(strain, replicate)))+
  geom_line()+
  geom_ribbon(data=summ_diam, aes(x=day, ymin=percentile_25, ymax=percentile_75, 
                                  fill=strain, group=interaction(strain, replicate)), alpha=0.2, linetype='blank')+
  theme_classic()+
  ylab("Median Cell Diameter")+
  xlab("Days")+
  NULL


grande_mean_diam=mean(cell_data[cell_data$strain=='grande' & cell_data$day==2,]$Major)
grande_mean_diam
# mean: 5.037635 um
# median: 5.01878 um
petite_mean_diam=mean(cell_data[cell_data$strain=='petite' & cell_data$day==2,]$Major)
petite_mean_diam
# mean: 4.844438 um
# median: 4.80983 um

ggplot(cell_data, aes(x=as.factor(day), y=Major, fill=strain))+
  geom_violin()+
  theme_classic()+
  ylab("Cell Diameter")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data, aes(x=as.factor(replicate), y=Major, fill=strain))+
  geom_violin()+
  theme_classic()+
  facet_grid(strain~day)+
  ylab("Cell Diameter")+
  xlab("Days")+
  stat_summary(fun='median', geom='crossbar', position=position_dodge(width=0.9), width=0.5)+
  NULL

ggplot(cell_data[cell_data$day==2,], 
       aes(x=Major, group=interaction(strain, replicate)))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic()+
  guides(alpha = "none")+
  xlab("Cell Diameter")+
  geom_vline(xintercept=grande_mean_diam, col='#F8766D', linetype='dashed')+
  geom_vline(xintercept=petite_mean_diam, col='#00BFC4', linetype='dashed')+
  # scale_x_continuous(trans='log10')+
  # facet_wrap(~day, ncol=1)+
  facet_wrap(~strain, ncol=1)+
  NULL


# Paper Figure ####

p_cluster_vol_log=ggplot(data_clust[data_clust$day==2,], aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
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

p_cluster_vol=ggplot(data_clust[data_clust$day==2,], aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic(base_size = 14)+
  guides(alpha = "none")+
  ylab(NULL)+
  xlab(expression("Cluster Volume"~(mu*m^3)))+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF"))+
  geom_vline(xintercept=grande_mean_volume, col="#AE93BEFF", linetype='dashed')+
  geom_vline(xintercept=petite_mean_volume, col="#F0D77BFF", linetype='dashed')+
  xlim(c(0, 250000))+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.line.y = element_blank()) +  # Remove y-axis line
  labs(fill = "Strain")+
  NULL
p_cluster_vol

p_ar=ggplot(cell_data[cell_data$day==2,], 
       aes(x=AR))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic(base_size = 14)+
  guides(alpha = "none", fill="none")+
  ylab(NULL)+
  theme(axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.line.y = element_blank()) +  # Remove y-axis line
  xlab("Aspect Ratio")+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF")) +
  geom_vline(xintercept=grande_mean_ar, col="#AE93BEFF", linetype='dashed')+
  geom_vline(xintercept=petite_mean_ar, col="#F0D77BFF", linetype='dashed')+
  NULL
p_ar


p_diam=ggplot(cell_data[cell_data$day==2,], 
            aes(x=Major))+
  geom_density(aes(fill=strain, alpha=0.5))+
  theme_classic(base_size = 14)+
  guides(alpha = "none", fill="none")+
  ylab(NULL)+
  theme(axis.text.y = element_blank(),  # Remove y-axis text
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.line.y = element_blank()) +  # Remove y-axis line
  xlab(expression("Cell Diameter"~(mu*m)))+
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

row_1=plot_grid(p_cluster_vol, col_clust_vol, 
                labels=c('A', 'B'), nrow=1, label_size=12, rel_widths=c(1, 0.5))


fig_1=plot_grid(p_cluster_vol, col_clust_vol, p_ar, col_ar, 
                labels=c('A', 'B','C','D'), nrow=2, label_size=16, rel_widths=c(1, 0.35), rel_heights=c(1,1))
fig_1


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_1_clust_size_ancestors_15may2025.svg',
#        plot=fig_1, dpi='retina', width=8, height=6, bg='white')
# note: image is saved as an svg to later add the strain labels for the images


#### version 2 ####

p_cluster_vol_log=ggplot(data_clust[data_clust$day==2,], aes(x=volume))+
  geom_density(aes(fill=strain, alpha=0.5))+
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

p_cluster_vol=ggplot(data_clust[data_clust$day==2,], aes(x=volume, fill=strain))+
  # geom_density(aes(fill=strain, alpha=0.5))+
  geom_histogram(position="identity", alpha=0.5, bins=30) +
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

p_ar=ggplot(cell_data[cell_data$day==2,], 
            aes(x=AR, fill=strain))+
  geom_histogram(alpha=0.5, position="identity", bins=30) +
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


p_diam=ggplot(cell_data[cell_data$day==2,], 
              aes(x=Major, fill=strain))+
  geom_histogram(alpha=0.5, position="identity", bins=30) +
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

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_1_clust_size_ancestors_15may2025_v2.svg',
       plot=fig_1_v2, dpi='retina', width=8, height=6, bg='white')
# note: image is saved as an svg to later add the strain labels for the images
