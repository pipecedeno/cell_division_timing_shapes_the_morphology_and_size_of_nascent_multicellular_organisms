

# Date 26Aug2024
# This code was made to analyze the results of using synthetic data to analyze the effects
# of the first division being faster or slower than the rest of the divisions


library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(ghibli)
library(png)
library(grid)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$YesterdayMedium)
syn_data_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("YesterdayMedium", direction = -1))


setwd("~/work_dir/observed_synchrony/evolution_results/fast_first_division_26aug2024/")

# Doubling time distributions ####

doub_t_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_8_fast_first_division/fast_first_div_doub_times_26aug2024.csv", header=TRUE)

doub_t_df$experiment_name=factor(doub_t_df$experiment_name, levels=c('Slow First Division', 'Synchronous',
                                                                     'Fast First Division'))

summary(doub_t_df)

ggplot(doub_t_df[doub_t_df$number_divisions<=2,], 
       aes(x=as.factor(number_divisions), y=minutes, fill=as.factor(number_divisions)))+
  geom_violin()+
  facet_wrap(~experiment_name)+
  xlab('Number of Divisions')+
  ylab("Minutes")+
  guides(fill = guide_legend(title = "div_num"))+
  NULL


# Difference in doubling time ####

diff_dt_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_8_fast_first_division/fast_first_div_difference_doub_times_26aug2024.csv", header=TRUE)

diff_dt_df$experiment_name=factor(diff_dt_df$experiment_name, levels=c('Slow First Division', 'Synchronous',
                                                                     'Fast First Division'))
summary(diff_dt_df)

#this plot uses all the data from all the generations so that is why it takes a long time to run
ggplot(data=diff_dt_df, 
       aes(x=percent_cell_cycle, col=experiment_name))+
  stat_ecdf(geom='step')+
  theme_minimal()+
  ylab("Fraction of mother-daughter cells that has divided")+
  xlab("Cell cycles")+
  syn_data_colors+
  NULL



# Fragmentation summary ####

frag_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_8_fast_first_division/fast_first_div_frag_df_26aug2024.csv", header=TRUE)

frag_df$experiment_name=factor(frag_df$experiment_name, levels=c('Slow First Division', 'Synchronous',
                                                                 'Fast First Division'))
summary(frag_df)

table(frag_df$generation)


#### Size at fracture ####

ggplot(frag_df[frag_df$generation%%10==0,],
       aes(x=as.factor(generation), y=cluster_size, fill=experiment_name))+
  geom_violin()+
  theme_bw()+
  ylab("Cluster size at fracture")+xlab("Generations")+
  syn_data_colors+
  NULL

mean_clust_size=frag_df %>%
  group_by(experiment_name, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))


ggplot(mean_clust_size, aes(x=generation, y=mean_size, col=experiment_name))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=experiment_name), 
              alpha=0.2, linetype='blank')+
  ggtitle("Mean Fracture Size")+
  ylab("Mean fracure size")+xlab('Generation')+
  syn_data_colors+
  NULL


ggplot(mean_clust_size, aes(x=experiment_name, y=mean_size, fill=experiment_name))+
  geom_violin()+
  # stat_summary(fun='mean', geom='crossbar')+
  syn_data_colors+
  guides(fill='none')+
  labs(x='Synthetic Strain',
       y='Mean Fracture Size')+
  NULL



#### Fracture proportion ####

summary(frag_df)

mean_frag_prop=frag_df %>%
  group_by(experiment_name, generation) %>%
  summarise(mean_propagule=mean(proportion_propagule),
            sd_propagule=sd(proportion_propagule))


ggplot(mean_frag_prop)+
  theme_bw()+
  geom_line(aes(x=generation, y=mean_propagule, col=experiment_name))+
  geom_ribbon(aes(x=generation, y=mean_propagule, ymin=mean_propagule-sd_propagule, ymax=mean_propagule+sd_propagule, fill=experiment_name), 
              alpha=0.2, linetype='blank')+
  ggtitle('Propagule proportion after fracture')+
  ylab("Mean fracture proportion")+xlab('Generation')+
  syn_data_colors+
  NULL


ggplot(frag_df[frag_df$generation%%10==0,], 
       aes(x=as.factor(generation), y=proportion_propagule, fill=experiment_name))+
  geom_violin()+
  theme_bw()+
  xlab("Generation")+ylab("Proportion offspring size")+
  syn_data_colors+
  NULL


ggplot(mean_frag_prop, aes(x=experiment_name, y=mean_propagule, fill=experiment_name))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  guides(fill='none')+
  labs(x='Strain',
       y='Mean Propagule Proportion')+
  syn_data_colors+
  NULL


# Network Diameter ####

diameter_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_8_fast_first_division/fast_first_div_diameter_df_26aug2024.csv", header=TRUE)

diameter_df$experiment_name=factor(diameter_df$experiment_name, levels=c('Slow First Division', 'Synchronous',
                                                                 'Fast First Division'))
summary(diameter_df)

summary(frag_df)

#To normalize the network diameter metric I need to first join the diameter and size at fracture information
#to later divide network diameter by cluster size

diam_size <- left_join(frag_df, diameter_df, 
                       by = c("sim_number", "generation", "experiment_name"))

summary(diam_size)

# diam_size$norm_diameter=diam_size$diameter/diam_size$cluster_size
diam_size$norm_diameter=diam_size$diameter/(6.64385619*log10(diam_size$cluster_size)-1)


# Plots of diameter by generation
summ_diam=diam_size %>%
  group_by(experiment_name, generation) %>%
  summarise(mean_diam=mean(norm_diameter),
            sd_diam=sd(norm_diameter),
            median_diam=median(norm_diameter))



ggplot(summ_diam, aes(x=generation, y=mean_diam, col=experiment_name))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_diam, ymin=mean_diam-sd_diam, ymax=mean_diam+sd_diam, fill=experiment_name), 
              alpha=0.2, linetype='blank')+
  ggtitle("Mean of Normalized Network Diameter")+
  xlab("Generation")+
  ylab("Normalized network diameter")+
  syn_data_colors+
  NULL

ggplot(summ_diam, aes(x=experiment_name, y=mean_diam, fill=experiment_name))+
  geom_violin()+
  guides(fill='none')+
  labs(x='Strain',
       y='Mean Normalized Network Diameter')+
  syn_data_colors+
  NULL


# Paper Figures ####


fig_a=ggplot(doub_t_df[doub_t_df$number_divisions<=2,], 
       aes(x=as.factor(number_divisions), y=minutes, fill=as.factor(number_divisions)))+
  geom_violin()+
  stat_summary(geom='crossbar', fun='mean')+
  facet_wrap(~experiment_name)+
  xlab('Number of Divisions')+
  ylab("Minutes")+
  guides(fill = 'none')+
  theme_classic(base_size = 10)+
  NULL
fig_a


fig_b=ggplot(data=diff_dt_df, 
             aes(x=percent_cell_cycle, col=experiment_name))+
  stat_ecdf(geom='step')+
  theme_minimal()+
  ylab("Fraction of mother-daughter\ncells that have divided")+
  xlab("Cell cycles")+
  syn_data_colors+
  guides(col = guide_legend(title = "Strain"))+
  guides(col='none')+
  theme_minimal(base_size = 10)+
  NULL
fig_b


fig_c=ggplot(mean_clust_size, aes(x=experiment_name, y=mean_size, fill=experiment_name))+
  geom_violin()+
  # stat_summary(fun='mean', geom='crossbar')+
  syn_data_colors+
  guides(fill='none')+
  labs(x='Strain',
       y='Mean Fracture Size')+
  theme_classic(base_size = 10)+
  NULL
fig_c

mean_clust_size %>%
  group_by(experiment_name) %>%
  summarise(mean_frac_size=mean(mean_size))
#   experiment_name         mean_frac_size
# 1 Slow First Division           120.8874
# 2 Synchronous                   216.1338
# 3 Fast First Division           542.2014


fig_d=ggplot(mean_frag_prop, aes(x=experiment_name, y=mean_propagule, fill=experiment_name))+
  geom_violin()+
  # stat_summary(fun='mean', geom='crossbar')+
  guides(fill='none')+
  labs(x='Strain',
       y='Mean Propagule Proportion')+
  syn_data_colors+
  theme_classic(base_size = 10)+
  NULL
fig_d


fig_e=ggplot(summ_diam, aes(x=experiment_name, y=mean_diam, fill=experiment_name))+
  geom_violin()+
  guides(fill='none')+
  labs(x='Strain',
       y='Mean Normalized\nNetwork Diameter')+
  syn_data_colors+
  theme_classic(base_size = 10)+
  NULL
fig_e


summ_diam %>%
  group_by(experiment_name) %>%
  summarise(mean_norm_diam=mean(mean_diam))
#   experiment_name     mean_norm_diam
# 1 Slow First Division          0.872
# 2 Synchronous                  0.990
# 3 Fast First Division          1.19 





img_sync <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/synchronous_div_filament_size_3_3oct2024.png")
img_plot_sync_net <- rasterGrob(img_sync, interpolate = TRUE)

img_fast_first <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/fast_first_div_filament_size_3_3oct2024.png")
img_plot_fast_first_net <- rasterGrob(img_fast_first, interpolate = TRUE)

# Create text annotations
text_sync <- textGrob("Synchronous", gp = gpar(fontsize = 10, fontface = "bold"))
text_fast_first <- textGrob("Fast First Division", gp = gpar(fontsize = 10, fontface = "bold"))

# Create ggplot objects for the images with annotations
p_sync_net <- ggplot() + 
  annotation_custom(img_plot_sync_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_sync, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_sync_net

p_fast_first <- ggplot() + 
  annotation_custom(img_plot_fast_first_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_fast_first, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_fast_first



figure_first_div_net=plot_grid(fig_a, fig_b, fig_c, fig_e, p_fast_first, p_sync_net,
                                 labels=c('A', 'B', 'C', 'D', 'E', 'F'), ncol=2, 
                                 align='hv', label_size=12)
figure_first_div_net

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_8_first_div_23apr2025_networks.png',
       plot=figure_first_div_net, dpi='retina', height=10, width=8, bg='white')


