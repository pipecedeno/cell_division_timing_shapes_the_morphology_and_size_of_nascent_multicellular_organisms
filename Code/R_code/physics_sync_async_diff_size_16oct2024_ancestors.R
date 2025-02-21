
#Date: 16 October 2024
# This code is used to analyze the results that are in:
# work_dir/observed_synchrony/evolution_results/physics_sim_26sep2024/clust_size_sync_vs_async_26sep2024
# the idea is to show that there is a difference in fragmentation sizes and how the strain is accumulated
# between sync and async divisions
# this version is the comparison in size between the petite and grande ancestors


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

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))


# Cluster size difference ####

clust_size=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/supp_fig_3_size_difference_physics_sim/physics_sim_frag_size_sync_vs_async_2024dec3.csv", header = TRUE)
clust_size$strain=factor(clust_size$strain, levels=c('Petite', 'Grande'))


summary(clust_size)

summ_clust_size=clust_size %>%
  group_by(strain, sim_number) %>%
  summarise(mean_size=mean(size_fracture),
            mean_volume=mean(volume),
            mean_diam=mean(gyration_diam))

ggplot(summ_clust_size, aes(x=strain, y=mean_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Mean Fracture Size (Number of Cells)')+
  NULL

ggplot(summ_clust_size, aes(x=strain, y=mean_volume, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Volume (μm^3)')+
  NULL

ggplot(summ_clust_size, aes(x=strain, y=mean_diam, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Cluster Diameter (μm)')+
  NULL

## All of the plots look the same

summ_clust_size %>%
  group_by(strain) %>%
  summarise(mean_fracture_size=mean(mean_size),
            mean_cluster_volume=mean(mean_volume),
            mean_gyration=mean(mean_diam))
#overlap threshold 30
#   strain     mean_fracture_size   mean_cluster_volume     mean_gyration
# 1 Petite               189.1327              25072.59          27.99673
# 2 Grande               207.3925              32022.67          29.86071
# 11 cell difference between grande and petite

# Strain accumulation ####


petite_overlap=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/supp_fig_3_size_difference_physics_sim/petite_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
petite_overlap$strain="Petite"
grande_overlap=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/supp_fig_3_size_difference_physics_sim/grande_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
grande_overlap$strain="Grande"
overlap_size_df=rbind(grande_overlap, petite_overlap)

overlap_size_df$strain=factor(overlap_size_df$strain, levels=c('Petite', 'Grande'))

summary(overlap_size_df)

summ_overlap_size=overlap_size_df %>%
  group_by(strain, file_num, cluster_size) %>%
  summarise(mean_overlap_vol=mean(overlap_vol),
            mean_not_added_nodes=mean(not_added_nodes))

ggplot(summ_overlap_size, 
       aes(x=cluster_size, y=mean_overlap_vol, col=strain, group=interaction(strain, file_num)))+
  geom_line()+
  petite_t200_colors+
  NULL


mean_summ_overlap_size = summ_overlap_size %>%
  group_by(strain, cluster_size) %>%
  summarise(mean_overlap = mean(mean_overlap_vol),
            mean_not_added=mean(mean_not_added_nodes),
            sd_overlap = sd(mean_overlap_vol),
            n = n(),  # Count number of observations
            sem = sd_overlap / sqrt(n),  # Standard error of mean
            ci_margin = qt(0.975, df = n-1) * sem,  # 95% CI margin of error
            ci_lower = mean_overlap - ci_margin,  # Lower bound of CI
            ci_upper = mean_overlap + ci_margin)   # Upper bound of CI

ggplot(mean_summ_overlap_size, 
       aes(x=cluster_size, y=mean_overlap, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=cluster_size, y=mean_overlap, ymin=ci_lower,
                  ymax=ci_upper, fill=strain), alpha=0.2, linetype='blank')+
  geom_hline(yintercept=30, linetype="dashed", color = "black")+
  petite_t200_colors+
  theme_bw()+
  xlab('Cluster size (number of cells)')+
  ylab('Mean Overlap')+
  NULL

mean_summ_overlap_size[mean_summ_overlap_size$cluster_size==200,]


# Overlap Positions ####

overlap_pos_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/supp_fig_3_size_difference_physics_sim/stats_overlap_pos_9dec2024.csv", header=TRUE)


overlap_pos_df$strain=factor(overlap_pos_df$strain, levels=c("Petite", "Grande"))
summary(overlap_pos_df)

tapply(overlap_pos_df$sd_pairwise_dist, overlap_pos_df$strain, summary)

summ_overlap_pos=overlap_pos_df %>%
  group_by(strain, file_num) %>%
  summarise(mean_dist=mean(avg_pairwise_dist),
            mean_hopkins=mean(hopkins_stat))

table(summ_overlap_pos$strain)

ggplot(summ_overlap_pos, aes(x=strain, y=mean_dist, fill=strain))+
  geom_violin()+
  petite_t200_colors+
  NULL

tapply(summ_overlap_pos$mean_dist, summ_overlap_pos$strain, summary)
# $Petite
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.51   13.38   13.62   13.65   13.92   15.22 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.98   13.83   14.08   14.10   14.33   15.60

t.test(summ_overlap_pos$mean_dist~summ_overlap_pos$strain)
# t = -18.222, df = 970.44, p-value < 2.2e-16
# 95 percent confidence interval:
#   -0.5008196 -0.4034362
# sample estimates:
# mean in group Petite mean in group Grande 
# 13.64515             14.09728 

anova_mean_dist=aov(summ_overlap_pos$mean_dist~summ_overlap_pos$strain)
summary(anova_mean_dist)
#                          Df Sum Sq Mean Sq F value Pr(>F)    
# summ_overlap_pos$strain   1   51.1   51.10     332 <2e-16 ***
# Residuals               998  153.6    0.15     
51.1/(51.1+153.6) #0.2496336
#25% of the variance is explained by the strain variable

ggplot(summ_overlap_pos, aes(x=strain, y=mean_hopkins, fill=strain))+
  geom_violin()+
  petite_t200_colors+
  NULL

tapply(summ_overlap_pos$mean_hopkins, summ_overlap_pos$strain, summary)
# $Petite
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9745  0.9809  0.9822  0.9822  0.9836  0.9889 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9767  0.9816  0.9829  0.9828  0.9843  0.9879

t.test(summ_overlap_pos$mean_hopkins~summ_overlap_pos$strain)
# t = -4.5006, df = 996.43, p-value = 7.574e-06
# 95 percent confidence interval:
#   -0.0008408220 -0.0003302254
# sample estimates:
# mean in group Petite mean in group Grande 
# 0.9822123            0.9827978



#### Paper Figure ####


img_petite <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_overlap_189cells_3overlap_cropped.png")
img_plot_petite <- rasterGrob(img_petite, interpolate = TRUE)

img_grande <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_overlap_201cells_3overlap_cropped.png")
img_plot_grande <- rasterGrob(img_grande, interpolate = TRUE)

# Create text annotations
text_petite <- textGrob("Petite", gp = gpar(fontsize = 14, fontface = "bold"))
text_grande <- textGrob("Grande", gp = gpar(fontsize = 14, fontface = "bold"))

# Create ggplot objects for the images with annotations
p_petite <- ggplot() + 
  annotation_custom(img_plot_petite, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_petite, xmin = 0.2, xmax = 0.2, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_petite

p_grande <- ggplot() + 
  annotation_custom(img_plot_grande, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_grande, xmin = 0.2, xmax = 0.2, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_grande



p_frac_size=ggplot(summ_clust_size, aes(x=strain, y=mean_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Mean Fracture Size\n(Number of Cells)')+
  NULL
p_frac_size

p_overlap=ggplot(mean_summ_overlap_size, 
                 aes(x=cluster_size, y=mean_overlap, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=cluster_size, y=mean_overlap, ymin=ci_lower,
                  ymax=ci_upper, fill=strain), alpha=0.5, linetype='blank')+
  petite_t200_colors+
  xlab('Cluster size (number of cells)')+
  ylab('Mean Overlap')+
  guides(col='none', fill='none')+
  NULL
p_overlap



p_mean_dist_overlaps=ggplot(summ_overlap_pos, aes(x=strain, y=mean_dist, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  labs(x='Strain', y='Mean Distance\nBetween Overlaps')+
  guides(fill='none')+
  NULL
p_mean_dist_overlaps




fig_physics_sim_mean_dist=plot_grid(p_petite, p_grande, p_overlap, p_frac_size, p_mean_dist_overlaps,
                          labels=c('A)', 'B)', 'C)', 'D)', 'E)'), ncol=2, label_size=16, rel_heights=c(1.2,1,1))
fig_physics_sim_mean_dist


ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig3_physics_sims_13nov2024_ancestors.png',
       plot=fig_physics_sim_mean_dist, dpi='retina', width=10, height=12, bg='white')


