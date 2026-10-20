
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


# setwd('~/Desktop/version_supreme_corrected_results/supp_fig4_biophysics_sim/')
# setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/corrections_supreme_version_sep2026/supp_fig4_biophysics_sim/')


# Cluster size difference ####

# clust_size=data.frame()
# 
# for(j in c("petite", "grande")){
#   temp_df=read.csv(paste(j, "_frac_size_30sim_100n_88overlap.csv", sep=""), header=TRUE)
# 
#   if(j=='petite'){
#     temp_df$strain='Petite'
#   } else {
#     temp_df$strain='Grande'
#   }
#   clust_size=rbind(clust_size, temp_df)
# }
# 
# clust_size$strain=factor(clust_size$strain, levels=c('Petite', 'Grande'))
# 
# write.csv(clust_size, file="physics_sim_frag_size_sync_vs_async_88overlap.csv", row.names=FALSE)

clust_size=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_5_size_difference_physics_sim/physics_sim_frag_size_sync_vs_async_88overlap.csv", header = TRUE)

clust_size$strain=factor(clust_size$strain, levels=c('Petite', 'Grande'))


summary(clust_size)

summ_clust_size=clust_size %>%
  group_by(strain, sim_number) %>%
  summarise(mean_size=mean(size_fracture),
            mean_volume=mean(volume))

summ_clust_size %>% group_by(strain) %>% summarise(mean_cluster_size=mean(mean_size))
# target size is 339
# 331.9747 cells for 85 overlap
# 336.7113 cells for 87 overlap
# 338.8387 cells for 88 overlap <---
# 345.4837 cells for 90 overlap
# 356.2283 cells for 95 overlap

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




summ_clust_size %>%
  group_by(strain) %>%
  summarise(median_fracture_size=median(mean_size),
            median_cluster_volume=median(mean_volume))
#overlap threshold 85
#   strain median_fracture_size   mean_cluster_volume
# 1 Petite               312.765              28976.86
# 2 Grande               339.385              35983.85
# 27 (26.6) cell difference between grande and petite

# Strain accumulation ####

petite_overlap=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_5_size_difference_physics_sim/petite_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
petite_overlap$strain="Petite"
grande_overlap=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_5_size_difference_physics_sim/grande_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
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
  geom_hline(yintercept=85, linetype="dashed", color = "black")+
  petite_t200_colors+
  theme_bw()+
  xlab('Cluster size (number of cells)')+
  ylab('Mean Overlap')+
  NULL

mean_summ_overlap_size[mean_summ_overlap_size$cluster_size==340,]


# Overlap Positions ####

# process_overlap_pos_physics_sim_1nov2024.py -i . -f _overlap_pos_30sim_100n_88overlap.csv -o stats_overlap_pos_88overlap.csv

overlap_pos_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_5_size_difference_physics_sim/stats_overlap_pos_88overlap.csv", header=TRUE)

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
# 12.51   13.00   13.26   13.25   13.48   14.26 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.18   13.71   13.96   13.92   14.17   14.65 

t.test(summ_overlap_pos$mean_dist~summ_overlap_pos$strain)
# t = -14.179, df = 197.31, p-value < 2.2e-16







#### Paper Figure ####


img_petite <- readPNG("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/matlab_plots/petite_overlap_302cells_cropped.png")
img_plot_petite <- rasterGrob(img_petite, interpolate = TRUE)

img_grande <- readPNG("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/matlab_plots/grande_overlap_329cells_cropped.png")
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



wilcox.test(summ_clust_size$mean_size~summ_clust_size$strain)
# p-value < 2.2e-16


p_frac_size=ggplot(summ_clust_size, aes(x=strain, y=mean_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  scale_fill_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  guides(fill='none')+
  labs(x='Strain', y='Mean Fracture Size\n(Number of Cells)')+
  NULL
p_frac_size

wilcox.test(summ_clust_size$mean_size~summ_clust_size$strain)
# p-value < 2.2e-16


p_volume=ggplot(summ_clust_size, aes(x=strain, y=mean_volume, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  scale_fill_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  guides(fill='none')+
  labs(x='Strain', y = expression(atop("Mean Cluster", "Volume"~(plain(µ)*m^3))))+
  NULL
p_volume

wilcox.test(summ_clust_size$mean_volume~summ_clust_size$strain)
# p-value < 2.2e-16


p_overlap=ggplot(mean_summ_overlap_size, 
                 aes(x=cluster_size, y=mean_overlap, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=cluster_size, y=mean_overlap, ymin=ci_lower,
                  ymax=ci_upper, fill=strain), alpha=0.5, linetype='blank')+
  scale_fill_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  scale_color_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  xlab('Cluster size (number of cells)')+
  ylab('Mean Overlap')+
  guides(col='none', fill='none')+
  NULL
p_overlap



p_mean_dist_overlaps=ggplot(summ_overlap_pos, aes(x=strain, y=mean_dist, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  scale_fill_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  labs(x='Strain', y='Mean Distance\nBetween Overlaps')+
  guides(fill='none')+
  NULL
p_mean_dist_overlaps




fig_physics_sim_mean_dist=plot_grid(p_petite, p_grande, p_overlap, p_frac_size, p_volume, p_mean_dist_overlaps,
                          labels=LETTERS[1:6], ncol=2, label_size=16, rel_heights=c(1.2,1,1))
fig_physics_sim_mean_dist


ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_5_physics_sims_14oct2026.png',
       plot=fig_physics_sim_mean_dist, dpi='retina', width=10, height=12, bg='white')


