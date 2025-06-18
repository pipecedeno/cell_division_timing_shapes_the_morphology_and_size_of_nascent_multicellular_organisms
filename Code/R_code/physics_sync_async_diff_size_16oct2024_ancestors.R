
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

setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/supp_fig3_physics_sim')


# Cluster size difference ####

# clust_size=data.frame()
# 
# for(j in c("petite", "grande")){
#   temp_df=read.csv(paste(j, "_frac_size_30sim_500n_1.2aspr_10attempts_70overlap_v2.csv", sep=""), header=TRUE)
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
# clust_size$gyration_diam=clust_size$gyration_diam*10**6
# 
# write.csv(clust_size, file="physics_sim_frag_size_sync_vs_async_10june2025.csv", row.names=FALSE)

clust_size=read.csv("physics_sim_frag_size_sync_vs_async_10june2025.csv", header = TRUE)
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
# 1 Petite               312.              39294.          30.7
# 2 Grande               335.              46741.          32.0
# 23 cell difference between grande and petite

# Strain accumulation ####

# overlap_size_df=data.frame()
# 
# for(j in c("petite", "grande")){
#   for(i in seq(1,5)){
#     temp_df=read.csv(paste("24oct2024_overlap_acum_100net_20n/", j, "_overlap_clust_size_1.2ar_10attempts_", i, ".csv", sep=""), header=TRUE)
#     if(j=='petite'){
#       temp_df$strain='Petite'
#     } else {
#       temp_df$strain='Grande'
#     }
#     overlap_size_df=rbind(overlap_size_df, temp_df)
#   }
# }


petite_overlap=read.csv("overlap_accumulation_2dec2024/petite_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
petite_overlap$strain="Petite"
grande_overlap=read.csv("overlap_accumulation_2dec2024/grande_overlap_clust_size_1.2ar_10attempts.csv", header=TRUE)
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

# process_overlap_pos_physics_sim_1nov2024.py -i . -f _overlap_pos_30sim_500n_1.2aspr_10attempts_70overlap.csv -o stats_overlap_pos_10jun2025.csv

overlap_pos_df=read.csv("stats_overlap_pos_10jun2025.csv", header=TRUE)

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
# 13.92   14.84   15.12   15.13   15.42   16.85 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.56   15.38   15.62   15.62   15.88   16.71

t.test(summ_overlap_pos$mean_dist~summ_overlap_pos$strain)
# t = -20.178, df = 971.8, p-value < 2.2e-16
# 95 percent confidence interval:
#   -0.5417887 -0.4457469


anova_mean_dist=aov(summ_overlap_pos$mean_dist~summ_overlap_pos$strain)
summary(anova_mean_dist)
#                          Df Sum Sq Mean Sq F value Pr(>F)    
# summ_overlap_pos$strain   1  60.95   60.95   407.2 <2e-16 ***
# Residuals               998 149.40    0.15 
60.95/(60.95+149.40) #0.2897552
#28% of the variance is explained by the strain variable





#### Paper Figure ####


img_petite <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/matlab_plots/petite_overlap_302cells_cropped.png")
img_plot_petite <- rasterGrob(img_petite, interpolate = TRUE)

img_grande <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/matlab_plots/grande_overlap_329cells_cropped.png")
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
  scale_fill_manual(values=c("#F0D77BFF", "#AE93BEFF"))+
  guides(fill='none')+
  labs(x='Strain', y='Mean Fracture Size\n(Number of Cells)')+
  NULL
p_frac_size

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




fig_physics_sim_mean_dist=plot_grid(p_petite, p_grande, p_overlap, p_frac_size, p_mean_dist_overlaps,
                          labels=c('A', 'B', 'C', 'D', 'E'), ncol=2, label_size=16, rel_heights=c(1.2,1,1))
fig_physics_sim_mean_dist


ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig3_physics_sims_10juune2025_ancestors.png',
       plot=fig_physics_sim_mean_dist, dpi='retina', width=10, height=12, bg='white')


