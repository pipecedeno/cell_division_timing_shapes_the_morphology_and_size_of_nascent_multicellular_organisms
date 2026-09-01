

#Date: 16oct2024
# This code is used to analyze the results using random edge selection for the edge degree fragmentation
# updated to use Petite, Petite w/o Delay and Grande

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

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))


# setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/supp_fig5_frag_clusters_random/')

# Fragmentation summary ####

# frag_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_15e_100g_random/fragmentation_inf.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   frag_df=rbind(frag_df, temp_df)
# }
# frag_df$strain <- ifelse(frag_df$strain=='petite', 'Petite', frag_df$strain)
# frag_df$strain <- ifelse(frag_df$strain=='petite-second-only', 'Petite w/o Delay', frag_df$strain)
# frag_df$strain <- ifelse(frag_df$strain=='grande', 'Grande', frag_df$strain)
# frag_df$strain=factor(frag_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# summary(frag_df)
# 
# table(frag_df$generation)
# 
# write.csv(frag_df, file="growth_frag_random_frag_df_25june2025.csv", row.names = FALSE)


frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_6_random_fragmentation/growth_frag_random_frag_df_25june2025.csv", 
                 header=TRUE)

frag_df$strain=factor(frag_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
summary(frag_df)

table(frag_df$generation)


#### Size at fracture ####

ggplot(frag_df[frag_df$generation%%20==0,],
       aes(x=as.factor(generation), y=cluster_size, fill=strain))+
  geom_violin()+
  ylab("Cluster size at fracture")+xlab("Generations")+
  petite_t200_colors+
  theme(legend.position="bottom")+
  NULL

ggplot(frag_df, aes(x=strain, y=cluster_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('Strain')+ylab('Mean Fracture Size')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

frag_df %>% 
  group_by(strain) %>%
  summarise(mean=mean(cluster_size))
#   strain            mean
# 1 Petite            217.8538
# 2 Petite w/o delay  354.8998
# 3 Grande            340.3699

pairwise.t.test(frag_df$cluster_size, frag_df$strain, method="bonferroni")
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16          
# P value adjustment method: holm 

mean_clust_size=frag_df %>%
  group_by(strain, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))


ggplot(mean_clust_size, aes(x=generation, y=mean_size, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=strain), 
              alpha=0.3, linetype='blank')+
  # ggtitle("Mean Fracture Size")+
  ylab("Mean fracure size")+xlab('Generation')+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL



#### Fracture proportion ####

summary(frag_df)

mean_frag_prop=frag_df %>%
  group_by(strain, generation) %>%
  summarise(mean_frac=mean(proportion_propagule),
            sd_frac=sd(proportion_propagule),
            min_frac=min(proportion_propagule),
            max_frac=max(proportion_propagule))

ggplot(mean_frag_prop)+
  geom_line(aes(x=generation, y=mean_frac, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_frac, ymin=mean_frac-sd_frac, ymax=mean_frac+sd_frac, fill=strain), 
              alpha=0.2, linetype='blank')+
  # ggtitle('Mean Propagule Proportion')+
  ylab("Mean propagule proportion")+xlab('Generation')+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL


ggplot(frag_df[frag_df$generation%%20==0,], 
       aes(x=as.factor(generation), y=proportion_propagule, fill=strain))+
  geom_violin()+
  xlab("Generation")+ylab("Proportion propagule size")+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL

ggplot(frag_df, aes(x=strain, y=proportion_propagule, fill=strain))+
  geom_violin()+
  xlab('Generation')+
  ylab('Proportion Propagule Size')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

frag_df %>%
  group_by(strain) %>%
  summarise(mean=mean(proportion_propagule))
#  strain            mean
# 1 Petite           0.398
# 2 Petite w/o delay 0.397
# 3 Grande           0.455

pairwise.t.test(frag_df$proportion_propagule, frag_df$strain)
#                   Petite  Petite w/o delay
# Petite w/o delay 6.2e-05 -               
# Grande           < 2e-16 < 2e-16 


# Network Diameter ####

# diameter_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_15e_100g_random/networks_diameter.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   diameter_df=rbind(diameter_df, temp_df)
# }
# diameter_df$strain <- ifelse(diameter_df$strain=='petite', 'Petite', diameter_df$strain)
# diameter_df$strain <- ifelse(diameter_df$strain=='petite-second-only', 'Petite w/o Delay', diameter_df$strain)
# diameter_df$strain <- ifelse(diameter_df$strain=='grande', 'Grande', diameter_df$strain)
# diameter_df$strain=factor(diameter_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# summary(diameter_df)
# summary(frag_df)
# 
# write.csv(diameter_df, file='growth_frag_random_diameter_df_26june2025.csv', row.names = FALSE)


diameter_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_6_random_fragmentation/growth_frag_random_diameter_df_26june2025.csv", header=TRUE)

diameter_df$strain=factor(diameter_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
summary(diameter_df)
summary(frag_df)

#To normalize the network diameter metric I need to first join the diameter and size at fracture information
#to later divide network diameter by cluster size

diam_size <- left_join(diameter_df, frag_df, 
                       by = c("sim_number", "generation", "strain"))

summary(diam_size)

# diam_size$norm_diameter=diam_size$diameter/diam_size$cluster_size
diam_size$norm_diameter=diam_size$diameter/(6.64385619*log10(diam_size$cluster_size)-1)
summary(diam_size)


ggplot(diam_size[diam_size$generation%%20==0,], 
       aes(x=as.factor(generation), y=norm_diameter, fill=strain))+
  geom_violin()+
  xlab("Generation")+ylab("Normalized Diameter")+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL


#normalizing diameter
ggplot(diam_size, aes(x=strain, y=norm_diameter, fill=strain))+
  geom_violin()+
  xlab('Strain')+ylab('Normalized Diameter')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

diam_size %>%
  group_by(strain) %>%
  summarise(mean=mean(norm_diameter, na.rm=TRUE))
#  strain            mean
# 1 Petite           0.912
# 2 Petite w/o delay 0.997
# 3 Grande           0.990

pairwise.t.test(diam_size$norm_diameter, diam_size$strain)
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16 


# Plots of diameter by generation
summ_diam=diam_size %>%
  group_by(strain, generation) %>%
  summarise(mean_diam=mean(norm_diameter),
            sd_diam=sd(norm_diameter),
            median_diam=median(norm_diameter))

mean_summ_diam=summ_diam %>%
  group_by(strain, generation) %>%
  summarise(across(mean_diam, list(mean = mean, sd = sd), .names = "{.fn}_diam"))

ggplot(summ_diam)+
  geom_line(aes(x=generation, y=mean_diam, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_diam, ymin=mean_diam-sd_diam, ymax=mean_diam+sd_diam, fill=strain), alpha=0.2)+
  # ggtitle("Mean of Normalized Network Diameter")+
  xlab("Generation")+
  ylab("Normalized Diameter")+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL



# Paper Figures ####


supp_p1_gen=ggplot(mean_clust_size, aes(x=generation, y=mean_size, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=strain), 
              alpha=0.3, linetype='blank')+
  ylab("Mean fracure\nsize")+xlab('Generation')+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p1_gen

supp_p2_gen=ggplot(mean_frag_prop)+
  geom_line(aes(x=generation, y=mean_frac, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_frac, ymin=mean_frac-sd_frac, ymax=mean_frac+sd_frac, fill=strain), 
              alpha=0.2, linetype='blank')+
  ylab("Mean propagule \nproportion")+xlab('Generation')+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p2_gen

supp_p3_gen=ggplot(summ_diam)+
  geom_line(aes(x=generation, y=mean_diam, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_diam, ymin=mean_diam-sd_diam, ymax=mean_diam+sd_diam, fill=strain), alpha=0.2)+
  xlab("Generation")+
  ylab("Mean Normalized \nDiameter")+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p3_gen


supp_p1_dist=ggplot(frag_df[frag_df$generation%%20==0,],
                    aes(x=as.factor(generation), y=cluster_size, fill=strain))+
  geom_violin()+
  ylab("Cluster size at \nfracture")+xlab("Generation")+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p1_dist

supp_p2_dist=ggplot(frag_df[frag_df$generation%%20==0,], 
                    aes(x=as.factor(generation), y=proportion_propagule, fill=strain))+
  geom_violin()+
  xlab("Generation")+ylab("Proportion\npropagule size")+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p2_dist

supp_p3_dist=ggplot(diam_size[diam_size$generation%%20==0,], 
                    aes(x=as.factor(generation), y =norm_diameter, fill=strain))+
  geom_violin()+
  xlab("Generation")+ylab("Normalized\nDiameter")+
  petite_t200_colors+
  # theme_classic(base_size = 10)+
  guides(col='none', fill='none')+
  NULL
supp_p3_dist


legend_labels <- unique(mean_clust_size$strain)
legend_colors <- mycolors
names(legend_colors) <- legend_labels

common_legend <- cowplot::get_legend(
  ggplot(diam_size, aes(x = as.factor(generation), y = cluster_size, fill = strain)) +
    geom_violin() +
    scale_fill_manual(values = legend_colors, labels = c("Petite", "Petite w/o Delay", "Grande")) +
    guides(fill = guide_legend(nrow = 1, title = NULL))
)


figure_supp_v2=plot_grid(supp_p1_dist,supp_p1_gen, supp_p2_dist,supp_p2_gen, supp_p3_dist,supp_p3_gen, 
                         labels=c('A', 'B', 'C', 'D', 'E', 'F'), ncol=2, align='hv', label_size=16)
figure_supp_v2

figure_supp_v2_with_legend <- plot_grid(figure_supp_v2, common_legend, ncol = 1, rel_heights = c(1, 0.05), label_size=11)
figure_supp_v2_with_legend

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig_6_network_trajectories_random_28aug2026.png',
       plot=figure_supp_v2_with_legend, dpi='retina', width=10, height=11, bg = 'white')




