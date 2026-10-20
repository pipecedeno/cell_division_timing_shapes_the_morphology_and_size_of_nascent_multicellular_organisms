


# Date: 12 sep 2025

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
library(see)
library(gridExtra)
library(metR)

# setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/test_death_prob/more_parameters_alive_status')
# setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/corrections_supreme_version_sep2026/supp_fig10_11_death_prob')

theme_set(theme_classic(base_size = 16))

# death prob variable, constant delay ####

default_mean=105

death_prob_name=sprintf("%02d", seq(0, 30, 5))
death_prob_value=as.character(seq(0, 0.3, 0.05))

delay_name=as.character(seq(-52.5, 52.5, length.out=21))
delay_percentage=as.numeric(delay_name)/default_mean*100

# only loading the sizes of interest
network_sizes=c(seq(250, 1250, 250), 1300)
# network_sizes=c(seq(200, 1200, 200), 1300)

# net_sim_df=data.frame()
# 
# for (j_cont in seq(length(delay_name))){
#   for (i_cont in seq(length(death_prob_value))){
#     temp_delay=delay_name[j_cont]
#     temp_percentages=delay_percentage[j_cont]
#     temp_death_prob=death_prob_value[i_cont]
#     temp_death_prob_name=death_prob_name[i_cont]
#     temp_df=read.csv(paste("growth_no_frag_results/prob_", temp_death_prob_name, "/test_15_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
# 
#     temp_df=temp_df[temp_df$num_nodes %in% network_sizes,]
# 
#     temp_df$variation=15
#     temp_df$diff_mean=temp_delay
#     temp_df$percentage_diff=temp_percentages
#     temp_df$death_percentage=temp_death_prob
#     net_sim_df=rbind(net_sim_df, temp_df)
#   }
# }
# summary(net_sim_df)
# 
# table(net_sim_df$num_nodes)


# write.csv(net_sim_df, file="growth_no_frag_death_prob_16sep2026.csv", row.names = FALSE)


# net_sim_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_11-14_cell_death_simulations/growth_no_frag_death_prob_16sep2026.csv", header=TRUE)

net_sim_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_11-14_cell_death_simulations/growth_no_frag_death_prob_0SD.csv", header=TRUE)

summary(net_sim_df)


#### Network Diameter ####

net_sim_df$norm_diameter=net_sim_df$diameter/(6.64385619*log10(net_sim_df$num_nodes)-1)

summ_diameter <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(norm_diameter),
    sd_diameter = sd(norm_diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

summary(summ_diameter[summ_diameter$num_nodes==1300,]$mean_diameter)


# Effect combination between death and delay

ggplot(summ_diameter[summ_diameter$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(breaks=seq(0.8, 1.7, 0.1)) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL



# Effect for synchronous strain

ggplot(summ_diameter[summ_diameter$num_nodes %in% network_sizes & summ_diameter$percentage_diff==0,],
       aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Network Diameter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL



#### Undivided cells ####

summ_undivided <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_undivided = mean(cases_mother_with_undivided_cells),
    sd_undivided = sd(cases_mother_with_undivided_cells),
    n = n(),
    se_undivided = sd_undivided / sqrt(n),
    ci_lower = mean_undivided - qt(0.975, n-1) * se_undivided,
    ci_upper = mean_undivided + qt(0.975, n-1) * se_undivided,
    .groups = 'drop'
  )

summary(summ_undivided$mean_undivided)


# Interaction of cell death and delay

ggplot(summ_undivided[summ_undivided$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_undivided)) +
  geom_contour_filled(breaks = seq(0, 180, 20)) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nDaughter\nCells")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL

# Effect for synchronous strain

ggplot(summ_undivided[summ_undivided$num_nodes %in% network_sizes & summ_undivided$percentage_diff==0,],
       aes(x=death_percentage, y=mean_undivided, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Mothers with\nDelayed Daughter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL



#### Filamentous branches ####


summ_filamentous <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_filament = mean(num_filamentous_branches),
    sd_filament = sd(num_filamentous_branches),
    n = n(),
    se_filament = sd_filament / sqrt(n),
    ci_lower = mean_filament - qt(0.975, n-1) * se_filament,
    ci_upper = mean_filament + qt(0.975, n-1) * se_filament,
    .groups = 'drop'
  )


summary(summ_filamentous$mean_filament)


# Interaction of cell death and delay

ggplot(summ_filamentous[summ_filamentous$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_filament)) +
  geom_contour_filled(breaks = seq(0, 180, 20)) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


# Effect for synchronous strain

ggplot(summ_filamentous[summ_filamentous$num_nodes %in% network_sizes & summ_filamentous$percentage_diff==0,],
       aes(x=death_percentage, y=mean_filament, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Filament Branches')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL



#### Motif difference ####

summ_motif_diff <- summ_filamentous %>%
  left_join(summ_undivided, 
            by = c("percentage_diff", "death_percentage", "num_nodes"),
            suffix = c("_filamentous", "_undivided")) %>%
  mutate(mean_motif_difference=mean_filament - mean_undivided,
         mean_motif_difference_norm = (mean_filament - mean_undivided)/num_nodes)

summary(summ_motif_diff[summ_motif_diff$num_nodes>=10,]$mean_motif_difference_norm)


# Interaction of cell death and delay

ggplot(summ_motif_diff[summ_motif_diff$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_motif_difference_norm)) +
  geom_contour_filled(breaks=seq(-0.10, 0.18, 0.04)) +
  scale_fill_brewer(palette = "Purples", name = "Mean\nMotif\nDifference\nNorm")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


ggplot(summ_motif_diff[summ_motif_diff$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_motif_difference)) +
  geom_contour_filled() + #breaks=seq(-0.10, 0.18, 0.04)
  scale_fill_brewer(palette = "Purples", name = "Mean\nMotif\nDifference\nNorm")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL

# Effect for synchronous strain

ggplot(summ_motif_diff[summ_motif_diff$num_nodes %in% network_sizes & summ_motif_diff$percentage_diff==0,],
       aes(x=death_percentage, y=mean_motif_difference_norm, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Normalized Mean Motif Difference')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL




#### Max edge degree ####

summ_max_edge <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_max_edge = mean(max_edge_degree),
    sd_max_edge = sd(max_edge_degree),
    n = n(),
    se_max_edge = sd_max_edge / sqrt(n),
    ci_lower = mean_max_edge - qt(0.975, n-1) * se_max_edge,
    ci_upper = mean_max_edge + qt(0.975, n-1) * se_max_edge,
    .groups = 'drop'
  )

summary(summ_max_edge[summ_max_edge$num_nodes>=10,]$mean_max_edge)


# Interaction of cell death and delay

ggplot(summ_max_edge[summ_max_edge$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_max_edge)) +
  geom_contour_filled(breaks = seq(3, 24, 3)) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


# Effect for synchronous strain

ggplot(summ_max_edge[summ_max_edge$num_nodes %in% network_sizes & summ_max_edge$percentage_diff==0,],
       aes(x=death_percentage, y=mean_max_edge, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Max Edge Degree')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL


#### Normalized max edge degree ####

net_sim_df$norm_max_edge_degree=net_sim_df$max_edge_degree/(6.64385619*log10(net_sim_df$num_nodes)-1)

summ_norm_max_edge <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_norm_max_edge = mean(norm_max_edge_degree),
    sd_norm_max_edge = sd(norm_max_edge_degree),
    n = n(),
    se_mean_max_edge = sd_norm_max_edge / sqrt(n),
    ci_lower = mean_norm_max_edge - qt(0.975, n-1) * se_mean_max_edge,
    ci_upper = mean_norm_max_edge + qt(0.975, n-1) * se_mean_max_edge,
    .groups = 'drop'
  )

summary(summ_norm_max_edge[summ_norm_max_edge$num_nodes>=10,]$mean_norm_max_edge)


# Interaction of cell death and delay

ggplot(summ_norm_max_edge[summ_norm_max_edge$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_norm_max_edge)) +
  geom_contour_filled(breaks=seq(0.6, 1.2, 0.075)) +
  scale_fill_brewer(palette = "Purples", name = "Normalized\nMax Edge\nDegree")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


# Effect for synchronous strain

ggplot(summ_norm_max_edge[summ_norm_max_edge$num_nodes %in% network_sizes & summ_norm_max_edge$percentage_diff==0,],
       aes(x=death_percentage, y=mean_norm_max_edge, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Normalized Max Edge Degree')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  NULL


# Supplementary Figure 12 ####

sizes_sf10=seq(250, 1250, 250)
# sizes_sf10=seq(200, 1200, 200)

# using not normalized
summ_diameter_not_norm <- net_sim_df %>%
  group_by(percentage_diff, death_percentage, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(diameter),
    sd_diameter = sd(diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )


p_diam=ggplot(summ_diameter[summ_diameter$num_nodes %in% sizes_sf10 & summ_diameter$percentage_diff==0,],
              aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Network Diameter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_diam

p_diam_not_norm=ggplot(summ_diameter_not_norm[summ_diameter_not_norm$num_nodes %in% sizes_sf10 & summ_diameter_not_norm$percentage_diff==0,],
              aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes), fill=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  # geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper))+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, color=NA)+
  theme_classic()+
  labs(x='Death Probability', y='Mean Network Diameter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  scale_fill_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_diam_not_norm

p_undivided=ggplot(summ_undivided[summ_undivided$num_nodes %in% sizes_sf10 & summ_undivided$percentage_diff==0,],
                   aes(x=death_percentage, y=mean_undivided, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Mothers with\nDelayed Daughter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_undivided

p_filament=ggplot(summ_filamentous[summ_filamentous$num_nodes %in% sizes_sf10 & summ_filamentous$percentage_diff==0,],
                  aes(x=death_percentage, y=mean_filament, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Filament Branches')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_filament

p_max_edge=ggplot(summ_max_edge[summ_max_edge$num_nodes %in% sizes_sf10 & summ_max_edge$percentage_diff==0,],
                  aes(x=death_percentage, y=mean_max_edge, col=as.factor(num_nodes), fill=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  # geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper))+
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, color=NA)+
  theme_classic()+
  labs(x='Death Probability', y='Mean Max Edge Degree')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  scale_fill_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_max_edge


pairwise.wilcox.test(net_sim_df[net_sim_df$num_nodes==1250 & net_sim_df$percentage_diff==0,]$max_edge_degree, 
                     net_sim_df[net_sim_df$num_nodes==1250 & net_sim_df$percentage_diff==0,]$death_percentage,
                     p.adjust.method = 'bonferroni')
# Results for 15 SD
#      0       0.05    0.1     0.15    0.2     0.25   
# 0.05 0.00028 -       -       -       -       -      
# 0.1  1.00000 0.01351 -       -       -       -      
# 0.15 1.00000 0.02198 1.00000 -       -       -      
# 0.2  0.00223 4.6e-06 0.30696 1.00000 -       -      
# 0.25 1.00000 0.01197 1.00000 1.00000 1.00000 -      
# 0.3  1.8e-05 9.7e-08 0.00585 0.06133 1.00000 0.27831
#significant differences
# 0.3 to 0, 0.05, 0.01
# 0.25 to 0.05
# 0.2 to 0, 0.05
# 0.05 is different from everything else

# Results for 0 SD (exponentially growing network)
#      0       0.05    0.1     0.15    0.2    0.25  
# 0.05 4.1e-06 -       -       -       -      -     
# 0.1  0.0020  1.0000  -       -       -      -     
# 0.15 0.0209  1.0000  1.0000  -       -      -     
# 0.2  1.0000  0.7911  0.4430  1.0000  -      -     
# 0.25 3.9e-06 3.2e-07 3.4e-06 5.5e-05 0.1138 -     
# 0.3  1.9e-10 7.0e-09 1.2e-07 1.6e-06 0.0031 1.0000

# temporary plot to extract the legend from
temp_plot <- ggplot(summ_diameter[summ_diameter$num_nodes %in% sizes_sf10 & summ_diameter$percentage_diff==0,],
                    aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  scale_color_brewer(palette="Purples", name="Network Size")+
  theme(legend.position = "bottom",           # Position at bottom
        legend.direction = "horizontal")      # Make it horizontal

# Extract the legend
shared_legend <- ggpubr::get_legend(temp_plot)


# Create the plot grid without legend
plots_grid <- plot_grid(p_diam_not_norm, p_max_edge, ncol = 2)

# Combine the plots with the shared legend
final_plot <- plot_grid(plots_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig12_cell_death_19sep2026_0SD.svg',
       plot=final_plot, dpi='retina', height=4, width=8, bg='white')


# Fragmentation results ####


#### Fragmentation size ####

# frag_df=data.frame()
# 
# for (j_cont in seq(length(delay_name))){
#   for (i_cont in seq(length(death_prob_value))){
#     temp_delay=delay_name[j_cont]
#     temp_percentages=delay_percentage[j_cont]
#     temp_death_prob=death_prob_value[i_cont]
#     temp_death_prob_name=death_prob_name[i_cont]
#     temp_df=read.csv(paste("fragmentation_simulations/prob_", temp_death_prob_name, "/test_15_var_", temp_delay, "_diff/fragmentation_inf.csv", sep=""), header=TRUE)
#     
#     # temp_df=temp_df[temp_df$num_nodes %in% network_sizes,]
#     
#     temp_df$variation=15
#     temp_df$diff_mean=temp_delay
#     temp_df$percentage_diff=temp_percentages
#     temp_df$death_percentage=temp_death_prob
#     frag_df=rbind(frag_df, temp_df)
#   }
# }
# summary(frag_df)
# 
# write.csv(frag_df, "fragmentation_data_death_prob_16sep2026.csv", row.names = FALSE)

frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_11-14_cell_death_simulations/fragmentation_data_death_prob_16sep2026.csv", header=TRUE)
summary(frag_df)

mean_clust_size=frag_df %>%
  group_by(death_percentage, percentage_diff, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size),
            n = n(),
            se = sd_size / sqrt(n),
            lower_ci = mean_size - qt(0.975, df = n - 1) * se,
            upper_ci = mean_size + qt(0.975, df = n - 1) * se)

mean_mean_clust_size <- mean_clust_size %>%
  group_by(death_percentage, percentage_diff) %>%
  summarize(mean = mean(mean_size),
            sd = sd(mean_size),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

summary(mean_mean_clust_size)

ggplot(mean_mean_clust_size, aes(x=death_percentage, y=mean, col=percentage_diff, group=percentage_diff)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  xlab('Death Probability') +
  ylab('Mean fracture size') +
  # color_diff_mean+
  NULL


ggplot(mean_clust_size[mean_clust_size$percentage_diff %in% c(-50, -25, 0, 25, 50),] %>%
         mutate(percentage_diff=factor(as.character(percentage_diff), levels=c("50", "25", "0", "-25", "-50"))), 
       aes(x=generation, y=mean_size, fill=percentage_diff, col=percentage_diff)) +
  # geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.4, color=NA)+
  facet_grid(percentage_diff~death_percentage)+
  stat_summary(fun="median", geom="crossbar")+
  xlab('Death Probability') +
  ylab('Mean fracture size') +
  guides(fill='none')+
  NULL



ggplot(mean_mean_clust_size, aes(x = death_percentage, y = percentage_diff, z = mean)) +
  geom_contour_filled(breaks = seq(0, 2400, 300)) +
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


ggplot(mean_mean_clust_size, aes(x = death_percentage, y = percentage_diff, z = mean)) +
  geom_contour_filled(breaks = c(2**seq(3, 12, 1)))+
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


#### Propagule Proportion ####

mean_propagule_prop=frag_df %>%
  group_by(death_percentage, percentage_diff, generation) %>%
  summarise(mean_propagule=mean(proportion_propagule))

mean_mean_propagule_prop <- mean_propagule_prop %>%
  group_by(death_percentage, percentage_diff) %>%
  summarize(mean = mean(mean_propagule),
            sd = sd(mean_propagule),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

summary(mean_mean_propagule_prop)

ggplot(mean_mean_propagule_prop, aes(x=death_percentage, y=mean, col=percentage_diff, group=percentage_diff)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  xlab('Death Probability') +
  ylab('Propagule Proportion') +
  # color_diff_mean+
  NULL



ggplot(mean_mean_propagule_prop, aes(x = death_percentage, y = percentage_diff, z = mean)) +
  geom_contour_filled() +
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Propagule\nProportion")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL




#### Alive cells in fragmentation ####

table(frag_df$alive_cells)

summ_alive_cells=frag_df %>%
  group_by(death_percentage, percentage_diff) %>%
  summarise(perc_none_alive = sum(alive_cells=='none')/n()*100,
            perc_one_alive = sum(alive_cells=='one')/n()*100,
            perc_both_alive = sum(alive_cells=='both')/n()*100) %>%
  pivot_longer(
    cols = c(perc_none_alive, perc_one_alive, perc_both_alive),
    names_to = "alive_group",
    values_to = "percentage"
  ) %>%
  mutate(alive_group = recode(alive_group,
                              perc_none_alive = "none",
                              perc_one_alive  = "one",
                              perc_both_alive = "both")) %>%
  mutate(alive_group = factor(alive_group, levels=c('both', 'one', 'none')))

# Probability of 15% of dying
#none alive
0.15*0.15*100 # 2.25

#one alive
0.85*0.15*100*2 # 25.5

#both alive
0.85*0.85*100 # 72.25

null_expectation_binomial <- data.frame(
  death_percentage = unique(summ_alive_cells$death_percentage)
) %>%
  mutate(
    none = death_percentage^2 * 100,
    one  = 2 * death_percentage * (1 - death_percentage) * 100,
    both = (1 - death_percentage)^2 * 100
  ) %>%
  pivot_longer(cols = c(none, one, both),
               names_to = "alive_group",
               values_to = "expected_percentage")


simulate_geometric_null <- function(death_prob, target = 15, n_sims = 200000) {
  if (death_prob == 0) {
    return(c(both = 100, one = 0, none = 0))
  }
  
  N_A <- rgeom(n_sims, death_prob) + 1
  N_B <- rgeom(n_sims, death_prob) + 1
  
  outcome <- vapply(seq_len(n_sims), function(i) {
    nA <- N_A[i]; nB <- N_B[i]
    seq_labels <- sample(c(rep("A", nA), rep("B", nB)))
    if (length(seq_labels) < target) return(NA_character_)
    
    at_target <- seq_labels[1:target]
    a_done <- sum(at_target == "A")
    b_done <- sum(at_target == "B")
    a_alive <- a_done < nA
    b_alive <- b_done < nB
    if (a_alive & b_alive) "both" else if (a_alive | b_alive) "one" else "none"
  }, character(1))
  
  valid <- outcome[!is.na(outcome)]
  c(both = mean(valid == "both"), one = mean(valid == "one"), none = mean(valid == "none")) * 100
}

simulate_geometric_null(0.1)

null_expectation_geometric <- lapply(round(seq(0, 0.3, 0.05), 2), function(d) {
  vals <- simulate_geometric_null(d)
  tibble(
    death_percentage = d,
    alive_group = names(vals),
    expected_percentage = vals  # already a percentage, no need to multiply again
  )
}) %>% bind_rows()



null_expectation_random_partner <- tibble(death_percentage = round(seq(0, 0.3, 0.05), 2)) %>%
  mutate(q    = 1 - 2*death_percentage,
         both = (1 - death_percentage)*q*100,
         one  = ((1 - death_percentage)*(1 - q) + death_percentage*q)*100,
         none = death_percentage*(1 - q)*100) %>%
  select(-q) %>%
  pivot_longer(c(both, one, none), names_to = "alive_group", values_to = "expected_percentage")




ggplot(summ_alive_cells[summ_alive_cells$death_percentage %in% c(0, 0.15, 0.3),], 
       aes(x=percentage_diff, y=percentage, col=alive_group))+
  geom_line()+
  # geom_hline(data = null_expectation_binomial[null_expectation_binomial$death_percentage %in% c(0, 0.15, 0.3),],
  #            aes(yintercept = expected_percentage, col = alive_group),
  #            linetype = "dashed", linewidth = 0.6) +
  # geom_hline(data = null_expectation_geometric[null_expectation_geometric$death_percentage %in% c(0, 0.15, 0.3),],
  #            aes(yintercept = expected_percentage, col = alive_group),
  #            linetype = "dashed", linewidth = 0.6) +
  geom_hline(data = null_expectation_random_partner[null_expectation_random_partner$death_percentage %in% c(0, 0.15, 0.3),],
             aes(yintercept = expected_percentage, col = alive_group),
             linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~death_percentage)+
  theme_classic()+
  scale_color_brewer(palette = "Purples", name = "Nodes\nAlive at\nFragmentation?")+
  labs(x="Delay %", y="Percentage alive nodes in fragmentation edge")+
  NULL

summ_alive_cells_no_delay = summ_alive_cells %>%
  group_by(death_percentage, alive_group) %>%
  summarise(mean_percentage=mean(percentage),
            min_percentage=min(percentage),
            max_percentage=max(percentage))

ggplot(summ_alive_cells_no_delay, aes(x=death_percentage, y=mean_percentage, col=alive_group, fill=alive_group))+
  geom_line()+
  geom_ribbon(aes(ymin=min_percentage, ymax=max_percentage), alpha=0.5, linetype = 0)+
  theme_classic()+
  # geom_line(data = null_expectation_binomial,
  #            aes(x=death_percentage, y = expected_percentage, col = alive_group),
  #            linetype = "dashed", linewidth = 0.6) +
  # geom_line(data = null_expectation_geometric,
  #           aes(x=death_percentage, y = expected_percentage, col = alive_group),
  #           linetype = "dashed", linewidth = 0.6) +
  geom_line(data = null_expectation_random_partner,
            aes(x=death_percentage, y = expected_percentage, col = alive_group),
            linetype = "dashed", linewidth = 0.6) +
  scale_color_brewer(palette = "Purples", name = "Nodes\nAlive at\nFragmentation?")+
  scale_fill_brewer(palette = "Purples", name = "Nodes\nAlive at\nFragmentation?")+
  labs(x="Death Probability", y="Mean Percentage")+
  NULL


#### Spatial constraints ####

# biophysics_frag=data.frame()
# 
# death_prob_value=c("0", "0.15", "0.3")
# death_prob_name=c("00", "15", "30")
# 
# for (i_cont in seq(length(death_prob_name))){
#   for (temp_cluster_selection in c("random")){ # "parent", "propagule", 
#     temp_death_prob_name=death_prob_name[i_cont]
#     temp_death_prob_value=death_prob_value[i_cont]
#     temp_df=read.csv(paste("fragmentation_simulations_spatial_constraints_test/prob_", temp_death_prob_name, "_", 
#                            temp_cluster_selection, "_only.csv", sep=""), header = TRUE)
#     
#     temp_df$variation=15
#     temp_df$death_percentage=temp_death_prob_value
#     temp_df$propagule_selection=str_to_sentence(temp_cluster_selection)
#     
#     biophysics_frag=rbind(biophysics_frag, temp_df)
#   }
# }
# 
# # extracting information from file_name column
# biophysics_frag=biophysics_frag %>%
#   mutate(
#     fname_clean = str_remove(file_name, "\\.csv$"),
#     parts = str_split(fname_clean, "_"),
#     diff_mean = as.numeric(sapply(parts, `[`, 3)),
#     percentage_diff = diff_mean/105*100,
#     sim_number_network_sim = as.numeric(sapply(parts, `[`, 6)),
#     cluster_generation = as.numeric(sapply(parts, `[`, 8)),
#     sim_number_physics_sim = sim_number
#   ) %>%
#   select(-parts, -fname_clean, -sim_number)
# 
# summary(biophysics_frag)
# 
# write.csv(biophysics_frag, "~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_11-14_cell_death_simulations/spatial_constraints_sim_8oct2026.csv", 
#           row.names=FALSE)

biophysics_frag=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_11-14_cell_death_simulations/spatial_constraints_sim_8oct2026.csv", 
                         header=TRUE)

summary(biophysics_frag)


# In how many cases is fragmentation triggered in the biophysics simulation before than the network simulation?

summ_biophysics_frag=biophysics_frag %>%
  group_by(death_percentage, percentage_diff, cluster_generation, propagule_selection) %>%
  summarise(prop_fragmentation = (n()-sum(is.na(fragmentation_size)))/n()*100)


  
ggplot(summ_biophysics_frag %>%
         mutate(percentage_diff=factor(as.character(percentage_diff), levels=c("50", "25", "0", "-25", "-50"))), 
       aes(x=cluster_generation, y=prop_fragmentation, col=propagule_selection))+
  geom_line()+
  facet_grid(percentage_diff~death_percentage)+
  theme_classic()+
  labs(x="Cluster Generation", 
       y="Percentage of Fragmentations in Biophysical Simulations\nTriggered before than in Network Simulations",
       col="Cluster\nSelected\nafter\nFragmentation")+
  NULL


# Around what size is fragmentation triggered according to the biophysical model?

max_network_sizes=biophysics_frag %>%
  group_by(death_percentage, percentage_diff, propagule_selection, cluster_generation) %>%
  filter(final_overlap<=85) %>%
  summarise(max_size=max(network_num_nodes))

ggplot(max_network_sizes%>%
         mutate(percentage_diff=factor(as.character(percentage_diff), levels=c("50", "25", "0", "-25", "-50"))), 
       aes(x=cluster_generation, y=max_size, col=propagule_selection))+
  geom_line()+
  facet_grid(percentage_diff~death_percentage)+
  theme_classic()+
  labs(x="Cluster Generation", 
       y="Maximum Cluster Size not Fractured\nin Physics Simulations",
       col="Cluster\nSelected\nafter\nFragmentation")+
  NULL


summary(biophysics_frag)

diff_fragmentation=biophysics_frag %>%
  group_by(death_percentage, percentage_diff, propagule_selection, cluster_generation) %>%
  mutate(diff_biophysics_networks=ifelse(is.na(fragmentation_size), 0, network_num_nodes-fragmentation_size))


ggplot(diff_fragmentation[diff_fragmentation$cluster_generation==20,]%>%
         mutate(percentage_diff=factor(as.character(percentage_diff), levels=c("50", "25", "0", "-25", "-50"))), 
       aes(x=propagule_selection, y=diff_biophysics_networks, fill=propagule_selection))+
  geom_violin()+
  facet_grid(percentage_diff~death_percentage)+
  theme_classic()+
  labs(x="Cluster Generation", 
       y="Difference in Size Prediction between\nNetwork and Biophysics Simulation",
       fill="Cluster\nSelected\nafter\nFragmentation")+
  scale_y_log10()+
  NULL
  
mean_diff_fragmentation=diff_fragmentation %>%
  group_by(death_percentage, percentage_diff, propagule_selection, cluster_generation)%>%
  summarise(mean_diff_biophysics_network_size=mean(diff_biophysics_networks))


ggplot(mean_diff_fragmentation[mean_diff_fragmentation$cluster_generation==20,] %>%
         mutate(death_percentage=factor(as.character(death_percentage), levels=c("0", "0.15", "0.3"))), 
       aes(x=percentage_diff, y=mean_diff_biophysics_network_size, shape=propagule_selection, 
           col=death_percentage, group=interaction(propagule_selection, death_percentage)))+
  geom_line()+
  geom_point()+
  # facet_grid(~death_percentage)+
  theme_classic()+
  labs(x="Delay (% Second Division)", 
       y="Mean Difference in Size Prediction between\nNetwork and Biophysics Simulation",
       shape="Cluster\nSelected\nafter\nFragmentation",
       col="Death\nProbability")+
  scale_y_log10()+
  scale_color_brewer(palette = "Purples")+
  NULL


diff_fragmentation_corrected=diff_fragmentation %>%
  mutate(fragmentation_size_biophysics=ifelse(is.na(fragmentation_size), network_num_nodes, fragmentation_size))


ggplot(diff_fragmentation_corrected[diff_fragmentation_corrected$cluster_generation==20,], 
       aes(x=fragmentation_size_biophysics, y=network_num_nodes, col=propagule_selection))+
  geom_point()+
  facet_grid(percentage_diff~death_percentage)+
  theme_classic()+
  scale_y_log10()+
  scale_x_log10()+
  labs(x="Fracture Size Biophysical Model", y="Fracture Size Network Model", 
       col="Cluster\nSelected\nafter\nFragmentation")+
  NULL



## To obtain the mean and confidence interval of the mean
# mean_frag_long <- diff_fragmentation_corrected %>%
#   filter(cluster_generation == 20) %>%
#   pivot_longer(
#     cols      = c(network_num_nodes, fragmentation_size_biophysics),
#     names_to  = "model",
#     values_to = "frag"
#   ) %>%
#   mutate(model = recode(model,
#                         network_num_nodes              = "Network",
#                         fragmentation_size_biophysics  = "Biophysic")) %>%
#   group_by(death_percentage, percentage_diff, propagule_selection, model) %>%
#   summarise(
#     mean_frag = mean(frag),
#     sd_frag   = sd(frag),
#     n         = n(),
#     se_frag   = sd_frag / sqrt(n),
#     ci_lower  = mean_frag - qt(0.975, n - 1) * se_frag,
#     ci_upper  = mean_frag + qt(0.975, n - 1) * se_frag,
#     .groups   = "drop"
#   )

## Median and confidence interval of the median
mean_frag_long <- diff_fragmentation_corrected %>%
  filter(cluster_generation == 20) %>%
  pivot_longer(
    cols      = c(network_num_nodes, fragmentation_size_biophysics),
    names_to  = "model",
    values_to = "frag"
  ) %>%
  mutate(model = recode(model,
                        network_num_nodes             = "Network",
                        fragmentation_size_biophysics = "Biophysic")) %>%
  group_by(death_percentage, percentage_diff, propagule_selection, model) %>%
  summarise(
    n           = n(),
    mean_frag = median(frag),
    sd_frag = sd(frag),
    se_frag = sd_frag / sqrt(n),
    q25 = quantile(frag, 0.25),
    q75 = quantile(frag, 0.75),
    k           = qbinom(0.025, n, 0.5),          # rank offset for a 95% CI
    ci_lower    = sort(frag)[max(k, 1)],
    ci_upper    = sort(frag)[min(n - k + 1, n)],
    .groups     = "drop"
  ) %>%
  select(-k)

ggplot(mean_frag_long, aes(x=percentage_diff, y=mean_frag, col=model, shape=propagule_selection))+
  geom_line()+
  geom_point()+
  facet_wrap(~death_percentage)+
  scale_y_log10()+
  labs(x="Delay", y="Fracture Size", col="Model", shape="Cluster\nSelected\nafter\nFragmentation")+
  # scale_shape_manual(values=c(1, 19))+
  # scale_color_manual(values=c("#9E8356FF", "#D1B79EFF"))+
  NULL



#### Supplementary Figure 13 ####

contour_diam=ggplot(summ_diameter[summ_diameter$num_nodes==1300,],
                    aes(x = death_percentage, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(breaks=seq(0.8, 1.7, 0.1)) +
  scale_fill_brewer(palette = "Purples", name = "Normalized\nNetwork\nDiameter")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_diam


contour_diam_not_norm=ggplot(summ_diameter_not_norm[summ_diameter_not_norm$num_nodes==1300,],
                             aes(x = death_percentage, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(breaks=seq(16, 34, 3)) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_diam_not_norm

contour_edge=ggplot(summ_max_edge[summ_max_edge$num_nodes==1300,],
                    aes(x = death_percentage, y = percentage_diff, z = mean_max_edge)) +
  geom_contour_filled(breaks = seq(13, 24, 1.2)) + # breaks = seq(13.5, 23.5, 1)
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_edge


contour_frag_size=ggplot(mean_mean_clust_size, aes(x = death_percentage, y = percentage_diff, z = mean)) +
  geom_contour_filled(bins=8)+ # breaks = c(2**seq(3, 12, 1))
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_frag_size

mean_mean_clust_size[mean_mean_clust_size$percentage_diff %in% c(-50, 0, 50) & 
                       mean_mean_clust_size$death_percentage %in% c(0, 0.3),]

contour_propagule_prop=ggplot(mean_mean_propagule_prop, aes(x = death_percentage, y = percentage_diff, z = mean)) +
  geom_contour_filled() +
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Propagule\nProportion")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_propagule_prop

contour_undivided_cells=ggplot(summ_undivided[summ_undivided$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_undivided)) +
  geom_contour_filled(breaks = seq(0, 180, 20)) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nDaughter\nCells")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_undivided_cells

contour_filament_motif=ggplot(summ_filamentous[summ_filamentous$num_nodes==1300,],
                                    aes(x = death_percentage, y = percentage_diff, z = mean_filament)) +
  geom_contour_filled(breaks = seq(0, 180, 20)) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_filament_motif


contour_motif_diff=ggplot(summ_motif_diff[summ_motif_diff$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_motif_difference_norm)) + #z=mean_motif_difference
  geom_contour_filled(breaks=seq(-0.10, 0.18, 0.04)) +
  scale_fill_brewer(palette = "Purples", name = "Mean\nMotif\nDifference\nNorm")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 9)+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  NULL
contour_motif_diff

## Not normalized diameter
# countour_plots_supp=plot_grid(contour_diam_not_norm, contour_edge,
#                               labels=c('A', 'B'), ncol=2, label_size=11)
# countour_plots_supp

countour_plots_supp=plot_grid(contour_diam_not_norm, contour_edge,
                              contour_frag_size, contour_propagule_prop,
                              labels=LETTERS[1:6], ncol=2, label_size = 10)
countour_plots_supp

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig13_death_prob_and_delay_17sep2026.png',
       plot=countour_plots_supp, dpi='retina', width=8.5, height=6, bg='white')





# Supp Fig 14 ####

size_bifurcation_random=ggplot(mean_frag_long[mean_frag_long$propagule_selection=='Random',], 
                               aes(x=percentage_diff, y=mean_frag, col=model, shape=model, fill=model)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, color=NA)+
  # geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.4, color=NA)+
  facet_wrap(~death_percentage) +
  scale_y_log10() +
  # labs(x="Delay", y="Fracture Size", col="Model", shape="Cluster\nSelected\nafter\nFragmentation") +
  labs(x="Delay (% Second Division)", y="Median Fracture Size", shape="Model", col="Model") +
  scale_shape_manual(values=c(19, 1)) +
  scale_color_manual(values=c("#833437FF", "#B7D9F2FF"))+
  scale_fill_manual(values=c("#833437FF", "#B7D9F2FF"))+
  guides(fill='none')+
  theme_classic(base_size = 9) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  NULL
size_bifurcation_random


size_bifurcation_random_apoptosis=ggplot(mean_frag_long[mean_frag_long$propagule_selection=='Random',], 
       aes(x=death_percentage, y=mean_frag, col=model, shape=model)) +
  geom_line() +
  geom_point() +
  facet_grid(~percentage_diff) +
  scale_y_log10() +
  # labs(x="Delay", y="Fracture Size", col="Model", shape="Cluster\nSelected\nafter\nFragmentation") +
  labs(x="Death Probability", y="Fracture Size", col="Model", shape="Model") +
  scale_shape_manual(values=c(19, 1)) +
  scale_color_manual(values=c("#833437FF", "#B7D9F2FF"))+
  theme_classic(base_size = 9) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  NULL
size_bifurcation_random_apoptosis


colnames(diff_fragmentation_corrected)

diff_frag_network=diff_fragmentation_corrected
diff_frag_network$mean_frag=diff_frag_network$network_num_nodes
diff_frag_network$model="Network"

diff_frag_biophysics=diff_fragmentation_corrected
diff_frag_biophysics$mean_frag=diff_frag_biophysics$fragmentation_size_biophysics
diff_frag_biophysics$model="Biophysic"

diff_frag_models_long=rbind(diff_frag_network, diff_frag_biophysics)

# Is apoptosis 0 from different from 30?
results_model_comparison_delay=data.frame()

for(temp_percentage_diff in unique(diff_frag_models_long$percentage_diff)){
  for(temp_model in c("Biophysic", "Network")){
    temp_df=diff_frag_models_long[diff_frag_models_long$death_percentage %in% c(0, 0.3) & 
                                    diff_frag_models_long$model==temp_model &
                                    diff_frag_models_long$percentage_diff==temp_percentage_diff,]
    
    temp_pvalue_wilcox=wilcox.test(data=temp_df, mean_frag ~ death_percentage, alternative = "less")$p.value
    temp_pvalue_t_test=t.test(data=temp_df, mean_frag ~ death_percentage, alternative = "less")$p.value
    temp_row=data.frame(model=temp_model,
                        percentage_diff=temp_percentage_diff,
                        p_value_wilcox=temp_pvalue_wilcox,
                        p_value_t_test=temp_pvalue_t_test)
    results_model_comparison_delay=rbind(results_model_comparison_delay, temp_row)
  }
}

results_model_comparison_delay

temp_df=diff_frag_models_long[diff_frag_models_long$death_percentage %in% c(0, 0.3) & 
                                diff_frag_models_long$model=='Biophysic' &
                                diff_frag_models_long$percentage_diff==-25,]

temp_df %>%
  group_by(model, death_percentage) %>%
  summarise(mean_size=mean(mean_frag))

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig14_apoptosis_comparison_models_19oct2026_median.png',
       plot=size_bifurcation_random, dpi='retina', width=5, height=3, bg='white')


# ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig13_apoptosis_comparison_models_15oct2026_apoptosis.png',
#        plot=size_bifurcation_random_apoptosis, dpi='retina', width=6.5, height=3, bg='white')


# Supp fig 11 ####


size_comp_no_apoptosis=ggplot(mean_frag_long[mean_frag_long$propagule_selection=='Random' & 
                                                mean_frag_long$death_percentage==0,], 
                               aes(x=percentage_diff, y=mean_frag, col=model, shape=model, fill=model)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.4, color=NA)+
  scale_y_log10() +
  labs(x="Delay (% Second Division)", y="Median Fracture Size", shape="Model", col="Model") +
  scale_shape_manual(values=c(19, 1)) +
  scale_color_manual(values=c("#833437FF", "#B7D9F2FF"))+
  scale_fill_manual(values=c("#833437FF", "#B7D9F2FF"))+
  guides(fill='none')+
  theme_classic(base_size = 9) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  NULL
size_comp_no_apoptosis



table(diff_frag_models_long$percentage_diff)

results_model_comparison_df=data.frame()

for(temp_percentage_diff in unique(diff_frag_models_long$percentage_diff)){
  for(temp_death_percentage in unique(diff_frag_models_long$death_percentage)){
    temp_df=diff_frag_models_long[diff_frag_models_long$death_percentage==temp_death_percentage & 
                                    diff_frag_models_long$percentage_diff==temp_percentage_diff &
                                    diff_frag_models_long$propagule_selection=='Random' & 
                                    diff_frag_models_long$cluster_generation==20,]
    
    temp_median_network=median(temp_df[temp_df$model=='Network',]$mean_frag)
    temp_median_biophysic=median(temp_df[temp_df$model=='Biophysic',]$mean_frag)
    
    temp_pvalue_wilcox=wilcox.test(data=temp_df, mean_frag ~ model)$p.value
    # temp_pvalue_t_test=t.test(data=temp_df, mean_frag ~ model)$p.value
    
    temp_row=data.frame(death_percentage=temp_death_percentage,
                        median_network=temp_median_network,
                        median_biophysic=temp_median_biophysic,
                        percentage_diff=temp_percentage_diff,
                        p_value_wilcox=temp_pvalue_wilcox)
                        # , p_value_t_test=temp_pvalue_t_test
    results_model_comparison_df=rbind(results_model_comparison_df, temp_row)
  }
}

results_model_comparison_df$not_significant=results_model_comparison_df$p_value_wilcox>(0.05/15)
results_model_comparison_df[results_model_comparison_df$death_percentage==0,]
# death prob = 0, 0.15 and 0.3; delays 0, 25 and 50 are not significant


pairwise.wilcox.test(diff_frag_models_long[diff_frag_models_long$model=="Biophysic" & 
                                                  diff_frag_models_long$percentage_diff<=0 & 
                                                  diff_frag_models_long$death_percentage==0,]$mean_frag,
                     diff_frag_models_long[diff_frag_models_long$model=="Biophysic" & 
                                             diff_frag_models_long$percentage_diff<=0 & 
                                             diff_frag_models_long$death_percentage==0,]$percentage_diff)
#     -50    -25   
# -25 <2e-16 -     
# 0   <2e-16 <2e-16

mean_frag_long[mean_frag_long$percentage_diff<=0 & 
                 mean_frag_long$death_percentage==0 &
                 mean_frag_long$propagule_selection=='Random',]

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig14_comparison_models_19oct2026_median.png',
       plot=size_comp_no_apoptosis, dpi='retina', width=4, height=3, bg='white')

