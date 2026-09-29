


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
setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/corrections_supreme_version_sep2026/supp_fig10_11_death_prob')

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

net_sim_df=data.frame()

for (j_cont in seq(length(delay_name))){
  for (i_cont in seq(length(death_prob_value))){
    temp_delay=delay_name[j_cont]
    temp_percentages=delay_percentage[j_cont]
    temp_death_prob=death_prob_value[i_cont]
    temp_death_prob_name=death_prob_name[i_cont]
    temp_df=read.csv(paste("growth_no_frag_results/prob_", temp_death_prob_name, "/test_15_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)

    temp_df=temp_df[temp_df$num_nodes %in% network_sizes,]

    temp_df$variation=15
    temp_df$diff_mean=temp_delay
    temp_df$percentage_diff=temp_percentages
    temp_df$death_percentage=temp_death_prob
    net_sim_df=rbind(net_sim_df, temp_df)
  }
}
summary(net_sim_df)

table(net_sim_df$num_nodes)


# write.csv(net_sim_df, file="growth_no_frag_death_prob_16sep2026.csv", row.names = FALSE)


net_sim_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_10-11_cell_death_simulations/growth_no_frag_death_prob_16sep2026.csv", header=TRUE)

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


# Supplementary Figure 10 ####

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
                     net_sim_df[net_sim_df$num_nodes==1250 & net_sim_df$percentage_diff==0,]$death_percentage)
#      0       0.05    0.1     0.15    0.2     0.25   
# 0.05 0.00024 -       -       -       -       -      
# 0.1  1.00000 0.00900 -       -       -       -      
# 0.15 1.00000 0.01361 1.00000 -       -       -      
# 0.2  0.00181 4.4e-06 0.14617 0.65271 -       -      
# 0.25 0.65271 0.00855 1.00000 1.00000 1.00000 -      
# 0.3  1.7e-05 9.7e-08 0.00446 0.03505 0.88322 0.14578


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

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig10_cell_death_16sep2026_ci.svg',
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

frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_10-11_cell_death_simulations/fragmentation_data_death_prob_16sep2026.csv", header=TRUE)
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



#### Supplementary Figure 11 ####

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
                              # contour_filament_motif, contour_undivided_cells,
                              labels=LETTERS[1:6], ncol=2, label_size = 10)
countour_plots_supp

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig11_death_prob_and_delay_17sep2026.png',
       plot=countour_plots_supp, dpi='retina', width=8.5, height=6, bg='white')



