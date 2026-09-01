


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

theme_set(theme_classic(base_size = 16))

# death prob variable, constant delay ####

death_prob_name=sprintf("%02d", seq(0, 30, 5))
death_prob_value=as.character(seq(0, 0.3, 0.05))

delay_name=seq(-45, 45, length.out=21)
delay_percentage=delay_name/90*100

# only loading the last size of the networks
network_sizes=seq(250, 1250, 250)

# net_sim_df=data.frame()
# 
# for (j_cont in seq(length(delay_name))){
#   for (i_cont in seq(length(death_prob_value))){
#     temp_delay=delay_name[j_cont]
#     temp_percentages=delay_percentage[j_cont]
#     temp_death_prob=death_prob_value[i_cont]
#     temp_death_prob_name=death_prob_name[i_cont]
#     temp_df=read.csv(paste("prob_", temp_death_prob_name, "/test_15_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
# 
#     # temp_df=temp_df[temp_df$num_nodes==network_size,]
# 
#     temp_df$variation=15
#     temp_df$diff_mean=temp_delay
#     temp_df$percentage_diff=temp_percentages
#     temp_df$death_percentage=temp_death_prob
#     net_sim_df=rbind(net_sim_df, temp_df)
#   }
# }
# summary(net_sim_df)


# write.csv(net_sim_df, file="growth_no_frag_more_parameters_12sep2025.csv", row.names = FALSE)


net_sim_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_10-11_cell_death_simulations/growth_no_frag_more_parameters_12sep2025.csv", header=TRUE)

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
  mutate(mean_motif_difference_norm = (mean_filament - mean_undivided)/num_nodes)

summary(summ_motif_diff[summ_motif_diff$num_nodes>=10,]$mean_motif_difference_norm)


# Interaction of cell death and delay

ggplot(summ_motif_diff[summ_motif_diff$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_motif_difference_norm)) +
  geom_contour_filled(breaks=seq(-0.10, 0.18, 0.04)) +
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


# Supplementary plot synchronous strain ####

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


p_diam=ggplot(summ_diameter[summ_diameter$num_nodes %in% network_sizes & summ_diameter$percentage_diff==0,],
              aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Network Diameter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_diam

p_diam_not_norm=ggplot(summ_diameter_not_norm[summ_diameter_not_norm$num_nodes %in% network_sizes & summ_diameter_not_norm$percentage_diff==0,],
              aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Network Diameter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_diam_not_norm

p_undivided=ggplot(summ_undivided[summ_undivided$num_nodes %in% network_sizes & summ_undivided$percentage_diff==0,],
                   aes(x=death_percentage, y=mean_undivided, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Mothers with\nDelayed Daughter')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_undivided

p_filament=ggplot(summ_filamentous[summ_filamentous$num_nodes %in% network_sizes & summ_filamentous$percentage_diff==0,],
                  aes(x=death_percentage, y=mean_filament, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Filament Branches')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_filament

p_max_edge=ggplot(summ_max_edge[summ_max_edge$num_nodes %in% network_sizes & summ_max_edge$percentage_diff==0,],
                  aes(x=death_percentage, y=mean_max_edge, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  labs(x='Death Probability', y='Mean Max Edge Degree')+
  scale_color_brewer(palette="Purples", name="Network\nSize")+
  theme(legend.position = "none")+
  NULL
p_max_edge


# temporary plot to extract the legend from
temp_plot <- ggplot(summ_diameter[summ_diameter$num_nodes %in% network_sizes & summ_diameter$percentage_diff==0,],
                    aes(x=death_percentage, y=mean_diameter, col=as.factor(num_nodes), group=as.factor(num_nodes)))+
  geom_line()+
  geom_point()+
  theme_classic()+
  scale_color_brewer(palette="Purples", name="Network Size")+
  theme(legend.position = "bottom",           # Position at bottom
        legend.direction = "horizontal")      # Make it horizontal

# Extract the legend
shared_legend <- ggpubr::get_legend(temp_plot)

# # Create the plot grid without legend
# plots_grid <- plot_grid(p_diam, p_max_edge, p_undivided, p_filament, nrow = 2)
# 
# # Combine the plots with the shared legend
# final_plot <- plot_grid(plots_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
# final_plot

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig11_cell_death_18sep2025.svg',
#        plot=final_plot, dpi='retina', height=6, width=8, bg='white')


# Create the plot grid without legend
plots_grid <- plot_grid(p_diam_not_norm, p_max_edge, ncol = 2)

# Combine the plots with the shared legend
final_plot <- plot_grid(plots_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
final_plot

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig10_cell_death_28aug2026.svg',
       plot=final_plot, dpi='retina', height=4, width=8, bg='white')

#### Supplementary plot combined effect ####

contour_diam=ggplot(summ_diameter[summ_diameter$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(breaks=seq(0.8, 1.7, 0.1)) +
  scale_fill_brewer(palette = "Purples", name = "Normalized\nNetwork\nDiameter")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 10)+
  NULL
contour_diam


contour_diam_not_norm=ggplot(summ_diameter_not_norm[summ_diameter_not_norm$num_nodes==1300,],
                    aes(x = death_percentage, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled() +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 10)+
  NULL
contour_diam_not_norm

contour_edge=ggplot(summ_max_edge[summ_max_edge$num_nodes==1300,],
       aes(x = death_percentage, y = percentage_diff, z = mean_max_edge)) +
  geom_contour_filled(breaks = seq(13.5, 23.5, 1)) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Death Probability", y = "Delay (% Second Division)") +
  theme_classic(base_size = 10)+
  NULL
contour_edge

## Normalized diameter
# countour_plots_supp=plot_grid(contour_diam, contour_edge, 
#                               labels=c('A', 'B'), ncol=2, label_size=11)
# countour_plots_supp
# 
# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig12_death_prob_and_delay_18sep2025.png',
#        plot=countour_plots_supp, dpi='retina', width=9, height=3.5, bg='white')

## Not normalized diameter
countour_plots_supp=plot_grid(contour_diam_not_norm, contour_edge,
                              labels=c('A', 'B'), ncol=2, label_size=11)
countour_plots_supp

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig11_death_prob_and_delay_28aug2026.png',
       plot=countour_plots_supp, dpi='retina', width=9, height=3.5, bg='white')
