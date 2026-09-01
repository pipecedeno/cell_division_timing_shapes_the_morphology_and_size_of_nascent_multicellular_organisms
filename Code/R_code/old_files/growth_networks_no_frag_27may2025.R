

#Date: 27may2025
# This code is used to analyze the results of growing the networks without fragmentation
# so it is only growing them until 1200 nodes in size

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

library(slider)


theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

setwd("~/work_dir/observed_synchrony/paper_results_2dec2024/growth_no_frag_27may2025/")

# Loading the data ####

# sim_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_300n_1200m/network_information.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   sim_df=rbind(sim_df, temp_df)
# }
# sim_df$strain <- ifelse(sim_df$strain=='petite', 'Petite', sim_df$strain)
# sim_df$strain <- ifelse(sim_df$strain=='petite-second-only', 'Petite w/o Delay', sim_df$strain)
# sim_df$strain <- ifelse(sim_df$strain=='grande', 'Grande', sim_df$strain)
# sim_df$strain=factor(sim_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# summary(sim_df)
# 
# write.csv(sim_df, file="growth_no_frag_diameter_27may2025.csv", row.names = FALSE)

sim_df=read.csv("growth_no_frag_diameter_27may2025.csv", header = TRUE)
sim_df$strain=factor(sim_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
summary(sim_df)


#### Smooth data calculation ####

# Function to calculate rolling statistics from raw data for multiple variables
calculate_rolling_stats <- function(df, window_size = 10) {
  num_nodes_vals <- sort(unique(df$num_nodes))
  
  map_dfr(num_nodes_vals, function(target_nodes) {
    # Define window boundaries
    lower <- target_nodes - window_size
    upper <- target_nodes + window_size
    
    # Get all raw data points within the window
    window_data <- df %>%
      filter(num_nodes >= lower & num_nodes <= upper)
    
    # Calculate statistics for each variable
    n_points <- nrow(window_data)
    t_val <- qt(0.975, df = max(n_points - 1, 1))
    
    # Diameter statistics
    mean_diameter <- mean(window_data$diameter, na.rm = TRUE)
    se_diameter <- sd(window_data$diameter, na.rm = TRUE) / sqrt(n_points)
    
    # Cases mother with undivided cells statistics
    mean_cases <- mean(window_data$cases_mother_with_undivided_cells, na.rm = TRUE)
    se_cases <- sd(window_data$cases_mother_with_undivided_cells, na.rm = TRUE) / sqrt(n_points)
    
    # Max edge degree statistics
    mean_edge <- mean(window_data$max_edge_degree, na.rm = TRUE)
    se_edge <- sd(window_data$max_edge_degree, na.rm = TRUE) / sqrt(n_points)
    
    tibble(
      num_nodes = target_nodes,
      # Diameter
      smooth_diameter = mean_diameter,
      smooth_diameter_se = se_diameter,
      smooth_diameter_lower = mean_diameter - t_val * se_diameter,
      smooth_diameter_upper = mean_diameter + t_val * se_diameter,
      # Cases mother with undivided cells
      smooth_cases = mean_cases,
      smooth_cases_se = se_cases,
      smooth_cases_lower = mean_cases - t_val * se_cases,
      smooth_cases_upper = mean_cases + t_val * se_cases,
      # Max edge degree
      smooth_edge = mean_edge,
      smooth_edge_se = se_edge,
      smooth_edge_lower = mean_edge - t_val * se_edge,
      smooth_edge_upper = mean_edge + t_val * se_edge,
      # Common
      smooth_n = n_points
    )
  })
}

# Apply to data
summ_smooth <- sim_df %>%
  group_by(strain) %>%
  group_modify(~ calculate_rolling_stats(.x, window_size = 10)) %>%
  ungroup()



#### Network Diameter

# Raw data
ggplot(sim_df, aes(x=num_nodes, y=diameter, col=strain, fill=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL

summ_diameter=sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(mean_diameter=mean(diameter))

ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=strain, fill=strain))+
  geom_line(alpha=0.75)+
  # geom_smooth()+
  petite_t200_colors+
  NULL


# Smooth data for diameter
ggplot(summ_smooth, aes(x = num_nodes, y = smooth_diameter, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_diameter_lower, ymax = smooth_diameter_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  labs(y = "Mean Diameter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL



# Undivided Daughters ####

#raw data
ggplot(sim_df, aes(x=num_nodes, y=cases_mother_with_undivided_cells, col=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL

summ_undivided=sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(mean_undivided=mean(cases_mother_with_undivided_cells))

ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=strain))+
  geom_line()+
  petite_t200_colors+
  NULL


# Smooth data for cases mother with undivided cells
ggplot(summ_smooth, aes(x = num_nodes, y = smooth_cases, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_cases_lower, ymax = smooth_cases_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  labs(y = "Mean Mothers with Undivided Nodes", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL



# Max Edge Degree ####

ggplot(sim_df, aes(x=num_nodes, y=max_edge_degree, col=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL

summ_max_edge=sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(mean_max_edge=mean(max_edge_degree))

ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=strain))+
  geom_line()+
  petite_t200_colors+
  NULL



# smooth data for max edge degree
ggplot(summ_smooth, aes(x = num_nodes, y = smooth_edge, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_edge_lower, ymax = smooth_edge_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  labs(y = "Mean max edge degree", x = "Number of nodes") +
  guides(col='none', fill='none')+
  NULL




# Paper Figure ####


# img_petite <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_diameter_undivided_23apr2025.png")
# img_plot_petite_net <- rasterGrob(img_petite, interpolate = TRUE)
# 
# img_grande <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_diameter_undivided_23apr2025.png")
# img_plot_grande_net <- rasterGrob(img_grande, interpolate = TRUE)
# 
# # Create text annotations
# text_petite <- textGrob("Petite", gp = gpar(fontsize = 14, fontface = "bold"))
# text_grande <- textGrob("Grande", gp = gpar(fontsize = 14, fontface = "bold"))
# 
# # Create ggplot objects for the images with annotations
# p_petite <- ggplot() + 
#   annotation_custom(img_plot_petite_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   annotation_custom(text_petite, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
#   theme_void() +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p_petite
# 
# p_grande <- ggplot() + 
#   annotation_custom(img_plot_grande_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   annotation_custom(text_grande, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
#   theme_void() +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p_grande



doubled=read.csv('~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/first_second_div_doubling_data_2024dec2.csv', header=TRUE)
doubled$division_number=factor(doubled$division_number)
doubled$timepoint=factor(doubled$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

petite_dt=doubled[doubled$strain=='petite' & as.numeric(doubled$division_number)<=2,]
petite_dt$name='Petite'

only_second_doubling=petite_dt[petite_dt$division_number==2,]
temp_first=only_second_doubling
temp_first$division_number=1
petite_second_doubling=rbind(only_second_doubling, temp_first)
petite_second_doubling$strain='second_doub'
petite_second_doubling$name='Petite w/o Delay'

petite_synthetic_data=rbind(petite_dt, petite_second_doubling)
summary(petite_synthetic_data)

grande_dt_temp=doubled[doubled$strain=='grande' & as.numeric(doubled$division_number)<=2,]
grande_dt_temp$name='Grande'

colnames(grande_dt_temp)

doub_t_data=rbind(petite_synthetic_data, grande_dt_temp)
doub_t_data$strain=factor(doub_t_data$strain, levels=c("petite", "second_doub", "grande"))
doub_t_data$name=factor(doub_t_data$name, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
table(doub_t_data$strain)

facet.labs=c('Ancestor', 'Ancestor \nw/o Delay', 'Evolved')
names(facet.labs)=c('Ancestor', 'Ancestor w/o Delay', 't200')



dt_distributions_v2=ggplot(doub_t_data, aes(x=division_number, y=hours, fill=name))+
  facet_wrap(~name)+
  geom_violin(adjust=2)+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none')+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  NULL
dt_distributions_v2


p_diam=ggplot(summ_smooth, aes(x = num_nodes, y = smooth_diameter, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_diameter_lower, ymax = smooth_diameter_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Diameter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_diam

p_undivided=ggplot(summ_smooth, aes(x = num_nodes, y = smooth_cases, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_cases_lower, ymax = smooth_cases_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Mothers with\nUndivided Nodes", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_undivided

p_edge_degree=ggplot(summ_smooth, aes(x = num_nodes, y = smooth_edge, col = strain, fill = strain)) +
  geom_ribbon(aes(ymin = smooth_edge_lower, ymax = smooth_edge_upper), alpha = 0.2, color = NA) +
  geom_line(alpha = 0.75) +
  petite_t200_colors +
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Max\nEdge Degree", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_edge_degree

fig=plot_grid(dt_distributions_v2, p_undivided, p_diam, p_edge_degree,
              labels=c('A', 'B', 'C', 'D'), ncol=2, label_size=12)
fig

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_3_growth_no_frag_28may2025.png',
       plot=fig, dpi='retina', width=6.5, height=4.5, bg='white')


