
# Date: 

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
mycolors <- rev(ghibli_palettes$YesterdayMedium)
syn_data_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("YesterdayMedium", direction = -1))

# color_diff_mean=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
#                                           breaks = delay_vector, limits = c(min(delay_vector), max(delay_vector))))
# color_variation=list(scale_color_gradient( name = "Variation", low = "pink", high = "darkred", 
#                                            breaks = std_vector, limits = c(min(std_vector), max(std_vector))))


setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig9_fast_first_div_7jul2025/")


# Function for Direct Parameter Conversion
#This formulas are in the wikipedia page of Log-normal distribution
normal_to_lognormal_direct <- function(mean, sd) {
  variance <- sd^2
  meanlog <- log(mean^2 / sqrt(variance + mean^2))
  sdlog <- sqrt(log(1 + variance / mean^2))
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}

# Create table of parameter ####


default_mean=90
default_std=15

delays=seq(default_mean-30, default_mean+30, 5)
# delays=seq(default_mean-30, default_mean+30, 2)

# for creating the plot of the distributions
n_times=1000
df_doub_times=data.frame()


table_parameters=data.frame()


for (delay_j in delays){
  
  #calculating parameters of the first division
  first_div_params=normal_to_lognormal_direct(delay_j, delay_j/default_mean*default_std)
  
  second_div_params=normal_to_lognormal_direct(default_mean, default_std)
  
  df_temp=data.frame(delay=delay_j-default_mean, std=default_std,
                     first_div_mean=delay_j, first_div_scaled_std=delay_j/default_mean*default_std,
                     first_div_log_mean=first_div_params$meanlog, first_div_log_sd=first_div_params$sdlog,
                     second_div_mean=default_mean, second_div_scaled_std=default_std,
                     second_div_log_mean=second_div_params$meanlog, second_div_log_sd=second_div_params$sdlog)
  
  table_parameters=rbind(table_parameters, df_temp)
  
  first_div_times=rlnorm(n_times, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
  second_div_times=rlnorm(n_times, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
  
  temp_df=data.frame(delay=as.character(delay_j-default_mean), std=as.character(default_std), div_num=c(rep("1", n_times), rep("2", n_times)), 
                     doub_t=c(first_div_times, second_div_times))
  
  df_doub_times=rbind(df_doub_times, temp_df)
  
}


summary(table_parameters)

# This files are used for simulations using this synthetic data distributions
# write.csv(table_parameters, file="params_distributions_mean_90_7july2025.csv", row.names=FALSE)
# write.csv(table_parameters, file="params_distributions_mean_90_var_10_7july2025.csv", row.names=FALSE)
# write.csv(table_parameters, file="params_distributions_mean_90_var_15_7july2025.csv", row.names=FALSE)
# write.csv(table_parameters, file="params_distributions_mean_90_var_10_step_2_7july2025.csv", row.names=FALSE)

df_doub_times$delay=factor(df_doub_times$delay, levels=unique(df_doub_times$delay))
ggplot(df_doub_times, aes(x=div_num, y=doub_t, fill=div_num))+
  geom_violin()+
  facet_wrap(~delay, nrow=1)+
  # stat_summary(fun='mean', geom='crossbar')+
  theme_classic()+
  labs(x="Number of Divisions", y="Doubling Time (min)")+
  guides(fill='none')+
  NULL



# No fragmentation sim data ####

# delay_values=c("-30","-25","-20","-15","-10","-5","0","5","10","15","20","25","30")
delay_values=as.character(delays-90)

# sim_df=data.frame()
# 
# for (temp_delay in delay_values){
#   # temp_df=read.csv(paste("no_frag_growth_7july2025/test_5_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
#   # temp_df=read.csv(paste("no_frag_growth_10_var_7july2025/test_10_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
#   # temp_df=read.csv(paste("no_frag_growth_15_var_7july2025/test_15_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
#   # temp_df=read.csv(paste("no_frag_growth_10_var_2_step_7july2025/test_10_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
# 
#   temp_df$variation=10
#   temp_df$diff_mean=temp_delay
#   sim_df=rbind(sim_df, temp_df)
# }
# summary(sim_df)
# 
# # write.csv(sim_df, file="growth_no_frag_inf_7july2025.csv", row.names = FALSE)
# # write.csv(sim_df, file="growth_no_frag_inf_10_var_7july2025.csv", row.names = FALSE)
# # write.csv(sim_df, file="growth_no_frag_inf_15_var_7july2025.csv", row.names = FALSE)
# # write.csv(sim_df, file="growth_no_frag_inf_10_var_2_step_7july2025.csv", row.names = FALSE)

# sim_df=read.csv("growth_no_frag_inf_7july2025.csv", header = TRUE)
# sim_df=read.csv("growth_no_frag_inf_10_var_7july2025.csv", header = TRUE)
sim_df=read.csv("growth_no_frag_inf_15_var_7july2025.csv", header = TRUE)
# sim_df=read.csv("growth_no_frag_inf_10_var_2_step_7july2025.csv", header = TRUE)
summary(sim_df)

#### Network Diameter ####

summ_diameter <- sim_df %>%
  group_by(diff_mean, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(diameter),
    sd_diameter = sd(diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL

ggplot(summ_diameter, aes(x = num_nodes, y = diff_mean, fill = mean_diameter))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=18,
                       name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(summ_diameter, aes(x = num_nodes, y = diff_mean, z = mean_diameter)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (min)") +
  theme_classic()+
  NULL




sim_df$norm_diameter=sim_df$diameter/(6.64385619*log10(sim_df$num_nodes)-1)

summ_norm_diameter <- sim_df %>%
  group_by(diff_mean, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(norm_diameter),
    sd_diameter = sd(norm_diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

ggplot(summ_norm_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # scale_x_continuous(trans = 'log2')+
  # petite_t200_colors+
  # geom_vline(xintercept = 2**seq(1,10), linetype='dashed')+
  NULL

ggplot(summ_norm_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  scale_x_continuous(trans = 'log2')+
  # petite_t200_colors+
  geom_vline(xintercept = 2**seq(1,10), linetype='dashed')+
  NULL



ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
       aes(x = num_nodes, y = diff_mean, fill = mean_diameter))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=1,
                       name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
       aes(x = num_nodes, y = diff_mean, z = mean_diameter)) +
  geom_contour_filled(bins = 7) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (min)") +
  theme_classic()+
  NULL


#### Max edge degree ####

summ_max_edge <- sim_df %>%
  group_by(diff_mean, variation, num_nodes) %>%
  summarise(
    mean_max_edge = mean(max_edge_degree),
    sd_max_edge = sd(max_edge_degree),
    n = n(),
    se_max_edge = sd_max_edge / sqrt(n),
    ci_lower = mean_max_edge - qt(0.975, n-1) * se_max_edge,
    ci_upper = mean_max_edge + qt(0.975, n-1) * se_max_edge,
    .groups = 'drop'
  )

ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL


ggplot(summ_max_edge, aes(x = num_nodes, y = diff_mean, fill = mean_max_edge))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=15,
                       name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(summ_max_edge, aes(x = num_nodes, y = diff_mean, z = mean_max_edge)) +
  geom_contour_filled(bins = 9) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (min)") +
  theme_classic()+
  NULL

#### Undivided cells ####

summ_undivided <- sim_df %>%
  group_by(diff_mean, variation, num_nodes) %>%
  summarise(
    mean_undivided = mean(cases_mother_with_undivided_cells),
    sd_undivided = sd(cases_mother_with_undivided_cells),
    n = n(),
    se_undivided = sd_undivided / sqrt(n),
    ci_lower = mean_undivided - qt(0.975, n-1) * se_undivided,
    ci_upper = mean_undivided + qt(0.975, n-1) * se_undivided,
    .groups = 'drop'
  )

ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL


ggplot(summ_undivided, aes(x = num_nodes, y = diff_mean, fill = mean_undivided))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=50,
                       name = "Delayed\nCells")+
  labs(x = "Network Size", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(summ_undivided, aes(x = num_nodes, y = diff_mean, z = mean_undivided)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nCells")+
  labs(x = "Network Size", y = "Delay (min)") +
  theme_classic()+
  NULL


#### Filamentous branches ####

summary(sim_df)

summ_filamentous <- sim_df %>%
  group_by(diff_mean, variation, num_nodes) %>%
  summarise(
    mean_filament = mean(num_filamentous_branches),
    sd_filament = sd(num_filamentous_branches),
    n = n(),
    se_filament = sd_filament / sqrt(n),
    ci_lower = mean_filament - qt(0.975, n-1) * se_filament,
    ci_upper = mean_filament + qt(0.975, n-1) * se_filament,
    .groups = 'drop'
  )

ggplot(summ_filamentous, aes(x=num_nodes, y=mean_filament, col=as.factor(diff_mean), fill=as.factor(diff_mean)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  geom_vline(xintercept = 100)+
  NULL


ggplot(summ_filamentous, aes(x = num_nodes, y = diff_mean, fill = mean_filament))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=25,
                       name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(summ_filamentous, aes(x = num_nodes, y = diff_mean, z = mean_filament)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (min)") +
  theme_classic()+
  NULL


# Fragmentation Size ####

# frag_df=data.frame()
# 
# for (temp_delay in delay_values){
#   for (temp_degree in as.character(seq(5, 15))){
#     # temp_df=read.csv(paste("edge_degree_sim/test_5_var_", temp_delay,"_diff_", temp_degree,"_deg/fragmentation_inf.csv", sep=""), header=TRUE)
#     # temp_df=read.csv(paste("edge_degree_sim_10_var/test_10_var_", temp_delay,"_diff_", temp_degree,"_deg/fragmentation_inf.csv", sep=""), header=TRUE)
#     # temp_df=read.csv(paste("edge_degree_sim_15_var/test_15_var_", temp_delay,"_diff_", temp_degree,"_deg/fragmentation_inf.csv", sep=""), header=TRUE)
#     # temp_df=read.csv(paste("edge_degree_sim_10_var_2_step/test_10_var_", temp_delay,"_diff_", temp_degree,"_deg/fragmentation_inf.csv", sep=""), header=TRUE)
# 
#     temp_df$variation=10
#     temp_df$diff_mean=temp_delay
#     temp_df$edge_degree=temp_degree
#     frag_df=rbind(frag_df, temp_df)
#   }
# }
# 
# summary(frag_df)
# 
# # write.csv(frag_df, file="edge_degree_sim_frag_7july2025.csv", row.names=FALSE)
# # write.csv(frag_df, file="edge_degree_sim_frag_10_var_7july2025.csv", row.names=FALSE)
# # write.csv(frag_df, file="edge_degree_sim_frag_15_var_7july2025.csv", row.names=FALSE)
# # write.csv(frag_df, file="edge_degree_sim_frag_10_var_2_step_7july2025.csv", row.names=FALSE)

# frag_df=read.csv("edge_degree_sim_frag_7july2025.csv", header=TRUE)
# frag_df=read.csv("edge_degree_sim_frag_10_var_7july2025.csv", header=TRUE)
frag_df=read.csv("edge_degree_sim_frag_15_var_7july2025.csv", header=TRUE)
# frag_df=read.csv("edge_degree_sim_frag_10_var_2_step_7july2025.csv", header=TRUE)

# frag_df$edge_degree=factor(frag_df$edge_degree, levels=seq(5,15))
# frag_df$diff_mean=factor(frag_df$diff_mean, levels=delay_values)
summary(frag_df)

mean_clust_size=frag_df %>%
  group_by(variation, edge_degree, diff_mean, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))



mean_mean_clust_size <- mean_clust_size %>%
  group_by(variation, edge_degree, diff_mean) %>%
  summarize(mean = mean(mean_size),
            sd = sd(mean_size),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

summary(mean_mean_clust_size)

ggplot(mean_mean_clust_size, aes(x=edge_degree, y=mean, col=diff_mean, group=diff_mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  xlab('Edge Degree') +
  ylab('Mean fracture size') +
  # color_diff_mean+
  NULL


ggplot(mean_mean_clust_size, aes(x = edge_degree, y = diff_mean, fill = mean))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=max(mean_mean_clust_size$mean)/3,
                       limits = c(0, 1000),
                       name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (min)")+
  theme_classic()+
  NULL

ggplot(mean_mean_clust_size, aes(x = edge_degree, y = diff_mean, z = mean)) +
  geom_contour_filled(bins = 9) +
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (min)") +
  theme_classic()+
  NULL


ggplot(mean_mean_clust_size, aes(x = edge_degree, y = diff_mean, z = mean)) +
  geom_contour_filled(breaks = c(0, 2**seq(3, 10, 1)))+
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (min)") +
  theme_classic()+
  NULL


# Paper Figures ####


img_sync <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/synchronous_div_combined_viz_delay_0min_9july2025.png")
img_plot_sync_net <- rasterGrob(img_sync, interpolate = TRUE)

img_fast_first <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/fast_first_div_combined_viz_delay_-30min_9july2025.png")
img_plot_fast_first_net <- rasterGrob(img_fast_first, interpolate = TRUE)

# Create text annotations
text_sync <- textGrob("Synchronous", gp = gpar(fontsize = 10, fontface = "bold", col="red"))
text_fast_first <- textGrob("Fast First Division", gp = gpar(fontsize = 10, fontface = "bold", col='blue'))

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


p_doub_times=ggplot(df_doub_times[df_doub_times$delay %in% c("-30", "-14", "0", "14", "30"),], 
                    aes(x=div_num, y=doub_t, fill=div_num))+
  geom_violin()+
  facet_wrap(~delay, nrow=1)+
  # stat_summary(fun='mean', geom='crossbar')+
  theme_classic()+
  labs(x="Number of Divisions", y="Doubling Time (min)")+
  guides(fill='none')+
  NULL
p_doub_times


p_diameter=ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
                  aes(x = num_nodes, y = diff_mean, z = mean_diameter)) +
  # geom_contour_filled(bins = 7) +
  geom_contour_filled(breaks=seq(0.85, 1.25, 0.05)) +
  scale_fill_brewer(palette = "Purples", name = "Normalized\nNetwork\nDiameter")+
  labs(x = "Network Size", y = "Delay (min)") +
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  theme_classic()+
  NULL
p_diameter

p_edge_degree=ggplot(summ_max_edge[summ_max_edge$num_nodes>=5,], 
                     aes(x = num_nodes, y = diff_mean, z = mean_max_edge)) +
  geom_contour_filled(breaks = seq(0, 21, 3)) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (min)") +
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  theme_classic()+
  NULL
p_edge_degree

p_delayed=ggplot(summ_undivided, aes(x = num_nodes, y = diff_mean, z = mean_undivided)) +
  # geom_contour_filled(bins = 8) +
  # geom_contour_filled(breaks=seq(0, 135, 15)) +
  geom_contour_filled(breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nDaughter\nCells")+
  labs(x = "Network Size", y = "Delay (min)") +
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  geom_point(data = data.frame(x = c(200, 200), y = c(0, -30)), 
             aes(x = x, y = y), 
             color = c("red", "blue"), 
             shape = 15,  # cross shape
             size = 2,   # adjust size as needed
             inherit.aes = FALSE) +  # don't inherit the z aesthetic
  theme_classic()+
  NULL
p_delayed

p_filament=ggplot(summ_filamentous, aes(x = num_nodes, y = diff_mean, z = mean_filament)) +
  # geom_contour_filled(bins = 8) +
  geom_contour_filled(breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128)) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  # scale_fill_brewer(palette = "Spectral", name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (min)") +
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  geom_point(data = data.frame(x = c(200, 200), y = c(0, -30)), 
             aes(x = x, y = y), 
             color = c("red", "blue"), 
             shape = 15,  # cross shape
             size = 2,   # adjust size as needed
             inherit.aes = FALSE) +  # don't inherit the z aesthetic
  theme_classic()+
  NULL
p_filament

p_fragmentation=ggplot(mean_mean_clust_size, aes(x = edge_degree, y = diff_mean, z = mean)) +
  geom_contour_filled(breaks = c(0, 2**seq(3, 10, 1)))+
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (min)") +
  scale_y_continuous(breaks = seq(-30, 30, 10)) +
  theme_classic()+
  NULL
p_fragmentation

figure_first_div_net=plot_grid(p_doub_times, p_diameter, 
                               p_edge_degree, p_delayed, 
                               p_filament, p_sync_net,
                               p_fast_first,
                               p_fragmentation,
                               labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), ncol=2, 
                               align='hv', label_size=12)
figure_first_div_net

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/temp_fig_9_first_div_7july2025_15_var_2_step_networks_log.png',
       plot=figure_first_div_net, dpi='retina', height=12, width=10, bg='white')
