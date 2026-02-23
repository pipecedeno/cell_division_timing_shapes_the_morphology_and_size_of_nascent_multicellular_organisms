
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
library(see)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$YesterdayMedium)
syn_data_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("YesterdayMedium", direction = -1))

# color_diff_mean=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
#                                           breaks = delay_vector, limits = c(min(delay_vector), max(delay_vector))))
# color_variation=list(scale_color_gradient( name = "Variation", low = "pink", high = "darkred", 
#                                            breaks = std_vector, limits = c(min(std_vector), max(std_vector))))


setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig9_fast_first_div_7jul2025/")


# # Function for Direct Parameter Conversion
# #This formulas are in the wikipedia page of Log-normal distribution
# normal_to_lognormal_direct <- function(mean, sd) {
#   variance <- sd^2
#   meanlog <- log(mean^2 / sqrt(variance + mean^2))
#   sdlog <- sqrt(log(1 + variance / mean^2))
#   
#   return(list(meanlog = meanlog, sdlog = sdlog))
# }
# 
# # Create table of parameter ####
# 
# 
# default_mean=90
# default_std=15
# 
# # Create delays as percentages: -50% to +50% in steps of 5%
# delay_percentages=seq(-50, 50, 5)
# delays=default_mean * (1 + delay_percentages/100)
# 
# # for creating the plot of the distributions
# n_times=1000
# 
# df_doub_times=data.frame()
# table_parameters=data.frame()
# 
# for (i_cont in seq(length(delays))){
#   
#   delay_j=delays[i_cont]
#   temp_percentages=delay_percentages[i_cont]
#   
#   #calculating parameters of the first division
#   first_div_params=normal_to_lognormal_direct(delay_j, delay_j/default_mean*default_std)
#   
#   second_div_params=normal_to_lognormal_direct(default_mean, default_std)
#   
#   df_temp=data.frame(delay=delay_j-default_mean, std=default_std,
#                      first_div_mean=delay_j, first_div_scaled_std=delay_j/default_mean*default_std,
#                      first_div_log_mean=first_div_params$meanlog, first_div_log_sd=first_div_params$sdlog,
#                      second_div_mean=default_mean, second_div_scaled_std=default_std,
#                      second_div_log_mean=second_div_params$meanlog, second_div_log_sd=second_div_params$sdlog)
#   
#   table_parameters=rbind(table_parameters, df_temp)
#   
#   first_div_times=rlnorm(n_times, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
#   second_div_times=rlnorm(n_times, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
#   
#   temp_df=data.frame(delay=as.character(delay_j-default_mean), std=as.character(default_std), div_num=c(rep("1", n_times), rep("2", n_times)), 
#                      doub_t=c(first_div_times, second_div_times), percentage_diff=temp_percentages)
#   
#   df_doub_times=rbind(df_doub_times, temp_df)
#   
# }
# summary(table_parameters)
# 
# 
# # This files are used for simulations using this synthetic data distributions
# # write.csv(table_parameters, file="params_distributions_50percent_15var_7july2025.csv", row.names=FALSE)
# 
# 
# df_doub_times$percentage_diff=factor(df_doub_times$percentage_diff, levels=unique(df_doub_times$percentage_diff))
# ggplot(df_doub_times, aes(x=div_num, y=doub_t, fill=div_num))+
#   geom_violin()+
#   facet_wrap(~percentage_diff, nrow=1)+
#   # stat_summary(fun='mean', geom='crossbar')+
#   theme_classic()+
#   labs(x="Number of Divisions", y="Doubling Time (min)")+
#   guides(fill='none')+
#   NULL


# Function for Direct Parameter Conversion
# This formulas are in the wikipedia page of Log-normal distribution
normal_to_lognormal_direct <- function(mean, sd) {
  variance <- sd^2
  meanlog <- log(mean^2 / sqrt(variance + mean^2))
  sdlog <- sqrt(log(1 + variance / mean^2))
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}

# Create table of parameter ####
default_mean <- 90
default_std <- 15
# Create delays as percentages: -50% to +50% in steps of 5%
delay_percentages <- seq(-50, 50, 5)
delays <- default_mean * (1 + delay_percentages/100)

# Function to find x limits where PDF drops below threshold
find_pdf_limits <- function(meanlog, sdlog, threshold = 0.0001) {
  # Start with a reasonable range to search
  x_search_min <- 1
  x_search_max <- 1000
  
  # Find approximate limits using quantiles first (more efficient)
  # Use very small and large quantiles to get approximate bounds
  approx_min <- qlnorm(0.001, meanlog = meanlog, sdlog = sdlog)
  approx_max <- qlnorm(0.999, meanlog = meanlog, sdlog = sdlog)
  
  # Refine the search around these approximate values
  x_test_min <- seq(max(0.1, approx_min * 0.1), approx_min * 2, length.out = 100)
  x_test_max <- seq(approx_max * 0.5, approx_max * 3, length.out = 100)
  
  # Find where PDF drops below threshold on the left side
  pdf_min <- dlnorm(x_test_min, meanlog = meanlog, sdlog = sdlog)
  x_min <- min(x_test_min[pdf_min >= threshold])
  
  # Find where PDF drops below threshold on the right side  
  pdf_max <- dlnorm(x_test_max, meanlog = meanlog, sdlog = sdlog)
  x_max <- max(x_test_max[pdf_max >= threshold])
  
  return(list(x_min = x_min, x_max = x_max))
}


# For creating the plot of the PDF distributions
# pdf_threshold <- 0.000001  # Threshold for PDF cutoff
n_points <- 500  # Number of points to evaluate PDF
y_min=25
y_max=225

df_pdf_data <- data.frame()
table_parameters <- data.frame()

for (i_cont in seq(length(delays))) {
  
  delay_j <- delays[i_cont]
  temp_percentages <- delay_percentages[i_cont]
  
  # Calculating parameters of the first division
  first_div_params <- normal_to_lognormal_direct(delay_j, delay_j/default_mean*default_std)
  
  second_div_params <- normal_to_lognormal_direct(default_mean, default_std)
  
  df_temp <- data.frame(delay = delay_j - default_mean, std = default_std,
                        first_div_mean = delay_j, first_div_scaled_std = delay_j/default_mean*default_std,
                        first_div_log_mean = first_div_params$meanlog, first_div_log_sd = first_div_params$sdlog,
                        second_div_mean = default_mean, second_div_scaled_std = default_std,
                        second_div_log_mean = second_div_params$meanlog, second_div_log_sd = second_div_params$sdlog)
  
  table_parameters <- rbind(table_parameters, df_temp)
  
  # # Find limits for each distribution
  # first_limits <- find_pdf_limits(first_div_params$meanlog, first_div_params$sdlog, pdf_threshold)
  # second_limits <- find_pdf_limits(second_div_params$meanlog, second_div_params$sdlog, pdf_threshold)
  # 
  # # Use the broader range to encompass both distributions
  # overall_x_min <- min(first_limits$x_min, second_limits$x_min)
  # overall_x_max <- max(first_limits$x_max, second_limits$x_max)
  # 
  # # Create x values for PDF evaluation using the dynamic range
  # x_values <- seq(overall_x_min, overall_x_max, length.out = n_points)
  
  x_values=seq(y_min, y_max, length.out=n_points)
  
  # Calculate PDF values for both distributions
  first_div_pdf <- dlnorm(x_values, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
  second_div_pdf <- dlnorm(x_values, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
  
  # Create dataframe with PDF values
  temp_df <- data.frame(
    delay = as.character(delay_j - default_mean), 
    std = as.character(default_std), 
    div_num = c(rep("1", n_points), rep("2", n_points)), 
    doub_t = c(x_values, x_values),
    pdf_value = c(first_div_pdf, second_div_pdf),
    percentage_diff = temp_percentages
  )
  
  df_pdf_data <- rbind(df_pdf_data, temp_df)
}

summary(table_parameters)

# Create the violin plot using PDF data
# We need to "expand" the PDF data to create the violin effect
df_violin_data <- data.frame()

for (delay_val in unique(df_pdf_data$delay)) {
  for (div_val in unique(df_pdf_data$div_num)) {
    
    subset_data <- df_pdf_data[df_pdf_data$delay == delay_val & df_pdf_data$div_num == div_val, ]
    
    if (nrow(subset_data) > 0) {
      # Normalize PDF values to create reasonable violin width
      max_pdf <- max(subset_data$pdf_value)
      if (max_pdf > 0) {
        normalized_pdf <- subset_data$pdf_value / max_pdf
        
        # Create expanded data points based on PDF density
        # Higher PDF values get more points (creating wider violin sections)
        n_samples_per_point <- pmax(1, round(normalized_pdf * 50))  # Scale factor for violin width
        
        expanded_data <- data.frame()
        for (i in 1:nrow(subset_data)) {
          if (n_samples_per_point[i] > 0) {
            temp_expanded <- data.frame(
              delay = rep(delay_val, n_samples_per_point[i]),
              div_num = rep(div_val, n_samples_per_point[i]),
              doub_t = rep(subset_data$doub_t[i], n_samples_per_point[i]),
              percentage_diff = rep(subset_data$percentage_diff[i], n_samples_per_point[i])
            )
            expanded_data <- rbind(expanded_data, temp_expanded)
          }
        }
        
        df_violin_data <- rbind(df_violin_data, expanded_data)
      }
    }
  }
}

# Convert percentage_diff to factor with correct levels
df_violin_data$percentage_diff <- factor(df_violin_data$percentage_diff, 
                                         levels = unique(df_violin_data$percentage_diff))


ggplot(df_violin_data[df_violin_data$percentage_diff %in% c('-50', '-25', '0', '25', '50'),], 
       aes(x = div_num, y = doub_t, col = div_num)) +
  geom_violinhalf(trim = FALSE) +
  facet_wrap(~percentage_diff, nrow = 1) +
  theme_classic() +
  labs(x = "Number of Divisions", y = "Doubling Time (min)") +
  guides(col = 'none') +
  NULL




# No fragmentation sim data ####

# delay_values=as.character(round(delays,digits=1)-90)
# 
# sim_df=data.frame()
# 
# for (i_cont in seq(length(delay_values))){
#   temp_delay=delay_values[i_cont]
#   temp_percentages=delay_percentages[i_cont]
#   temp_df=read.csv(paste("no_frag_growth_50percent_15_var/test_15_var_", temp_delay, "_diff/network_information.csv", sep=""), header=TRUE)
#   
#   temp_df$variation=15
#   temp_df$diff_mean=temp_delay
#   temp_df$percentage_diff=temp_percentages
#   sim_df=rbind(sim_df, temp_df)
# }
# summary(sim_df)
# 
# 
# # write.csv(sim_df, file="growth_no_frag_inf_15_var_pencentages_7july2025.csv", row.names = FALSE)


sim_df=read.csv("growth_no_frag_inf_15_var_pencentages_7july2025.csv", header = TRUE)

summary(sim_df)

#### Network Diameter ####

summ_diameter <- sim_df %>%
  group_by(percentage_diff, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(diameter),
    sd_diameter = sd(diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL

ggplot(summ_diameter, aes(x = num_nodes, y = percentage_diff, fill = mean_diameter))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=18,
                       name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(summ_diameter, aes(x = num_nodes, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL




sim_df$norm_diameter=sim_df$diameter/(6.64385619*log10(sim_df$num_nodes)-1)

summ_norm_diameter <- sim_df %>%
  group_by(percentage_diff, variation, num_nodes) %>%
  summarise(
    mean_diameter = mean(norm_diameter),
    sd_diameter = sd(norm_diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

ggplot(summ_norm_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # scale_x_continuous(trans = 'log2')+
  # petite_t200_colors+
  # geom_vline(xintercept = 2**seq(1,10), linetype='dashed')+
  NULL

ggplot(summ_norm_diameter, aes(x=num_nodes, y=mean_diameter, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  scale_x_continuous(trans = 'log2')+
  # petite_t200_colors+
  geom_vline(xintercept = 2**seq(1,10), linetype='dashed')+
  NULL



ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
       aes(x = num_nodes, y = percentage_diff, fill = mean_diameter))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=1,
                       name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
       aes(x = num_nodes, y = percentage_diff, z = mean_diameter)) +
  geom_contour_filled(bins = 7) +
  scale_fill_brewer(palette = "Purples", name = "Network\nDiameter")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


#### Max edge degree ####

summ_max_edge <- sim_df %>%
  group_by(percentage_diff, variation, num_nodes) %>%
  summarise(
    mean_max_edge = mean(max_edge_degree),
    sd_max_edge = sd(max_edge_degree),
    n = n(),
    se_max_edge = sd_max_edge / sqrt(n),
    ci_lower = mean_max_edge - qt(0.975, n-1) * se_max_edge,
    ci_upper = mean_max_edge + qt(0.975, n-1) * se_max_edge,
    .groups = 'drop'
  )

ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL


ggplot(summ_max_edge, aes(x = num_nodes, y = percentage_diff, fill = mean_max_edge))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=15,
                       name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(summ_max_edge, aes(x = num_nodes, y = percentage_diff, z = mean_max_edge)) +
  geom_contour_filled(bins = 9) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL

#### Undivided cells ####

summ_undivided <- sim_df %>%
  group_by(percentage_diff, variation, num_nodes) %>%
  summarise(
    mean_undivided = mean(cases_mother_with_undivided_cells),
    sd_undivided = sd(cases_mother_with_undivided_cells),
    n = n(),
    se_undivided = sd_undivided / sqrt(n),
    ci_lower = mean_undivided - qt(0.975, n-1) * se_undivided,
    ci_upper = mean_undivided + qt(0.975, n-1) * se_undivided,
    .groups = 'drop'
  )

ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  NULL


ggplot(summ_undivided, aes(x = num_nodes, y = percentage_diff, fill = mean_undivided))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=50,
                       name = "Delayed\nCells")+
  labs(x = "Network Size", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(summ_undivided, aes(x = num_nodes, y = percentage_diff, z = mean_undivided)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nCells")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


#### Filamentous branches ####

summary(sim_df)

summ_filamentous <- sim_df %>%
  group_by(percentage_diff, variation, num_nodes) %>%
  summarise(
    mean_filament = mean(num_filamentous_branches),
    sd_filament = sd(num_filamentous_branches),
    n = n(),
    se_filament = sd_filament / sqrt(n),
    ci_lower = mean_filament - qt(0.975, n-1) * se_filament,
    ci_upper = mean_filament + qt(0.975, n-1) * se_filament,
    .groups = 'drop'
  )

ggplot(summ_filamentous, aes(x=num_nodes, y=mean_filament, col=as.factor(percentage_diff), fill=as.factor(percentage_diff)))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  # petite_t200_colors+
  geom_vline(xintercept = 100)+
  NULL


ggplot(summ_filamentous, aes(x = num_nodes, y = percentage_diff, fill = mean_filament))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=25,
                       name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(summ_filamentous, aes(x = num_nodes, y = percentage_diff, z = mean_filament)) +
  geom_contour_filled(bins = 8) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


#### Motif Difference ####

summ_motif_diff <- summ_filamentous %>%
  left_join(summ_undivided, 
            by = c("percentage_diff", "num_nodes"),
            suffix = c("_filamentous", "_undivided")) %>%
  mutate(mean_motif_difference = mean_filament - mean_undivided,
         mean_motif_difference_norm = (mean_filament - mean_undivided)/num_nodes)

summary(summ_motif_diff$mean_motif_difference)


ggplot(summ_motif_diff,
       aes(x = num_nodes, y = percentage_diff, z = mean_motif_difference)) +
  geom_contour_filled() +
  scale_fill_brewer(palette = "Purples", name = "Mean\nMotif\nDifference")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL

ggplot(summ_motif_diff,
       aes(x = num_nodes, y = percentage_diff, z = mean_motif_difference_norm)) +
  geom_contour_filled() +
  scale_fill_brewer(palette = "Purples", name = "Mean\nMotif\nDifference\nNorm")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL

# Fragmentation Size ####

# frag_df=data.frame()
# 
# for (i_cont in seq(length(delay_values))){
#   temp_delay=delay_values[i_cont]
#   temp_percentages=delay_percentages[i_cont]
#   for (temp_degree in as.character(seq(5, 15))){
#     temp_df=read.csv(paste("edge_degree_sim_50percent_15_var/test_15_var_", temp_delay,"_diff_", temp_degree,"_deg/fragmentation_inf.csv", sep=""), header=TRUE)
#     
#     temp_df$variation=15
#     temp_df$diff_mean=temp_delay
#     temp_df$percentage_diff=temp_percentages
#     temp_df$edge_degree=temp_degree
#     frag_df=rbind(frag_df, temp_df)
#   }
# }
# 
# summary(frag_df)
# 
# # write.csv(frag_df, file="edge_degree_sim_frag_15_var_percentages_7july2025.csv", row.names=FALSE)

frag_df=read.csv("edge_degree_sim_frag_15_var_percentages_7july2025.csv", header=TRUE)


summary(frag_df)

mean_clust_size=frag_df %>%
  group_by(variation, edge_degree, percentage_diff, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))



mean_mean_clust_size <- mean_clust_size %>%
  group_by(variation, edge_degree, percentage_diff) %>%
  summarize(mean = mean(mean_size),
            sd = sd(mean_size),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

summary(mean_mean_clust_size)

ggplot(mean_mean_clust_size, aes(x=edge_degree, y=mean, col=percentage_diff, group=percentage_diff)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  xlab('Edge Degree') +
  ylab('Mean fracture size') +
  # color_diff_mean+
  NULL


ggplot(mean_mean_clust_size, aes(x = edge_degree, y = percentage_diff, fill = mean))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=max(mean_mean_clust_size$mean)/3,
                       name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (% Second Division)")+
  theme_classic()+
  NULL

ggplot(mean_mean_clust_size, aes(x = edge_degree, y = percentage_diff, z = mean)) +
  geom_contour_filled(bins = 9) +
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


ggplot(mean_mean_clust_size, aes(x = edge_degree, y = percentage_diff, z = mean)) +
  geom_contour_filled(breaks = c(2**seq(3, 12, 1)))+
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (% Second Division)") +
  theme_classic()+
  NULL


# Paper Figures ####


img_sync <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/synchronous_div_combined_viz_delay_0min_9july2025.png")
img_plot_sync_net <- rasterGrob(img_sync, interpolate = TRUE)

img_fast_first <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/fast_first_div_combined_viz_delay_-30min_9july2025.png")
img_plot_fast_first_net <- rasterGrob(img_fast_first, interpolate = TRUE)

# Create text annotations
text_sync <- textGrob("Synchronous (0% Delay)", gp = gpar(fontsize = 10, fontface = "bold", col="black"))
text_fast_first <- textGrob("Fast First Division (-30% Delay)", gp = gpar(fontsize = 10, fontface = "bold", col='blue'))

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


# p_doub_times=ggplot(df_doub_times[df_doub_times$percentage_diff %in% c("-50", "-25", "0", "25", "50"),], 
#                     aes(x=div_num, y=doub_t, fill=div_num))+
#   geom_violin()+
#   facet_wrap(~percentage_diff, nrow=1)+
#   # stat_summary(fun='mean', geom='crossbar')+
#   theme_classic()+
#   labs(x="Division Number", y="Doubling Time (min)")+
#   guides(fill='none')+
#   NULL
# p_doub_times

p_doub_times=ggplot(df_violin_data[df_violin_data$percentage_diff %in% c('-50', '-25', '0', '25', '50'),], 
                    aes(x = div_num, y = doub_t, col = div_num)) +
  geom_violinhalf(trim = FALSE) +
  facet_wrap(~percentage_diff, nrow = 1) +
  theme_classic() +
  scale_color_manual(values=c("#BBA78CFF", "#333544FF"))+
  labs(x = "Number of Divisions", y = "Doubling Time (min)") +
  guides(col = 'none') +
  NULL
p_doub_times

p_diameter=ggplot(summ_norm_diameter[summ_norm_diameter$num_nodes>=10,], 
                  aes(x = num_nodes, y = percentage_diff, z = mean_diameter)) +
  # geom_contour_filled(bins = 7) +
  geom_contour_filled(breaks=seq(0.8, 1.4, 0.1)) +
  scale_fill_brewer(palette = "Purples", name = "Normalized\nNetwork\nDiameter")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  theme_classic()+
  NULL
p_diameter

p_edge_degree=ggplot(summ_max_edge[summ_max_edge$num_nodes>=5,], 
                     aes(x = num_nodes, y = percentage_diff, z = mean_max_edge)) +
  geom_contour_filled(breaks=seq(0, 24, 3)) +
  scale_fill_brewer(palette = "Purples", name = "Max Edge\nDegree")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  theme_classic()+
  NULL
p_edge_degree

p_delayed=ggplot(summ_undivided, aes(x = num_nodes, y = percentage_diff, z = mean_undivided)) +
  # geom_contour_filled(bins = 8) +
  # geom_contour_filled(breaks=seq(0, 135, 15)) +
  geom_contour_filled(breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  scale_fill_brewer(palette = "Purples", name = "Delayed\nDaughter\nCells")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  geom_point(data = data.frame(x = c(200, 200), y = c(0, -30)), 
             aes(x = x, y = y), 
             color = c("black", "blue"), 
             shape = 15,  # cross shape
             size = 2,   # adjust size as needed
             inherit.aes = FALSE) +  # don't inherit the z aesthetic
  theme_classic()+
  NULL
p_delayed

p_filament=ggplot(summ_filamentous, aes(x = num_nodes, y = percentage_diff, z = mean_filament)) +
  # geom_contour_filled(bins = 8) +
  geom_contour_filled(breaks = c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256)) +
  scale_fill_brewer(palette = "Purples", name = "Filament\nBranches")+
  labs(x = "Network Size", y = "Delay (% Second Division)") +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  geom_point(data = data.frame(x = c(200, 200), y = c(0, -30)), 
             aes(x = x, y = y), 
             color = c("black", "blue"), 
             shape = 15,  # cross shape
             size = 2,   # adjust size as needed
             inherit.aes = FALSE) +  # don't inherit the z aesthetic
  theme_classic()+
  NULL
p_filament

p_fragmentation=ggplot(mean_mean_clust_size, aes(x = edge_degree, y = percentage_diff, z = mean)) +
  geom_contour_filled(breaks = c(2**seq(3, 12, 1)))+
  # scale_fill_viridis_d(name = "Fracture\nSize") +  # Use discrete scale
  scale_fill_brewer(palette = "Purples", name = "Fracture\nSize")+
  labs(x = "Edge Degree", y = "Delay (% Second Division)") +
  scale_y_continuous(breaks = seq(-50, 50, 10)) +
  theme_classic()+
  NULL
p_fragmentation

# figure_first_div_net=plot_grid(p_doub_times, p_diameter, 
#                                p_edge_degree, p_delayed, 
#                                p_filament, p_sync_net,
#                                p_fast_first,
#                                p_fragmentation,
#                                labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), ncol=2, 
#                                align='hv', label_size=12)
# figure_first_div_net

figure_first_div_net=plot_grid(p_doub_times, p_diameter,
                               p_edge_degree, p_fragmentation,
                               p_filament, p_delayed,
                               p_fast_first, p_sync_net,
                               labels=c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'), ncol=2, 
                               align='hv', label_size=12)
figure_first_div_net

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_9_first_div_3oct2025.png',
       plot=figure_first_div_net, dpi='retina', height=12, width=10, bg='white')
