
# Date: 8Jun2024


library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(MASS)
library(grid)
library(cowplot)
library(grid)
library(gtable)
library(see)


# Function for Direct Parameter Conversion
#This formulas are in the wikipedia page of Log-normal distribution
normal_to_lognormal_direct <- function(mean, sd) {
  variance <- sd^2
  meanlog <- log(mean^2 / sqrt(variance + mean^2))
  sdlog <- sqrt(log(1 + variance / mean^2))
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}


# Estimating parameters for populations ####

default_mean=60
means <- c(60, 90, 120, 150)
std_devs <- c(0, 5, 10, 15)

results <- data.frame()

for (mu in means) {
  for (sigma in std_devs) {
    direct_params <- normal_to_lognormal_direct(mu, mu/default_mean*sigma)
    
    results <- rbind(results, 
                     data.frame(normal_mean = mu, 
                                scaled_sd = mu/default_mean*sigma, 
                                meanlog = direct_params$meanlog, 
                                sdlog = direct_params$sdlog))
  }
}

print(results)


# Create table of parameter ####


default_mean=60

#to create only 4 levels
delays=seq(default_mean, default_mean+60, 20)
std_devs=seq(0,15,5)


table_parameters=data.frame()

# iterating standard deviations
for (std_i in std_devs){
  
  first_sim_std=TRUE #flag to set second divisions distribution paramaters
  
  #iterating different delays
  for (delay_j in delays){
    
    #calculating parameters of the first division
    first_div_params=normal_to_lognormal_direct(delay_j, delay_j/default_mean*std_i)
    
    #saving second division distribution parameters and not changing it until the
    #std is changed
    if(first_sim_std){
      temp_second_div_mean=delay_j
      temp_second_div_var=std_i
      second_div_params=first_div_params
      first_sim_std=FALSE
    }
    
    df_temp=data.frame(delay=delay_j-default_mean, std=std_i,
                       first_div_mean=delay_j, first_div_scaled_std=delay_j/default_mean*std_i,
                       first_div_log_mean=first_div_params$meanlog, first_div_log_sd=first_div_params$sdlog,
                       second_div_mean=temp_second_div_mean, second_div_scaled_std=temp_second_div_var,
                       second_div_log_mean=second_div_params$meanlog, second_div_log_sd=second_div_params$sdlog)
    
    table_parameters=rbind(table_parameters, df_temp)
    
  }
}

summary(table_parameters)

# This files are used for simulations using this synthetic data distributions
# write.csv(table_parameters, file="~/work_dir/observed_synchrony/evolution_results/syn_diff_mean_var_24apr2024/4_diffs_4_stds_17oct2024/params_4diffs_4stds_17oct2024.csv", row.names=FALSE)







#### Making plot of PDF distributions ####
# default_mean <- 60
# means <- c(60, 80, 100, 120)
# delay_names <- c("0", "20", "40", "60")
# std_devs <- c(0, 5, 10, 15)
# std_names <- c("0", "5", "10", "15")
# delay_percentages=c("0", "33", "66", "100")
# 
# # Define range for PDF evaluation
# x_min <- 10  # Minimum doubling time to consider
# x_max <- 300 # Maximum doubling time to consider
# n_points <- x_max-x_min+1  # Number of points to evaluate PDF
# 
# df_pdf_data <- data.frame()
# 
# for (i in seq(1, length(std_devs))) {
#   
#   first_std <- TRUE
#   
#   for (j in seq(1, length(means))) {
#     first_div_params <- normal_to_lognormal_direct(means[j], means[j]/default_mean*std_devs[i])
#     
#     if(first_std) {
#       second_div_params <- first_div_params
#       first_std <- FALSE
#     }
#     
#     # Create x values for PDF evaluation
#     x_values <- seq(x_min, x_max, length.out = n_points)
#     
#     # Calculate PDF values for both distributions
#     first_div_pdf <- dlnorm(x_values, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
#     second_div_pdf <- dlnorm(x_values, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
#     
#     # Create dataframe with PDF values
#     temp_df <- data.frame(
#       delay = delay_names[j], 
#       std = std_names[i], 
#       delay_perc = delay_percentages[j],
#       div_num = c(rep("1", n_points), rep("2", n_points)), 
#       doub_t = c(x_values, x_values),
#       pdf_value = c(first_div_pdf, second_div_pdf)
#     )
#     
#     df_pdf_data <- rbind(df_pdf_data, temp_df)
#   }
# }
# 
# # Set factor levels
# df_pdf_data$std <- factor(df_pdf_data$std, levels = c("0", "5", "10", "15"))
# df_pdf_data$delay_perc <- factor(df_pdf_data$delay_perc, levels = c("0", "33", "66", "100"))
# 
# # Create the PDF plot
# p_pdf=ggplot(df_pdf_data, aes(x = doub_t, y = pdf_value, color = div_num)) +
#   geom_line(alpha = 0.75) +
#   facet_grid(std ~ delay_perc) +
#   theme_classic() +
#   xlim(c(0, 225))+
#   scale_y_continuous(breaks = c(0, 0.05)) +
#   scale_x_continuous(breaks = c(100, 200))+
#   coord_flip() +
#   labs(x = "Doubling Time (min)", y = "Probability Density", color = "Division") +
#   scale_color_discrete(name = "Number of\nDivisions") +
#   NULL
# p_pdf
# 
# # Convert to grob
# g_pdf <- ggplotGrob(p_pdf)
# 
# # Add top label (reduced height)
# top_margin <- unit(0.5, "cm")  # Adjust this value to change the top gap
# g_pdf <- gtable_add_rows(g_pdf, heights = top_margin, pos = 0)
# g_pdf <- gtable_add_grob(g_pdf, 
#                      grob = textGrob("Delay (% Second Division)", gp = gpar(fontsize = 11)), 
#                      t = 1, l = 3, r = ncol(g_pdf) - 1)
# 
# # Add right label (reduced width)
# right_margin <- unit(0.5, "cm")  # Adjust this value to change the right gap
# g_pdf <- gtable_add_cols(g_pdf, widths = right_margin, pos = -1)
# g_pdf <- gtable_add_grob(g_pdf, 
#                      grob = textGrob("Standard Deviation", rot = -90, gp = gpar(fontsize = 11)), 
#                      t = 3, b = nrow(g_pdf) - 1, l = ncol(g_pdf))
# 
# # Draw the plot
# grid.newpage()
# grid.draw(g_pdf)





# Your original parameters
default_mean <- 60
means <- c(60, 80, 100, 120)
delay_names <- c("0", "20", "40", "60")
std_devs <- c(0, 5, 10, 15)
std_names <- c("0", "5", "10", "15")
delay_percentages <- c("0", "33", "66", "100")

# Define range for PDF evaluation
x_min <- 10
x_max <- 300
n_points <- x_max - x_min + 1

# Create PDF data (modified to handle std = 0 case)
df_pdf_data <- data.frame()

for (i in seq(1, length(std_devs))) {
  
  first_std <- TRUE
  
  for (j in seq(1, length(means))) {
    
    # Special handling for std_dev = 0 case
    if (std_devs[i] == 0) {
      # When std = 0, the distribution is deterministic (single point)
      # Create just 3 points very close together to make an extremely thin violin
      
      # The exact value where all probability mass is located
      exact_value <- means[j]
      
      if(first_std) {
        second_div_params_exact <- exact_value
        first_std <- FALSE
      }
      
      # Create only 3 points: center and two very close neighbors
      tiny_offset <- 0.01  # Extremely small offset
      spike_x <- c(exact_value - tiny_offset, exact_value, exact_value + tiny_offset)
      spike_pdf <- c(0.001, 1, 0.001)  # Center has all the mass, neighbors almost nothing
      
      # Same for second division
      second_spike_x <- c(second_div_params_exact - tiny_offset, 
                          second_div_params_exact, 
                          second_div_params_exact + tiny_offset)
      second_spike_pdf <- c(0.001, 1, 0.001)
      
      # Create dataframe for the spike
      temp_df <- data.frame(
        delay = delay_names[j], 
        std = std_names[i], 
        delay_perc = delay_percentages[j],
        div_num = c(rep("1", 3), rep("2", 3)), 
        doub_t = c(spike_x, second_spike_x),
        pdf_value = c(spike_pdf, second_spike_pdf),
        is_deterministic = TRUE
      )
      
    } else {
      # Normal case with std > 0
      first_div_params <- normal_to_lognormal_direct(means[j], means[j]/default_mean*std_devs[i])
      
      if(first_std) {
        second_div_params <- first_div_params
        first_std <- FALSE
      }
      
      # Create x values for PDF evaluation
      x_values <- seq(x_min, x_max, length.out = n_points)
      
      # Calculate PDF values for both distributions
      first_div_pdf <- dlnorm(x_values, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
      second_div_pdf <- dlnorm(x_values, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
      
      # Create dataframe with PDF values
      temp_df <- data.frame(
        delay = delay_names[j], 
        std = std_names[i], 
        delay_perc = delay_percentages[j],
        div_num = c(rep("1", n_points), rep("2", n_points)), 
        doub_t = c(x_values, x_values),
        pdf_value = c(first_div_pdf, second_div_pdf),
        is_deterministic = FALSE
      )
    }
    
    df_pdf_data <- rbind(df_pdf_data, temp_df)
  }
}

# Set factor levels
df_pdf_data$std <- factor(df_pdf_data$std, levels = c("0", "5", "10", "15"))
df_pdf_data$delay_perc <- factor(df_pdf_data$delay_perc, levels = c("0", "33", "66", "100"))

# Create violin data by expanding PDF data points
df_violin_data <- data.frame()

for (std_val in unique(df_pdf_data$std)) {
  for (delay_val in unique(df_pdf_data$delay_perc)) {
    for (div_val in unique(df_pdf_data$div_num)) {
      
      subset_data <- df_pdf_data[df_pdf_data$std == std_val & 
                                   df_pdf_data$delay_perc == delay_val & 
                                   df_pdf_data$div_num == div_val, ]
      
      if (nrow(subset_data) > 0) {
        
        # Check if this is a deterministic case (std = 0)
        if (unique(subset_data$is_deterministic)) {
          # For deterministic case, create very few points with minimal expansion
          max_pdf <- max(subset_data$pdf_value)
          if (max_pdf > 0) {
            normalized_pdf <- subset_data$pdf_value / max_pdf
            
            # For deterministic cases, use very minimal expansion
            # Only expand the center point significantly
            n_samples_per_point <- ifelse(normalized_pdf >= 0.9, 50, 1)  # Only center gets many points
            
            expanded_data <- data.frame()
            for (k in 1:nrow(subset_data)) {
              if (n_samples_per_point[k] > 0) {
                temp_expanded <- data.frame(
                  std = rep(std_val, n_samples_per_point[k]),
                  delay_perc = rep(delay_val, n_samples_per_point[k]),
                  div_num = rep(div_val, n_samples_per_point[k]),
                  doub_t = rep(subset_data$doub_t[k], n_samples_per_point[k])
                )
                expanded_data <- rbind(expanded_data, temp_expanded)
              }
            }
            
            df_violin_data <- rbind(df_violin_data, expanded_data)
          }
          
        } else {
          # Normal case with std > 0
          max_pdf <- max(subset_data$pdf_value)
          if (max_pdf > 0) {
            normalized_pdf <- subset_data$pdf_value / max_pdf
            
            # Create expanded data points based on PDF density
            n_samples_per_point <- pmax(1, round(normalized_pdf * 50))
            
            expanded_data <- data.frame()
            for (k in 1:nrow(subset_data)) {
              if (n_samples_per_point[k] > 0) {
                temp_expanded <- data.frame(
                  std = rep(std_val, n_samples_per_point[k]),
                  delay_perc = rep(delay_val, n_samples_per_point[k]),
                  div_num = rep(div_val, n_samples_per_point[k]),
                  doub_t = rep(subset_data$doub_t[k], n_samples_per_point[k])
                )
                expanded_data <- rbind(expanded_data, temp_expanded)
              }
            }
            
            df_violin_data <- rbind(df_violin_data, expanded_data)
          }
        }
      }
    }
  }
}

# Ensure factor levels are maintained
df_violin_data$std <- factor(df_violin_data$std, levels = c("0", "5", "10", "15"))
df_violin_data$delay_perc <- factor(df_violin_data$delay_perc, levels = c("0", "33", "66", "100"))

# Create the violin plot
p_violin <- ggplot(df_violin_data, aes(x = div_num, y = doub_t, color = div_num)) +
  geom_violinhalf(trim = FALSE, alpha = 0.7) +
  facet_grid(std ~ delay_perc) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 225), breaks = c(100, 200)) +
  labs(x = "Number of Divisions", y = "Doubling Time (min)") +
  guides(col='none') +
  scale_color_discrete(name = "Number of\nDivisions") +
  theme(strip.text = element_text(size = 10), axis.text = element_text(size = 9)) +
  NULL
p_violin


g_violin <- ggplotGrob(p_violin)

# Add top label
top_margin <- unit(0.5, "cm")
g_violin <- gtable_add_rows(g_violin, heights = top_margin, pos = 0)
g_violin <- gtable_add_grob(g_violin, 
                            grob = textGrob("Delay (% Second Division)", gp = gpar(fontsize = 11)), 
                            t = 1, l = 3, r = ncol(g_violin) - 1)

# Add right label
right_margin <- unit(0.5, "cm")
g_violin <- gtable_add_cols(g_violin, widths = right_margin, pos = -1)
g_violin <- gtable_add_grob(g_violin, 
                            grob = textGrob("Standard Deviation", rot = -90, gp = gpar(fontsize = 11)), 
                            t = 3, b = nrow(g_violin) - 1, l = ncol(g_violin))

# Draw the final plot
grid.newpage()
grid.draw(g_violin)





# Heatmap of delay times ####

default_mean=60
delays=seq(default_mean, default_mean+60)
# delays=seq(default_mean, default_mean+120)

std_devs=seq(0,15,length.out=50)

n_differences=10000

# In what percentage of the cell cycle are we interested of quantifying the differences in doubling times?
percentage_cell_cycle=0.25

df_divisions_diff=data.frame()

# iterating standard deviations
for (std_i in std_devs){
  
  first_sim_std=TRUE #flag to set second divisions distribution paramaters
  
  #iterating different delays
  for (delay_j in delays){
    
    #calculating parameters of the first division
    first_div_params=normal_to_lognormal_direct(delay_j, delay_j/default_mean*std_i)
    
    #saving second division distribution parameters and not changing it until the
    #std is changed
    if(first_sim_std){
      second_div_params=first_div_params
      first_sim_std=FALSE
    }
    
    first_div_times=rlnorm(n_differences, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
    
    second_div_times=rlnorm(n_differences, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
    
    mean_dt=mean(c(first_div_times, second_div_times))
    
    norm_first_times=first_div_times/mean_dt
    norm_second_times=second_div_times/mean_dt
    
    difference=abs(norm_first_times-norm_second_times)
    
    temp_ans=mean(difference<percentage_cell_cycle)
    
    temp_median_diff=median(difference)
    
    temp_correlation=cor(first_div_times, second_div_times)
    
    df_temp=data.frame(delay=delay_j-default_mean, std=std_i, percentage=temp_ans, ld50=temp_median_diff,
                       first_div_mean=first_div_params$meanlog, first_div_sd=first_div_params$sdlog,
                       second_div_mean=second_div_params$meanlog, second_div_sd=second_div_params$sdlog,
                       correlation=temp_correlation)
    
    df_divisions_diff=rbind(df_divisions_diff, df_temp)
    
  }
}

summary(df_divisions_diff)

heat_map=ggplot(df_divisions_diff, aes(x = delay, y = std, fill = ld50))+
  geom_tile()+
  scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint=0, limits = c(0, 1),
                       name = "Median\nDifference")+
  labs(x="First Division Delay (min)", y="Doubling Time\nStandard Deviation")+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL
heat_map

contour_plot=ggplot(df_divisions_diff, aes(x = delay, y = std, z = ld50))+
  geom_contour_filled(breaks = seq(0, 1, by = 0.1))+
  scale_fill_brewer(palette = "Purples")+
  labs(x = "First Division Delay (min)", 
       y = "Doubling Time Standard Deviation",
       fill = "Median\nDoubling\nTime\nDifference")+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  NULL
contour_plot

contour_plot=ggplot(df_divisions_diff, aes(x = delay/60*100, y = std, z = ld50))+
  geom_contour_filled(breaks = seq(0, 1, by = 0.1))+
  scale_fill_brewer(palette = "Purples")+
  labs(x = "Delay (% Second Division)", 
       y = "Doubling Time Standard Deviation",
       fill = "Median\nDoubling\nTime\nDifference")+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  NULL
contour_plot


# test correlation plot

ggplot(df_divisions_diff, aes(x = delay/60*100, y = std, z = correlation))+
  geom_contour_filled()+
  scale_fill_brewer(palette = "Purples")+
  labs(x = "Delay (% Second Division)", 
       y = "Doubling Time Standard Deviation",
       fill = "Correlation")+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  NULL



summary(df_divisions_diff)

# df_divisions_diff[df_divisions_diff$delay%%20==0 & df_divisions_diff$delay%%5==0,]
divisions_diff_var=df_divisions_diff[df_divisions_diff$delay==0 | 
                                       df_divisions_diff$delay==20 |
                                       df_divisions_diff$delay==40 |
                                       df_divisions_diff$delay==60,]

divisions_diff_delay=df_divisions_diff[df_divisions_diff$std==0 | 
                                         df_divisions_diff$std==unique(df_divisions_diff$std)[18] |
                                         df_divisions_diff$std==unique(df_divisions_diff$std)[34] |
                                         df_divisions_diff$std==unique(df_divisions_diff$std)[50],]

p_test_delay=ggplot(divisions_diff_delay, aes(x=delay/60*100, y=ld50, col=as.factor(round(std,0)), group=as.factor(round(std,0))))+
  geom_line()+
  theme_classic()+
  labs(color="Std", x="Delay (% Second Division)", y="Median Doubling\nTime Difference")+
  NULL
p_test_delay

p_test_var=ggplot(divisions_diff_var, aes(x=std, y=ld50, col=as.factor(floor(delay/60*100)), group=as.factor(round(delay/60*100))))+
  geom_line()+
  theme_classic()+
  labs(color="Delay", x="Standard Deviation", y="Median Doubling\nTime Difference")+
  NULL
p_test_var

plot_grid(p_test_delay, p_test_var)



# Paper Figure ####


figure_2_v3=plot_grid(g_violin, p_test_delay, p_test_var,
                      labels=c('A', 'B', 'C'), ncol=1, label_size=12, rel_heights=c(2.5,2,2))
figure_2_v3

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_causes_asynchrony_median_diff_6aug2025.svg',
       plot=figure_2_v3, dpi='retina', width=4.2, height=7, bg='white')
# Note: this image is saved as an svg to later modify the label of Delay to center align it correctly in the plot

# Supplementary figure of the contour map
ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig2_causes_async_contour_map_30july2025.png',
       plot=contour_plot, width=4, height=3)

