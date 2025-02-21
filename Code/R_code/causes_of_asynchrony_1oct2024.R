
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
# used to create 5 means and standard deviations
# delays=seq(default_mean, default_mean+120, 30)
# std_devs=seq(0,20,5)

#to create only 4 levels
# delays=seq(default_mean, default_mean+60, 20)
# std_devs=seq(0,15,5)

#to create 15 levels every 5 minutes
# delays=seq(default_mean, default_mean+75, 5)
# std_devs=seq(0,15,1)

#to create 15 levels every 2 minutes
delays=seq(default_mean, default_mean+30, 2)
std_devs=seq(0,15,1)

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

# write.csv(table_parameters, file="~/work_dir/observed_synchrony/evolution_results/syn_diff_mean_var_24apr2024/5_diffs_5_stds_16oct2024/params_5diffs_5stds_16oct2024.csv", row.names=FALSE)
# write.csv(table_parameters, file="~/work_dir/observed_synchrony/evolution_results/syn_diff_mean_var_24apr2024/4_diffs_4_stds_17oct2024/params_4diffs_4stds_17oct2024.csv", row.names=FALSE)
# write.csv(table_parameters, file="~/work_dir/observed_synchrony/evolution_results/syn_diff_mean_var_24apr2024/15_diffs_15_stds_17oct2024_every_5_minutes/params_15diffs_15stds_17oct2024.csv", row.names=FALSE)
# write.csv(table_parameters, file="~/work_dir/observed_synchrony/evolution_results/syn_diff_mean_var_24apr2024/15_diffs_15_stds_17oct2024_every_2_minutes/params_15diffs_15stds_2min_17oct2024.csv", row.names=FALSE)


#### Making plot of distributions ####

default_mean=60
means <- c(60, 80, 100, 120)
delay_names=c("0 min", "20 min", "40 min", "60 min")
std_devs <- c(0, 5, 10, 15)
std_names=c("0", "5", "10", "15")

n_times=1000

df_doub_times=data.frame()

for (i in seq(1,length(std_devs))){
  
  first_std=TRUE
  
  for (j in seq(1,length(means))){
    first_div_params=normal_to_lognormal_direct(means[j], means[j]/default_mean*std_devs[i])
    
    if(first_std){
      second_div_params=first_div_params
      first_std=FALSE
    }
    
    first_div_times=rlnorm(n_times, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
    
    second_div_times=rlnorm(n_times, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
    
    temp_df=data.frame(delay=delay_names[j], std=std_names[i], div_num=c(rep("1", n_times), rep("2", n_times)), 
                       doub_t=c(first_div_times, second_div_times))
    
    df_doub_times=rbind(df_doub_times, temp_df)
  }
}

df_doub_times$std=factor(df_doub_times$std, levels=c("0", "5", "10", "15"))

summary(df_doub_times)

ggplot(df_doub_times, aes(x=div_num, y=doub_t, fill=div_num))+
  geom_violin()+
  facet_grid(std~delay)+
  # stat_summary(fun='mean', geom='crossbar')+
  theme_classic()+
  labs(x="Number of Divisions", y="Doubling Time (min)")+
  guides(fill='none')+
  NULL



# Create the base plot
p <- ggplot(df_doub_times, aes(x=div_num, y=doub_t, fill=div_num)) +
  geom_violin() +
  facet_grid(std~delay) +
  theme_classic() +
  labs(x="Number of Divisions", y="Doubling Time (min)") +
  guides(fill='none') +
  theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))+ # Reduced plot margins
  scale_y_continuous(breaks = c(100, 200))+
  NULL

# Convert to grob
g <- ggplotGrob(p)

# Add top label (reduced height)
top_margin <- unit(0.5, "cm")  # Adjust this value to change the top gap
g <- gtable_add_rows(g, heights = top_margin, pos = 0)
g <- gtable_add_grob(g, 
                     grob = textGrob("Delay", gp = gpar(fontsize = 12, fontface = "bold")), 
                     t = 1, l = 3, r = ncol(g) - 1)

# Add right label (reduced width)
right_margin <- unit(0.5, "cm")  # Adjust this value to change the right gap
g <- gtable_add_cols(g, widths = right_margin, pos = -1)
g <- gtable_add_grob(g, 
                     grob = textGrob("Standard Deviation", rot = -90, gp = gpar(fontsize = 12, fontface = "bold")), 
                     t = 3, b = nrow(g) - 1, l = ncol(g))

# Draw the plot
grid.newpage()
grid.draw(g)



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
    
    df_temp=data.frame(delay=delay_j-default_mean, std=std_i, percentage=temp_ans, ld50=temp_median_diff,
                       first_div_mean=first_div_params$meanlog, first_div_sd=first_div_params$sdlog,
                       second_div_mean=second_div_params$meanlog, second_div_sd=second_div_params$sdlog)
    
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

p_test_delay=ggplot(divisions_diff_delay, aes(x=delay, y=ld50, col=as.factor(round(std,2)), group=as.factor(round(std,2))))+
  geom_line()+
  theme_classic()+
  labs(color="Std", x="Delay", y="Median Doubling\nTime Difference")+
  NULL
p_test_delay

p_test_var=ggplot(divisions_diff_var, aes(x=std, y=ld50, col=as.factor(delay), group=as.factor(delay)))+
  geom_line()+
  theme_classic()+
  labs(color="Delay", x="Standard Deviation", y="Median Doubling\nTime Difference")+
  NULL
p_test_var

plot_grid(p_test_delay, p_test_var)



# Paper Figure ####


figure_2_v3=plot_grid(g, p_test_delay, p_test_var,
                      labels=c('A)', 'B)', 'C)'), ncol=1, label_size=12, rel_heights=c(2.5,2,2))

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_causes_asynchrony_median_diff_20feb2025.svg',
       plot=figure_2_v3, dpi='retina', width=4.2, height=7, bg='white')
# Note: this image is saved as an svg to later modify the label of Delay to center align it correctly in the plot

# Supplementary figure of the contour map
ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig2_causes_async_contour_map_14feb2025.png',
       plot=contour_plot, width=4, height=3)

