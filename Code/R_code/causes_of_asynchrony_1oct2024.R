
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
std_devs=seq(0,30,10)


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



