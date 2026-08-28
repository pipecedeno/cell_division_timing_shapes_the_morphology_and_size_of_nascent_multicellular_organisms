
library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggridges)
library(ghibli)
# library(ggpubr)
# library(tidyverse)
library(ggnewscale)
library(cowplot)
library(png)
library(grid)
library(effsize)
library(ggh4x)
library(MASS)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

cell_sync=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/cell_sync_mother_daughter_2025may5.csv", header=TRUE)

cell_sync$ordered_timepoint=cell_sync$timepoint
cell_sync$ordered_timepoint=ifelse(cell_sync$strain=='petite', 'petite',
                                   ifelse(cell_sync$strain=='grande', 'grande',
                                          cell_sync$timepoint))
cell_sync$timepoint=factor(cell_sync$timepoint, levels=c('t0', 't200', 't400', 't600', 't1000'))
cell_sync$ordered_timepoint=factor(cell_sync$ordered_timepoint, levels=c('grande', 'petite', 't200', 't400', 't600', 't1000'))

unique(cell_sync[cell_sync$strain=='grande',]$id_file)

cell_sync=cell_sync[! cell_sync$id_file %in% c("2022nov22_gob8_6", "2022nov22_gob8_7", "2022nov23_gob8_7"),]

summary(cell_sync)



ggplot(cell_sync, aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.1)+
  facet_wrap(~ordered_timepoint)+
  geom_abline()+
  xlab('Mother doubling time (min)')+
  ylab('Daughter doubling time (min)')+
  NULL


# Branch synchrony ####

df_all=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/timelapse_doubling_time_inf_2025may5.csv", header = TRUE)

unique(df_all[df_all$strain=='grande',]$id_file)
# "2024nov5_gob8_4"   "2022nov22_gob8_6"  "2024oct31_gob8_2"  "2024oct30_gob8_11" "2022nov22_gob8_7"  "2024nov5_gob8_12" "2022nov23_gob8_7"  "2024oct30_gob8_6"

df_all=df_all[! df_all$id_file %in% c("2022nov22_gob8_6", "2022nov22_gob8_7", "2022nov23_gob8_7"),]

df_all$division_number=as.factor(df_all$division_number)
df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

summary(df_all)


doubled=df_all[df_all$N_PREDECESSORS==1 & df_all$N_SUCCESSORS==2,]


table(doubled[doubled$strain=='grande',]$division_number)
#   1   2   3   4   5   6   7   8 
# 147  53   3   0   0   0   0   0

table(doubled[doubled$strain=='petite',]$division_number)
#   1   2   3   4   5   6   7   8 
# 135  48   6   2   0   0   0   0 

# normalizing the frames ####
# The idea is to make the doubling time of the cell 
normal_frames=data.frame(matrix(ncol = dim(doubled)[2], nrow = 0))
colnames(normal_frames)=colnames(doubled)
for (i in unique(doubled$unique_id_branch)){
  temp_df=doubled[doubled$unique_id_branch==i & doubled$FRAME!=180,]
  
  if(dim(temp_df)[1]>2){
    #normalizing the FRAME value by making it start at 0
    first_frame=min(temp_df$FRAME)
    temp_df$FRAME=temp_df$FRAME-first_frame
    normal_frames=rbind(normal_frames,temp_df)
  }
}


# Selected time lapse plots ####

summary(normal_frames)

table(normal_frames[normal_frames$strain_timepoint=='petite_t0',]$id_file)

ggplot(normal_frames[normal_frames$id_file=='2022nov22_gob440_2',], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram()+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL


ggplot(normal_frames[normal_frames$id_file=='2022sep30_gob21_6',], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram()+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL




# Relationship between delay and SD ####

# Function for Direct Parameter Conversion
#This formulas are in the wikipedia page of Log-normal distribution
normal_to_lognormal_direct <- function(mean, sd) {
  variance <- sd^2
  meanlog <- log(mean^2 / sqrt(variance + mean^2))
  sdlog <- sqrt(log(1 + variance / mean^2))
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}


# Loading percentage of cells in their first division according to python numerical simulations
# first_div_percentages=read.csv('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/percentage_cells_in_first_division_23sep2025.csv', header=TRUE)


# Parameters
default_mean <- 60
delays <- seq(default_mean, default_mean + 60, 3)
std_devs <- seq(0, 30, length.out = 21)
n_samples <- 10000

# DATASET 1: Fixed variation at 0, iterate through delays
df_delay_effect <- data.frame()

fixed_variation <- 0  # No variation in synchrony

for (delay_j in delays) {
  
  # First division parameters (with delay)
  first_div_params <- normal_to_lognormal_direct(delay_j, fixed_variation)
  
  # Second division parameters (no delay, no variation)
  second_div_params <- normal_to_lognormal_direct(default_mean, fixed_variation)
  
  # Commenting this to not use the calculations of the number of cells in each division number from the python calculations
  # first_div_samples=(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100)*(n_samples*2)
  # second_div_samples=(1-(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100))*(n_samples*2)
  
  first_div_samples=n_samples
  second_div_samples=n_samples
  
  # Generate samples
  first_div_times <- rlnorm(first_div_samples, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
  first_div_times[first_div_times < 0] <- 0
  second_div_times <- rlnorm(second_div_samples, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
  second_div_times[second_div_times < 0] <- 0
  
  # Combine both distributions
  combined_times <- c(first_div_times, second_div_times)
  
  # Calculate standard deviation of combined distribution
  combined_sd <- sd(combined_times)
  
  # Store results
  df_temp <- data.frame(
    delay = delay_j - default_mean,
    variation = fixed_variation,
    combined_sd = combined_sd,
    first_div_mean = mean(first_div_times),
    second_div_mean = mean(second_div_times),
    dataset = "delay_effect"
  )
  
  df_delay_effect <- rbind(df_delay_effect, df_temp)
}

# DATASET 2: Fixed delay at 0, iterate through variations
df_variation_effect <- data.frame()

fixed_delay <- default_mean  # No delay

for (std_i in std_devs) {
  
  # First division parameters (no delay, with variation)
  first_div_params <- normal_to_lognormal_direct(fixed_delay, std_i)
  
  # Second division parameters (no delay, with variation)
  second_div_params <- normal_to_lognormal_direct(fixed_delay, std_i)
  
  # Generate samples
  first_div_times <- rlnorm(n_samples, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
  first_div_times[first_div_times < 0] <- 0
  second_div_times <- rlnorm(n_samples, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
  second_div_times[second_div_times < 0] <- 0
  
  # Combine both distributions
  combined_times <- c(first_div_times, second_div_times)
  
  # Calculate standard deviation of combined distribution
  combined_sd <- sd(combined_times)
  
  # Store results
  df_temp <- data.frame(
    delay = fixed_delay - default_mean,
    variation = std_i,
    combined_sd = combined_sd,
    first_div_mean = mean(first_div_times),
    second_div_mean = mean(second_div_times),
    dataset = "variation_effect"
  )
  
  df_variation_effect <- rbind(df_variation_effect, df_temp)
}

# Combine both datasets for easier analysis
df_combined_analysis <- rbind(df_delay_effect, df_variation_effect)




# Combined plot for comparison
ggplot(df_combined_analysis, aes(x = ifelse(dataset == "delay_effect", delay*100/60, variation), 
                                 y = combined_sd, 
                                 color = dataset)) +
  geom_line() +
  geom_point() +
  labs(title = "Comparison of Delay vs Standard Deviation",
       x = "Parameter Value (delay % & variation in SD)",
       y = "Combined Standard Deviation",
       color = "Effect Type") +
  theme_classic()+
  NULL


ggplot(df_combined_analysis, aes(x = ifelse(dataset == "delay_effect", delay*100/60, variation*(10/3)), 
                                 y = combined_sd, 
                                 color = dataset)) +
  geom_line() +
  geom_point() +
  labs(title = "Comparison of Delay vs Standard Deviation*3.33",
       x = "Parameter Value (delay % & variation in SD*3.33)",
       y = "Combined Standard Deviation",
       color = "Effect Type") +
  theme_classic()+
  NULL


relationship_df=data.frame(sd_variation=df_combined_analysis[df_combined_analysis$dataset=='variation_effect',]$variation,
                           sd_delay=df_combined_analysis[df_combined_analysis$dataset=='delay_effect',]$combined_sd)

ggplot(relationship_df, aes(x=sd_delay, y=sd_variation))+
  geom_line(color='lightblue')+
  geom_point(color='blue')+
  theme_classic()+
  labs(x='Combined Standard Deviation Effect by Delay',
       y='Doubling Time Distribution Standard Deviation')+
  NULL



# Interaction between delay and standard deviation ####


df_grid_analysis <- data.frame()

delays_grid <- seq(0, 60, 3)  # 0 to 60 minutes delay
std_devs_grid <- seq(0, 30, length.out = 31)  # 0 to 30 SD units

# Create grid of all combinations
for (delay_i in delays_grid) {
  for (std_j in std_devs_grid) {
    
    # First division parameters (with both delay and variation)
    first_div_params <- normal_to_lognormal_direct(default_mean + delay_i, std_j)
    
    # Second division parameters (no delay, with variation)
    second_div_params <- normal_to_lognormal_direct(default_mean, std_j)
    
    # Commenting this to not use the calculations of the number of cells in each division number from the python calculations
    # Obtain percentage of cells in the first cell division and the rest
    # first_div_samples=(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100)*(n_samples*2)
    # second_div_samples=(1-(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100))*(n_samples*2)
    
    first_div_samples=n_samples
    second_div_samples=n_samples
    
    # Generate samples
    first_div_times <- rlnorm(first_div_samples, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
    first_div_times[first_div_times < 0] <- 0
    
    second_div_times <- rlnorm(second_div_samples, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
    second_div_times[second_div_times < 0] <- 0
    
    # Combine both distributions
    combined_times <- c(first_div_times, second_div_times)
    
    # Calculate standard deviation of combined distribution
    combined_sd <- sd(combined_times)
    
    # Difference in doubling times
    # Resampling division times to avoid an issue of different sample sizes
    first_div_times <- rlnorm(n_samples, meanlog = first_div_params$meanlog, sdlog = first_div_params$sdlog)
    first_div_times[first_div_times < 0] <- 0
    
    second_div_times <- rlnorm(n_samples, meanlog = second_div_params$meanlog, sdlog = second_div_params$sdlog)
    second_div_times[second_div_times < 0] <- 0
    
    median_difference=median(first_div_times-second_div_times)
    mean_difference=mean(first_div_times-second_div_times)
    
    # Store results
    df_temp <- data.frame(
      delay = delay_i,
      std = std_j,
      combined_sd = combined_sd,
      first_div_mean = mean(first_div_times),
      second_div_mean = mean(second_div_times),
      median_diff=median_difference,
      mean_diff=mean_difference
    )
    
    df_grid_analysis <- rbind(df_grid_analysis, df_temp)
  }
  
  # Progress indicator
  if (delay_i %% 10 == 0) {
    cat("Completed delay =", delay_i, "minutes\n")
  }
}

# Create the contour plot
ggplot(df_grid_analysis, aes(x = delay*100/60, y = std, z = combined_sd)) +
  geom_contour_filled(breaks=seq(0, 45, 5)) +  # You can adjust bins or use breaks = seq(min, max, by = step)
  # scale_fill_viridis_d(name = "Combined\nStandard\nDeviation") +  # Using viridis, but you can change to scale_fill_brewer()
  scale_fill_brewer(palette = "Purples", name = "Combined\nStandard\nDeviation")+
  labs(x = "Delay (% Second Division)", 
       y = "Doubling Time Standard Deviation") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL



# Figure 2 ####


# img_petite_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_2hr_5f_scale_masks_v2.png")
img_petite_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_45min-step_50-86f_modified.png")
img_plot_petite_v2 <- rasterGrob(img_petite_v2, interpolate = TRUE)

# img_grande_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_49f_1hr_scale_1x_masks_complete.png")
img_grande_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_30min-step_start-2hr_24-48f_modified.png")
img_plot_grande_v2 <- rasterGrob(img_grande_v2, interpolate = TRUE)


ancestor_tl_dt=normal_frames[normal_frames$strain=='petite' | normal_frames$strain=='grande',]
table(ancestor_tl_dt$strain)
ancestor_tl_dt$strain=factor(ancestor_tl_dt$strain, levels=c("grande", "petite"))


p2_ancestors_v2 = ggplot(ancestor_tl_dt, 
                         aes(x = 5*FRAME, fill = strain)) +
  xlim(c(-15, 350))+
  geom_histogram(binwidth = 10, position = 'identity')+
  facet_wrap(~strain, ncol=1, scales='free')+
  xlab('Time from First Division (min)')+
  ylab('Number of Divisions')+
  labs(fill = "Strain")+
  scale_fill_manual(labels = c("Grande", "Petite"),
                    values = c("#AE93BEFF", "#F0D77BFF"))+
  theme_classic(base_size = 10) +  # Set base theme with font size FIRST
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        legend.key.size = unit(0.2, 'cm'),      # Size of the colored squares
        legend.text = element_text(size = 8),   # Size of legend text
        legend.title = element_text(size = 8))+ # Size of legend title) +  
  NULL
p2_ancestors_v2

table(ancestor_tl_dt$strain)


ks.test(ancestor_tl_dt[ancestor_tl_dt$strain == 'petite', ]$FRAME,
        ancestor_tl_dt[ancestor_tl_dt$strain == 'grande', ]$FRAME)
# D = 0.10345, p-value = 0.2711

wilcox.test(FRAME * 5 ~ strain, data = ancestor_tl_dt) # p-value = 0.2806



ggplot(ancestor_tl_dt, aes(x = 5*FRAME, color = strain)) +
  stat_ecdf(linewidth = 1) +
  xlim(c(-15, 350)) +
  xlab('Time from First Division (min)') +
  ylab('Cumulative Probability') +
  labs(color = "Strain") +
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.2, 0.8),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  NULL



cell_sync_ancestors=cell_sync[cell_sync$strain=='grande' | cell_sync$strain=='petite',]
cell_sync_ancestors$strain=ifelse(cell_sync_ancestors$strain=='petite', 'Petite', 'Grande')

table(cell_sync_ancestors$strain)
# Grande Petite 
# 53     46

p3_ancestors=ggplot(cell_sync_ancestors, 
                    aes(x=mother_minutes/60, y=daughter_minutes/60, col=strain))+
  geom_point(alpha=0.75)+
  geom_abline(linetype='dashed')+
  # geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, linewidth = 0.7) +
  facet_wrap(~strain, ncol=5)+
  guides(col='none')+
  xlab('Mother Doubling Time (hours)')+
  ylab('Daughter Doubling\nTime (hours)')+
  # scale_x_continuous(breaks = seq(100, 300, by = 100))+
  scale_color_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  theme_classic(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p3_ancestors


rlm_results=rlm(daughter_minutes ~ mother_minutes, data=cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',])

lm_results=lm(daughter_minutes ~ mother_minutes, data=cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',])


ggplot(cell_sync_ancestors, 
       aes(x=daughter_minutes-mother_minutes, col=strain))+
  geom_histogram()+
  facet_wrap(~strain)+
  NULL


#signed rank wilcoxon test
result_petite=wilcox.test(y=cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',]$daughter_minutes, x=cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',]$mother_minutes, paired=TRUE)
result_petite
# V = 101, p-value = 2.624e-06
median(cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',]$daughter_minutes-cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',]$mother_minutes)
# 35

Z <- qnorm(result_petite$p.value / 2)
n_pairs <- nrow(cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',])
abs(Z) / sqrt(n_pairs)
# 0.6927234

result_grande=wilcox.test(y=cell_sync_ancestors[cell_sync_ancestors$strain=='Grande',]$daughter_minutes, x=cell_sync_ancestors[cell_sync_ancestors$strain=='Grande',]$mother_minutes, paired=TRUE)
result_grande
# V = 223.5, p-value = 0.006635
median(cell_sync_ancestors[cell_sync_ancestors$strain=='Grande',]$daughter_minutes-cell_sync_ancestors[cell_sync_ancestors$strain=='Grande',]$mother_minutes)
# 5

Z <- qnorm(result_grande$p.value / 2)
n_pairs <- nrow(cell_sync_ancestors[cell_sync_ancestors$strain=='Petite',])
abs(Z) / sqrt(n_pairs)
# 0.4002517

# Filter for grande and petite strains
strains_of_interest <- c("Grande", "Petite")

test_slope_equals_1 <- function(data) {
  
  model <- lm(daughter_minutes ~ mother_minutes -1 , data = data)
  
  # --- Intercept ---
  # intercept_est <- coef(model)["(Intercept)"]
  # intercept_se  <- summary(model)$coefficients["(Intercept)", "Std. Error"]
  # intercept_ci  <- confint(model)["(Intercept)", ]
  df            <- df.residual(model)
  
  # Test: intercept == 0
  # t_stat_int  <- (intercept_est - 0) / intercept_se
  # p_value_int <- 2 * pt(abs(t_stat_int), df = df, lower.tail = FALSE)
  
  # --- Slope ---
  slope_est <- coef(model)["mother_minutes"]
  slope_se  <- summary(model)$coefficients["mother_minutes", "Std. Error"]
  slope_ci  <- confint(model)["mother_minutes", ]
  
  # Test: slope == 1
  t_stat_slope  <- (slope_est - 1) / slope_se
  p_value_slope <- 2 * pt(abs(t_stat_slope), df = df, lower.tail = FALSE)
  
  tibble(
    # intercept        = round(intercept_est, 4),
    # intercept_ci_lower = round(intercept_ci["2.5 %"], 4),
    # intercept_ci_upper = round(intercept_ci["97.5 %"], 4),
    # intercept_se     = round(intercept_se, 4),
    # intercept_t_stat = round(t_stat_int, 4),
    # intercept_p      = p_value_int,
    
    slope            = round(slope_est, 4),
    slope_ci_lower   = round(slope_ci["2.5 %"], 4),
    slope_ci_upper   = round(slope_ci["97.5 %"], 4),
    slope_se         = round(slope_se, 4),
    slope_t_stat     = round(t_stat_slope, 4),
    slope_p          = p_value_slope,
    
    df               = df
  )
}

# Run test for each strain
results <- cell_sync_ancestors %>%
  filter(strain %in% strains_of_interest) %>%
  group_by(strain) %>%
  group_modify(~ test_slope_equals_1(.x)) %>%
  ungroup()

# print(results |> select(strain, contains("slope")))
#   strain slope slope_ci_lower slope_ci_upper slope_se slope_t_stat slope_p
#   <chr>  <dbl>          <dbl>          <dbl>    <dbl>        <dbl>   <dbl>
# 1 Grande 0.940          0.732          1.15     0.104       -0.580  0.564 
# 2 Petite 0.636          0.309          0.964    0.162       -2.24   0.0303


# print(results |> select(strain, contains("intercept")))
#   strain intercept intercept_ci_lower intercept_ci_upper intercept_se intercept_t_stat intercept_p
# 1 Grande      10.7              -11.6               33.0         11.1            0.966    0.339   
# 2 Petite      81.3               38.8              123.8          21.1            3.86     0.000369


petite_dt=doubled[doubled$strain=='petite',]
petite_dt$division_number_12=ifelse(as.numeric(petite_dt$division_number)>=2, "2+", "1")
petite_dt$name='Petite'

grande_dt=doubled[doubled$strain=='grande',]
grande_dt$division_number_12=ifelse(as.numeric(grande_dt$division_number)>=2, "2+", "1")
grande_dt$name='Grande'

p4_ancestors_dt=rbind(petite_dt, grande_dt)
p4_ancestors_dt$division_number_12=factor(p4_ancestors_dt$division_number_12, levels=c("1", "2+"))
table(p4_ancestors_dt$strain)
summary(p4_ancestors_dt)


p4_doubling_times_violin_ancestors=ggplot(p4_ancestors_dt, 
                                          aes(x=division_number_12, y=hours, fill=timepoint, col=name))+
  facet_wrap(~name, ncol=5)+
  # geom_boxplot()+
  geom_jitter(alpha=0.75)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  xlab('Division Number')+
  ylab('Cell Doubling Time\n(hours)')+
  guides(fill='none', col='none')+
  # petite_t200_colors+
  scale_color_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  # theme_classic(base_size = 11)+
  theme_classic(base_size=10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times_violin_ancestors

p4_ancestors_dt %>%
  group_by(name, division_number_12) %>%
  summarize(mean_min=median(hours*60),
            sd_min=sd(hours*60)/(mean_min/106))
# Formula to get the standard deviation parameter for simulations:
# sd_dist = delay_j/default_mean*std_param
# we want so it is: std_param = sd_dist / (delay_j/default_mean)
#   sd_dist: is the measured std from the distribution (empirical data)
#   delay_j: measured mean of the distribution
#   default_mean: mean of the smallest distribution (base mean)
#   name   division_number mean_min sd_min
# 1 Grande 1                   107.   15.6
# 2 Grande 2                   106.   14.7
# 3 Petite 1                   153.   28.0
# 4 Petite 2                   123.   28.8
# Grande: mean(c(15.6, 14.7)) = 15.15
# Petite: mean(c(28, 28.8)) = 28.4

p4_ancestors_dt[p4_ancestors_dt$name=='Petite',] %>%
  group_by(name, division_number) %>%
  summarize(mean_min=mean(hours*60),
            sd_min=sd(hours*60)/(mean_min/123))

table(p4_ancestors_dt$strain, p4_ancestors_dt$division_number_12)
  #          1  2+
  # grande 147  56
  # petite 135  56


p4_ancestors_dt %>%
  group_by(name, division_number_12) %>%
  summarize(median_minutes=median(hours*60))
#   name   division_number_12 median_minutes
# 1 Grande 1                            105 
# 2 Grande 2+                           105 
# 3 Petite 1                            150 
# 4 Petite 2+                           118.

wilcox.test(p4_ancestors_dt[p4_ancestors_dt$name=='Petite',]$hours ~ p4_ancestors_dt[p4_ancestors_dt$name=='Petite',]$division_number_12)
# p-value = 3.771e-07

cliff.delta(p4_ancestors_dt[p4_ancestors_dt$name=='Petite',]$hours ~ p4_ancestors_dt[p4_ancestors_dt$name=='Petite',]$division_number_12)
# delta estimate: 0.4670635 (medium)

wilcox.test(p4_ancestors_dt[p4_ancestors_dt$name=='Grande',]$hours ~ p4_ancestors_dt[p4_ancestors_dt$name=='Grande',]$division_number_12)
# p-value = 0.7969

cliff.delta(p4_ancestors_dt[p4_ancestors_dt$name=='Grande',]$hours ~ p4_ancestors_dt[p4_ancestors_dt$name=='Grande',]$division_number_12)
# delta estimate: 0.02332362 (negligible)

# Create the first column with two plots
column1 <- plot_grid(img_plot_grande_v2, img_plot_petite_v2, 
                     labels = c('A', 'B'), 
                     ncol = 1, 
                     label_size = 16, 
                     rel_heights = c(0.8, 1))

column2=plot_grid(p3_ancestors, p4_doubling_times_violin_ancestors, 
                  labels=c('D', 'E'),
                  ncol=1,
                  label_size = 16,
                  rel_heights = c(1,1))


row_v2=plot_grid(p2_ancestors_v2, column2,
                 labels = c('C', ''),
                 ncol = 2,
                 label_size = 16,
                 rel_widths = c(1.4, 1.4))

figure_1_paper_ancestors_dt_v2 <- plot_grid(column1, row_v2, 
                                         ncol = 1, 
                                         rel_heights = c(0.7, 1))
figure_1_paper_ancestors_dt_v2


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_synchronization_ancestors_dt_26sep2025.png',
#        plot=figure_1_paper_ancestors_dt_v2, dpi='retina', width=10, height=12, bg='white')



#### version 2 - 9dec2025

# Doubling time standard deviation
# Grande: 15.15
# Petite: 28.4
# Doubling time delay (% second division):
# Grande: 0.009433962 (0.9 %)
# Petite: 0.2439024 (24.3 %)
# Grande color, "#AE93BEFF"
# Petite color, "#F0D77BFF"


contour_sd=ggplot(df_grid_analysis, aes(x = delay*100/60, y = std, z = combined_sd)) +
  geom_contour_filled(breaks=seq(0, 45, 5)) +  # You can adjust bins or use breaks = seq(min, max, by = step)
  # scale_fill_viridis_d(name = "Combined\nStandard\nDeviation") +  # Using viridis, but you can change to scale_fill_brewer()
  scale_fill_brewer(palette = "Purples", name = "Asynchrony\n(Combined\nStandard\nDeviation)")+
  geom_point(data = data.frame(x = 0.9, y = 15.15), 
             aes(x = x, y = y, z = NULL), 
             shape = 8, size = 1, color = "#AE93BEFF", stroke = 1) +  # Grande
  geom_point(data = data.frame(x = 24.3, y = 28.4), 
             aes(x = x, y = y, z = NULL), 
             shape = 8, size = 1, color = "#F0D77BFF", stroke = 1) +  # Petite
  labs(x = "Delay (% Second Division)", 
       y = "Doubling Time\nStandard Deviation") +
  theme_classic(base_size=10) +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL
contour_sd


# Create the first column with two plots
column1 <- plot_grid(img_plot_grande_v2, img_plot_petite_v2, 
                     labels = c('A', 'B'), 
                     ncol = 1, 
                     label_size = 11, 
                     rel_heights = c(1, 0.8))

column2=plot_grid(p3_ancestors, contour_sd,
                  labels=c('D', 'F'),
                  ncol=1,
                  label_size = 11,
                  rel_heights = c(1,1))

column3=plot_grid(p2_ancestors_v2, p4_doubling_times_violin_ancestors,
                  labels = c('C', 'E'),
                  ncol = 1,
                  label_size = 11,
                  rel_widths = c(1.4, 1.4))


row1=plot_grid(column3, column2,
                 labels = c('C', ''),
                 ncol = 2,
                 label_size = 11,
                 rel_widths = c(1.4, 1.4))


figure_2_paper_ancestors_dt_v3 <- plot_grid(column1, row1, 
                                            ncol = 1, 
                                            rel_heights = c(0.7, 1))
figure_2_paper_ancestors_dt_v3


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_synchronization_ancestors_dt_23feb2026.png',
#        plot=figure_2_paper_ancestors_dt_v3, dpi='retina', width=6, height=7.2, bg='white')


# Updated version April 1st 2026

p4_doubling_times_violin_ancestors_v2=ggplot(p4_ancestors_dt, 
                                          aes(x=name, y=hours, fill=timepoint, col=name))+
  # facet_wrap(~name, ncol=5)+
  # geom_boxplot()+
  geom_jitter(alpha=0.75)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='median', geom='crossbar', col='black')+
  xlab('Strain')+
  ylab('Cell Doubling Time\n(hours)')+
  guides(fill='none', col='none')+
  # petite_t200_colors+
  scale_color_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  # theme_classic(base_size = 11)+
  theme_classic(base_size=10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times_violin_ancestors_v2

fligner.test(hours ~ name, data = p4_ancestors_dt)
# Fligner-Killeen:med chi-squared = 99.527, df = 1, p-value < 2.2e-16

table(p4_ancestors_dt$name)


# Create the first column with two plots
column1 <- plot_grid(img_plot_grande_v2, img_plot_petite_v2, 
                     labels = c('A', 'B'), 
                     ncol = 1, 
                     label_size = 11, 
                     rel_heights = c(1, 0.8))

column2=plot_grid(p4_doubling_times_violin_ancestors_v2, contour_sd,
                  labels=c('D', 'F'),
                  ncol=1,
                  label_size = 11,
                  rel_heights = c(1,1))

column3=plot_grid(p2_ancestors_v2, p3_ancestors,
                  labels = c('C', 'E'),
                  ncol = 1,
                  label_size = 11,
                  rel_widths = c(1.4, 1.4))


row1=plot_grid(column3, column2,
               labels = c('C', ''),
               ncol = 2,
               label_size = 11,
               rel_widths = c(1.4, 1.4))


figure_2_paper_ancestors_dt_v4 <- plot_grid(column1, row1, 
                                            ncol = 1, 
                                            rel_heights = c(0.7, 1))
figure_2_paper_ancestors_dt_v4


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_synchronization_ancestors_dt_1apr2026.png',
#        plot=figure_2_paper_ancestors_dt_v4, dpi='retina', width=6, height=7.2, bg='white')



#### Version 28Apr2026 ####



p4_doubling_times_ancestors_v3=ggplot(p4_ancestors_dt, 
                                             aes(x=name, y=hours, fill=name, col=name))+
  # facet_wrap(~name, ncol=5)+
  geom_boxplot()+
  # geom_jitter(alpha=0.75)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='median', geom='crossbar', col='black')+
  xlab('Strain')+
  ylab('Cell Doubling Time\n(hours)')+
  guides(fill='none', col='none')+
  # petite_t200_colors+
  scale_fill_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  scale_color_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  theme_classic(base_size=10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) + # , axis.text.x = element_blank()
  NULL
p4_doubling_times_ancestors_v3



column1 <- plot_grid(img_plot_grande_v2, img_plot_petite_v2, 
                     labels = c('A', 'B'), 
                     ncol = 1, 
                     label_size = 10, 
                     rel_heights = c(1, 0.8))

row_cde=plot_grid(p2_ancestors_v2, p4_doubling_times_ancestors_v3, p3_ancestors,
                  labels=c('', 'D', 'E'),
                  nrow=1,
                  label_size = 10,
                  rel_widths = c(0.9, 0.4, 0.9))



row_fg=plot_grid(p4_doubling_times_violin_ancestors, contour_sd,
                  labels = c('F', 'G'),
                  nrow = 1,
                  label_size = 10,
                 rel_widths = c(1, 1))


row_cde_fg=plot_grid(row_cde, row_fg,
               ncol = 1,
               labels=c('C',''),
               label_size = 10,
               rel_heights = c(1.4, 1.4))


figure_2_paper_ancestors_dt_v5 <- plot_grid(column1, row_cde_fg, 
                                            ncol = 1, 
                                            rel_heights = c(0.7, 1))
figure_2_paper_ancestors_dt_v5


ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_synchronization_ancestors_dt_29apr2026.png',
       plot=figure_2_paper_ancestors_dt_v5, dpi='retina', width=6, height=7.2, bg='white')


# Figure 7 ####

p3_paper=ggplot(cell_sync[cell_sync$strain!='grande',], 
                aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.5)+
  geom_abline(linetype='dashed')+
  facet_wrap(~timepoint, ncol=5)+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter Doubling\nTime (minutes)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p3_paper

p4_doubling_times_violin=ggplot(doubled[as.numeric(doubled$division_number)<=2 & doubled$strain!='grande',], 
                                aes(x=division_number, y=hours, fill=timepoint, col=timepoint))+
  facet_wrap(~timepoint, ncol=5)+
  # geom_boxplot()+
  geom_jitter(alpha=0.2)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  xlab('Division Number')+
  ylab('Cell Doubling Time\n(hours)')+
  guides(fill='none', col='none')+
  petite_t200_colors+
  scale_fill_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times_violin


figure_evol_sync=plot_grid(p3_paper, p4_doubling_times_violin, 
                           labels=c('A', 'B'), ncol=1, label_size=16)
figure_evol_sync

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_evolution_synchrony_PA_line_5may2025.png',
#        plot=figure_evol_sync, dpi='retina', width=10, height=6, bg='white')



doubled$division_number_12=ifelse(as.numeric(doubled$division_number)>=2, "2+", "1")
doubled$division_number_12=factor(doubled$division_number_12, levels=c("1", "2+"))

# not using timelapse 2022nov23_gob385_2 as it has higher doubling times that differ from the 2 other time lapses
# of PA4 t200
p_all_timepoints=ggplot(doubled[doubled$strain!='grande' & doubled$id_file!='2022nov23_gob385_2',], 
                        aes(x=strain, y=hours, fill=division_number_12, col=division_number_12))+
  facet_wrap(~timepoint, ncol=5, scales="free_x") +
  force_panelsizes(cols = c(1, 5, 5, 5, 5)) +
  geom_boxplot(outlier.size=0.25)+
  # geom_jitter(alpha=0.2)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='median', geom='crossbar', col='black', position='dodge')+
  xlab('Replicate Population')+
  ylab('Cell Doubling Time (hours)')+
  guides(fill='none', col='none')+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  scale_color_manual(values=c("lightgreen","darkgreen"))+
  # scale_fill_manual(values=c("#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # scale_color_manual(values=c("#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  theme_classic(base_size = 10)+
  scale_x_discrete(labels = function(x) ifelse(grepl("^PA", x), gsub("PA", "", x), ""))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p_all_timepoints


doubled_no_pa4_t200_outlier=doubled[doubled$id_file!='2022nov23_gob385_2',]
for (temp_timepoint in unique(doubled_no_pa4_t200_outlier$timepoint)){
  for (temp_strain in unique(doubled_no_pa4_t200_outlier[doubled_no_pa4_t200_outlier$timepoint==temp_timepoint,]$strain)){
    temp_strain_timepoint=paste(temp_strain, temp_timepoint, sep="_")
    temp_df=doubled_no_pa4_t200_outlier[doubled_no_pa4_t200_outlier$strain_timepoint==temp_strain_timepoint,]
    
    temp_median_diff=median(temp_df[temp_df$division_number_12=='1',]$minutes)-median(temp_df[temp_df$division_number_12=='2+',]$minutes)
    
    temp_p_value=wilcox.test(temp_df$hours~temp_df$division_number_12)$p.value
    
    print(paste(temp_strain_timepoint,": median=", temp_median_diff, " p=",  temp_p_value, "  ", temp_p_value<(0.05/21), sep=""))
  }
}
# [1] "PA3_t400: median=0 p=0.00678514968297457  FALSE"
# [1] "PA4_t400: median=0 p=0.61750863593883  FALSE"
# [1] "PA5_t400: median=0 p=0.951532969843358  FALSE"
# [1] "PA1_t400: median=5 p=0.676230215555684  FALSE"
# [1] "PA2_t400: median=10 p=0.00362611207968799  FALSE"
# [1] "PA1_t1000: median=0 p=0.741526421293436  FALSE"
# [1] "PA4_t1000: median=-5 p=0.22018384973057  FALSE"
# [1] "PA3_t1000: median=5 p=0.0698526294296691  FALSE"
# [1] "PA5_t1000: median=0 p=0.172722821193342  FALSE"
# [1] "PA2_t1000: median=5 p=0.336121861072508  FALSE"
# [1] "PA3_t600: median=0 p=0.039877749135946  FALSE"
# [1] "PA4_t600: median=0 p=0.1832091820291  FALSE"
# [1] "PA5_t600: median=0 p=0.569423136689105  FALSE"
# [1] "PA2_t600: median=0 p=0.125565239983413  FALSE"
# [1] "PA1_t600: median=0 p=0.689629388955025  FALSE"
# [1] "PA4_t200: median=10 p=0.041127523553486  FALSE"
# [1] "PA5_t200: median=5 p=0.0928361420904245  FALSE"
# [1] "PA3_t200: median=0 p=0.519740510434724  FALSE"
# [1] "PA1_t200: median=0 p=0.843582984142905  FALSE"
# [1] "PA2_t200: median=0 p=0.503230809022216  FALSE"
# [1] "petite_t0: median=32.5 p=3.77108895395933e-07  TRUE"
# [1] "grande_t0: median=0 p=0.796906120057638  FALSE"


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_evolution_synchrony_PA_line_all_23feb2026.png',
#        plot=p_all_timepoints, dpi='retina', width=6.5, height=3.5, bg='white')

# not using timelapse 2022nov23_gob385_2 as it has higher doubling times that differ from the 2 other time lapses
# of PA4 t200
p_all_timepoints_no_outliers=ggplot(doubled[doubled$strain!='grande' & doubled$id_file!='2022nov23_gob385_2',], 
                        aes(x=strain, y=hours, fill=division_number_12, col=division_number_12))+
  facet_wrap(~timepoint, ncol=5, scales="free_x") +
  force_panelsizes(cols = c(1, 5, 5, 5, 5)) +
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun='median', geom='crossbar', col='black', position='dodge')+
  xlab('Replicate Population')+
  ylab('Cell Doubling Time (hours)')+
  guides(fill='none', col='none')+
  scale_fill_manual(values=c("lightgreen","darkgreen"))+
  scale_color_manual(values=c("lightgreen","darkgreen"))+
  # scale_fill_manual(values=c("#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # scale_color_manual(values=c("#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  theme_classic(base_size = 10)+
  scale_y_continuous(limits=c(1, 4.5))+
  scale_x_discrete(labels = function(x) ifelse(grepl("^PA", x), gsub("PA", "", x), ""))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p_all_timepoints_no_outliers


cell_sync_no_grande=cell_sync[cell_sync$strain!='grande',]
temp_cell_sync=data.frame()
for(i in seq(nrow(cell_sync_no_grande))){
  temp_cell_sync_row=cell_sync_no_grande[i,]
  temp_strain=temp_cell_sync_row$strain
  temp_timepoint=temp_cell_sync_row$timepoint
  
  temp_hour_daughter=temp_cell_sync_row$daughter_minutes
  temp_division_number='1'
  temp_row=data.frame(strain=temp_strain, timepoint=temp_timepoint, 
                      division_number_12=temp_division_number, hours=temp_hour_daughter/60)
  temp_cell_sync=rbind(temp_cell_sync, temp_row)
  
  temp_hour_mother=temp_cell_sync_row$mother_minutes
  temp_division_number='2'
  temp_row=data.frame(strain=temp_strain, timepoint=temp_timepoint, 
                      division_number_12=temp_division_number, hours=temp_hour_mother/60)
  temp_cell_sync=rbind(temp_cell_sync, temp_row)
}

ggplot(cell_sync[cell_sync$timepoint=='t200' & cell_sync$strain=='PA4',], aes(x=daughter_minutes-mother_minutes))+
  geom_histogram()+
  NULL

ggplot(cell_sync[cell_sync$timepoint!='t0',], aes(x=daughter_minutes-mother_minutes))+
  geom_histogram(binwidth=5)+
  facet_grid(timepoint~strain)+
  geom_vline(xintercept=0, linetype='dashed')+
  NULL

ggplot(cell_sync_ancestors, 
       aes(x=mother_minutes, y=daughter_minutes-mother_minutes))+
  geom_point()+
  # facet_wrap(~strain)+
  facet_grid(strain~date)+
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, linewidth = 0.7)+
  NULL

ggplot(cell_sync[cell_sync$timepoint!='t0',], 
       aes(x=mother_minutes, y=daughter_minutes-mother_minutes))+
  geom_point()+
  # facet_wrap(~strain)+
  facet_grid(strain~timepoint)+
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, linewidth = 0.7)+
  NULL



# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_evolution_synchrony_PA_line_all_no_outliers_27apr2026.png',
#        plot=p_all_timepoints_no_outliers, dpi='retina', width=6.5, height=3, bg='white')

# View(doubled %>%
#   group_by(strain_timepoint) %>%
#   summarise(mean_doubling_time=mean(minutes)))

# Supplementary Figure ####

p_mother_daughter=ggplot(cell_sync[cell_sync$timepoint!='t0' & cell_sync$id_file!='2022nov23_gob385_2',], 
                            aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.5, size=1)+
  geom_abline(linetype='dashed')+
  # geom_smooth(method = 'lm', formula = y ~ x - 1, se = TRUE, linewidth = 0.7) +
  facet_grid2(strain ~ timepoint, axes = "all", remove_labels = "all")+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter DoublingTime (minutes)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  scale_colour_ghibli_d("PonyoMedium", direction = -1)+
  theme_classic(base_size = 10)+
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  )+
  NULL
p_mother_daughter


cell_sync_no_pa4_t200_outlier=cell_sync[cell_sync$id_file!='2022nov23_gob385_2',]
for (temp_timepoint in unique(cell_sync_no_pa4_t200_outlier$timepoint)){
  for (temp_strain in unique(cell_sync_no_pa4_t200_outlier[cell_sync_no_pa4_t200_outlier$timepoint==temp_timepoint,]$strain)){
    temp_strain_timepoint=paste(temp_strain, temp_timepoint, sep="_")
    temp_df=cell_sync_no_pa4_t200_outlier[cell_sync_no_pa4_t200_outlier$strain_timepoint==temp_strain_timepoint,]
    
    temp_samples=nrow(temp_df)
    temp_mean=median(temp_df$daughter_minutes-temp_df$mother_minutes)
    
    temp_p_value=wilcox.test(y=temp_df$daughter_minutes, x=temp_df$mother_minutes, paired=TRUE)$p.value
    
    print(paste(temp_strain_timepoint,": median=",temp_mean, " n=", temp_samples, " p-value= ", temp_p_value, "  ", temp_p_value<(0.05/21), sep=""))
  }
}

# [1] "PA1_t200: median=5 p-value= 0.00929154518699802  FALSE"
# [1] "PA2_t200: median=0 p-value= 0.0918346064781602  FALSE"
# [1] "PA3_t200: median=0 p-value= 0.505792044277484  FALSE"
# [1] "PA4_t200: median=5 p-value= 0.00398595954631529  FALSE"
# [1] "PA5_t200: median=0 p-value= 0.0381508024479592  FALSE"

# [1] "PA1_t400: median=0 p-value= 0.00076705865400928  TRUE" ###
# [1] "PA2_t400: median=5 p-value= 1.01481497948934e-07  TRUE" ###
# [1] "PA3_t400: median=0 p-value= 0.0320020646359707  FALSE"
# [1] "PA4_t400: median=0 p-value= 0.0353336996104896  FALSE"
# [1] "PA5_t400: median=0 p-value= 0.266157991908009  FALSE"

# [1] "PA1_t600: median=0 p-value= 0.0169693577434589  FALSE"
# [1] "PA2_t600: median=0 p-value= 1.39115784345679e-05  TRUE" ###
# [1] "PA3_t600: median=0 p-value= 1.27031909192806e-06  TRUE" ###
# [1] "PA4_t600: median=0 p-value= 0.0650362919807907  FALSE"
# [1] "PA5_t600: median=0 p-value= 0.00025804651334413  TRUE" ###


# [1] "PA1_t1000: median=0 p-value= 0.000447569467883027  TRUE" ###
# [1] "PA2_t1000: median=0 p-value= 0.570581277021137  FALSE"
# [1] "PA3_t1000: median=0 p-value= 0.0276433269884977  FALSE"
# [1] "PA4_t1000: median=0 p-value= 0.420573936928138  FALSE"
# [1] "PA5_t1000: median=0 p-value= 0.471935910252018  FALSE"

# [1] "petite_t0: median=35 p-value= 2.62363453085217e-06  TRUE"
# [1] "grande_t0: median=5 p-value= 0.00663480591712819  FALSE"


ggplot(cell_sync_no_pa4_t200_outlier[cell_sync_no_pa4_t200_outlier$strain_timepoint=='petite_t0',],
       aes(x=diff_minutes))+
  geom_histogram()+
  NULL

for(temp_strain_timepoint in sort(unique(cell_sync_no_pa4_t200_outlier$strain_timepoint))){
  if (temp_strain_timepoint == 'petite_t0'){
    next
  }
  petite_cell_sync_df=cell_sync_no_pa4_t200_outlier[cell_sync_no_pa4_t200_outlier$strain_timepoint=='petite_t0',]
  temp_cell_sync_df=cell_sync_no_pa4_t200_outlier[cell_sync_no_pa4_t200_outlier$strain_timepoint==temp_strain_timepoint,]
  
  temp_result=wilcox.test(temp_cell_sync_df$diff_minutes, petite_cell_sync_df$diff_minutes, alternative = "less")$p.value
  print(paste(temp_strain_timepoint, ": p=", temp_result, sep=""))
}
# [1] "grande_t0: p=5.19013084054496e-13"
# [1] "PA1_t1000: p=1.1436428987509e-16"
# [1] "PA1_t200: p=5.19049193319374e-09"
# [1] "PA1_t400: p=1.50099697147939e-22"
# [1] "PA1_t600: p=2.03115264148666e-31"
# [1] "PA2_t1000: p=3.14151240466191e-24"
# [1] "PA2_t200: p=4.95636358879076e-19"
# [1] "PA2_t400: p=1.17047478773858e-13"
# [1] "PA2_t600: p=3.2179058817524e-28"
# [1] "PA3_t1000: p=5.85563586229025e-17"
# [1] "PA3_t200: p=3.18438709250436e-16"
# [1] "PA3_t400: p=1.94486086214801e-26"
# [1] "PA3_t600: p=2.04678584914912e-25"
# [1] "PA4_t1000: p=2.06259001024433e-21"
# [1] "PA4_t200: p=1.54845713405484e-11"
# [1] "PA4_t400: p=3.06076431446822e-26"
# [1] "PA4_t600: p=1.91248353052381e-24"
# [1] "PA5_t1000: p=1.50240382118024e-12"
# [1] "PA5_t200: p=1.34808887359538e-17"
# [1] "PA5_t400: p=1.95876347455339e-24"
# [1] "PA5_t600: p=1.29221208563382e-30"


comparison_differences=data.frame()

for (temp_timepoint in unique(cell_sync$timepoint)){
  for (temp_strain in unique(cell_sync[cell_sync$timepoint==temp_timepoint,]$strain)){
    temp_strain_timepoint=paste(temp_strain, temp_timepoint, sep="_")
    temp_df=cell_sync[cell_sync$strain_timepoint==temp_strain_timepoint,]
    
    temp_mean_cell_sync=mean(temp_df$daughter_minutes-temp_df$mother_minutes)
    temp_median_cell_sync=median(temp_df$daughter_minutes-temp_df$mother_minutes)
    
    temp_doubled_df=doubled[doubled$strain_timepoint==temp_strain_timepoint,]
    
    temp_mean_doubled=mean(temp_doubled_df[temp_doubled_df$division_number_12=='1',]$minutes)-mean(temp_doubled_df[temp_doubled_df$division_number_12=='2+',]$minutes)
    temp_median_doubled=median(temp_doubled_df[temp_doubled_df$division_number_12=='1',]$minutes)-median(temp_doubled_df[temp_doubled_df$division_number_12=='2+',]$minutes)
    
    temp_row=data.frame(strain_timepoint=temp_strain_timepoint, mean_all=temp_mean_doubled, 
                        median_all=temp_median_doubled, mean_pairs=temp_mean_cell_sync, median_pairs=temp_median_cell_sync)
    comparison_differences=rbind(comparison_differences, temp_row)
  }
}


ggplot()+
  geom_point(data=comparison_differences, aes(x=strain_timepoint, mean_all), shape=1, col='lightblue')+
  geom_point(data=comparison_differences, aes(x=strain_timepoint, mean_pairs), shape=2, col='blue')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  NULL

# Run test for each strain
results <- cell_sync[cell_sync$timepoint!='t0',] %>%
  group_by(strain, timepoint) %>%
  group_modify(~ test_slope_equals_1(.x)) %>%
  ungroup()


# results |> filter(slope_p <= 0.05/40) |> select(strain, timepoint, contains("slope")) |> print()
#   strain timepoint slope slope_ci_lower slope_ci_upper slope_se slope_t_stat  slope_p
# 1 PA1    t400      0.632          0.514          0.750   0.0598        -6.16 4.93e- 9
# 2 PA1    t1000     0.640          0.507          0.774   0.0675        -5.33 4.16e- 7
# 3 PA2    t600      0.855          0.796          0.913   0.0296        -4.91 1.52e- 6
# 4 PA2    t1000     0.896          0.836          0.956   0.0303        -3.42 8.15e- 4
# 5 PA3    t400      0.837          0.769          0.906   0.0348        -4.68 5.48e- 6
# 6 PA3    t1000     0.325          0.209          0.441   0.0584       -11.6  3.36e-19
# 7 PA4    t600      0.680          0.618          0.742   0.0314       -10.2  1.69e-20
# 8 PA5    t400      0.430          0.300          0.559   0.0656        -8.69 2.39e-15


# results |> filter(slope_p <= 0.05/40) |> select(strain, timepoint, contains("intercept")) |> print()
# strain timepoint intercept intercept_ci_lower intercept_ci_upper intercept_se intercept_t_stat intercept_p
# 1 PA1    t400           38.4              26.6                50.3         6.01             6.40    1.45e- 9
# 2 PA1    t1000          44.6              28.6                60.7         8.11             5.50    1.94e- 7
# 3 PA2    t600           14.3               8.97               19.6         2.70             5.29    2.35e- 7
# 4 PA2    t1000          10.1               4.38               15.9         2.91             3.48    6.66e- 4
# 5 PA3    t400           15.0               8.34               21.6         3.36             4.46    1.41e- 5
# 6 PA3    t1000          66.6              54.8                78.5         5.95            11.2     1.85e-18
# 7 PA4    t600           32.3              25.6                38.9         3.36             9.61    1.04e-18
# 8 PA5    t400           54.9              42.3                67.5         6.40             8.58    4.71e-15



ggplot(cell_sync[cell_sync$timepoint=='t200' & cell_sync$strain=='PA5',], 
       aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.5, size=1)+
  geom_abline(linetype='dashed')+
  # geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, linewidth = 0.7) +
  geom_smooth(method = 'lm', formula = y ~ x - 1, se = TRUE, linewidth = 0.7) +
  facet_grid2(strain ~ timepoint, axes = "all", remove_labels = "all")+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter DoublingTime (minutes)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  xlim(c(0, 300))+
  ylim(c(0,300))+
  scale_colour_ghibli_d("PonyoMedium", direction = -1)+
  theme_classic(base_size = 10)+
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  )+
  NULL

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig_mother_daughter_synchrony_PA_line_all_27apr2026.png',
#               plot=p_mother_daughter, dpi='retina', width=6.5, height=5.5, bg='white')

