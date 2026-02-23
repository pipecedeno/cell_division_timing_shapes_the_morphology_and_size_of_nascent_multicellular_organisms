
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
first_div_percentages=read.csv('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/percentage_cells_in_first_division_23sep2025.csv', header=TRUE)


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
  
  first_div_samples=(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100)*(n_samples*2)
  second_div_samples=(1-(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100))*(n_samples*2)
  
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
    
    # Obtain percentage of cells in the first cell division and the rest
    first_div_samples=(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100)*(n_samples*2)
    second_div_samples=(1-(first_div_percentages[first_div_percentages$delay_time==delay_j,]$fraction_first_cycle/100))*(n_samples*2)
    
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
ancestor_tl_dt$strain=factor(ancestor_tl_dt$strain, levels=c("petite", "grande"))


p2_ancestors_v2 = ggplot(p-value = 0.2711, 
                         aes(x = 5*FRAME, fill = strain)) +
  xlim(c(-15, 350))+
  geom_histogram(binwidth = 10, position = 'identity')+
  facet_wrap(~strain, ncol=1, scales='free')+
  xlab('Time from First Division (min)')+
  ylab('Number of Divisions')+
  labs(fill = "Strain")+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#F0D77BFF", "#AE93BEFF"))+
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


petite_time <- ancestor_tl_dt[ancestor_tl_dt$strain == 'petite' & ancestor_tl_dt$FRAME!=0, ]$FRAME * 5
ks.test(petite_time, "punif", min = min(petite_time), max = max(petite_time))
# p-value = 0.000328

grande_time <- ancestor_tl_dt[ancestor_tl_dt$strain == 'grande' & ancestor_tl_dt$FRAME!=0, ]$FRAME * 5
ks.test(grande_time, "punif", min = min(grande_time), max = max(grande_time))
# p-value = 8.32e-08

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
                    aes(x=mother_minutes, y=daughter_minutes, col=strain))+
  geom_point(alpha=0.75)+
  geom_abline(linetype='dashed')+
  facet_wrap(~strain, ncol=5)+
  guides(col='none')+
  xlab('Mother Doubling Time (min)')+
  ylab('Daughter Doubling\nTime (min)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  scale_color_manual(values=c("#AE93BEFF", "#F0D77BFF"))+
  theme_classic(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p3_ancestors

petite_dt=doubled[doubled$strain=='petite' & as.numeric(doubled$division_number)<=2,]
petite_dt$name='Petite'

grande_dt=doubled[doubled$strain=='grande' & as.numeric(doubled$division_number)<=2,]
grande_dt$name='Grande'

p4_ancestors_dt=rbind(petite_dt, grande_dt)
table(p4_ancestors_dt$strain)
summary(p4_ancestors_dt)


p4_doubling_times_violin_ancestors=ggplot(p4_ancestors_dt, 
                                          aes(x=division_number, y=hours, fill=timepoint, col=name))+
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
  group_by(name, division_number) %>%
  summarize(mean_min=mean(hours*60),
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

table(p4_ancestors_dt$strain, p4_ancestors_dt$division_number)
  #          1   2   3   4   5   6   7   8
  # grande 147  53   0   0   0   0   0   0
  # petite 135  48   0   0   0   0   0   0


p4_ancestors_dt %>%
  group_by(name) %>%
  summarize(mean_minutes=mean(hours*60))
#   name   mean_minutes
# 1 Grande         107.57
# 2 Petite         145.30

t.test(p4_ancestors_dt$hours*60~p4_ancestors_dt$name)
# p-value < 2.2e-16

cohen.d(p4_ancestors_dt$hours*60~p4_ancestors_dt$name)
# d estimate: -1.276712 (large)

# Create the first column with two plots
column1 <- plot_grid(img_plot_petite_v2, img_plot_grande_v2, 
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
  scale_fill_brewer(palette = "Purples", name = "Combined\nStandard\nDeviation")+
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
        legend.title = element_text(size = 7)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL
contour_sd


# Create the first column with two plots
column1 <- plot_grid(img_plot_petite_v2, img_plot_grande_v2, 
                     labels = c('A', 'B'), 
                     ncol = 1, 
                     label_size = 11, 
                     rel_heights = c(0.8, 1))

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


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_2_synchronization_ancestors_dt_9dec2025.png',
#        plot=figure_2_paper_ancestors_dt_v3, dpi='retina', width=6, height=7.2, bg='white')


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




