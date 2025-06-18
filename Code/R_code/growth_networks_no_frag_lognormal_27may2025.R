
#Date: 29may2025
# This code is used to analyze the results of growing the networks without fragmentation
# so it is only growing them until 1200 nodes in size, but instead of using the empirical data
# we are going to use a smoothed data using a log normal distribution

library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(ggridges)
library(stringi)
library(glue)
library(MASS)
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
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))



#test fitting lognormal distribution ####
petite_dt=read.csv("~/work_dir/observed_synchrony/data/petite_doubling_time_dist_2023may30.csv", header=TRUE)
grande_dt=read.csv("~/work_dir/observed_synchrony/data/grande_doubling_time_dist_2024dec2.csv", header=TRUE)
petite_delay_dt=read.csv("~/work_dir/observed_synchrony/data/petite_only_second_doubling_8feb2024.csv", header=TRUE)

table(petite_dt$division_number)

petite_div1=petite_dt[petite_dt$division_number==0,]

fit=fitdistr(petite_div1$minutes, densfun = "lognormal")
fit

ggplot(petite_div1, aes(x=minutes))+
  geom_histogram(aes(y = ..density..),fill = "lightblue", color = "black")+
  # geom_density(alpha = 0.3, fill = "lightgreen")+
  stat_function(fun = dlnorm, args = list(meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"]),
                color = "red")+
  labs(x = "Minutes", y = "Density", title = "Histogram with Fitted Log-Normal Distribution")+
  NULL

#simulating data for testing the distribution ####

#fitting data for petite
petite_div1=petite_dt[petite_dt$division_number==0,]
p_fit_div1=fitdistr(petite_div1$minutes, densfun = "lognormal")
p_mu_div1=as.numeric(p_fit_div1$estimate["meanlog"])
p_sigma_div1=as.numeric(p_fit_div1$estimate["sdlog"])

petite_div2=petite_dt[petite_dt$division_number==1,]
p_fit_div2=fitdistr(petite_div2$minutes, densfun = "lognormal")
p_mu_div2=as.numeric(p_fit_div2$estimate["meanlog"])
p_sigma_div2=as.numeric(p_fit_div2$estimate["sdlog"])

#simulating data for each distribution
sim_obs_df=petite_dt[,c('division_number','minutes')]
sim_obs_df$type='observed'
sim_obs_df$division_number=sim_obs_df$division_number+1
sim_obs_df$strain='Petite'
table(sim_obs_df$division_number)

temp_df=data.frame(minutes=rlnorm(100, meanlog = p_mu_div1, sdlog = p_sigma_div1))
temp_df$type='simulated'
temp_df$division_number=1
temp_df$strain='Petite'

sim_obs_df=rbind(sim_obs_df, temp_df)

temp_df=data.frame(minutes=rlnorm(100, meanlog = p_mu_div2, sdlog = p_sigma_div2))
temp_df$type='simulated'
temp_df$division_number=2
temp_df$strain='Petite'

sim_obs_df=rbind(sim_obs_df, temp_df)

#fitting data for grande
temp_df=grande_dt[,c('division_number','minutes')]
temp_df$type='observed'
temp_df$division_number=temp_df$division_number+1
temp_df$strain='Grande'
sim_obs_df=rbind(sim_obs_df, temp_df)

grande_div1=grande_dt[grande_dt$division_number==0,]
g_fit_div1=fitdistr(grande_div1$minutes, densfun = "lognormal")
g_mu_div1=as.numeric(g_fit_div1$estimate["meanlog"])
g_sigma_div1=as.numeric(g_fit_div1$estimate["sdlog"])

grande_div2=grande_dt[grande_dt$division_number==1,]
g_fit_div2=fitdistr(grande_div2$minutes, densfun = "lognormal")
g_mu_div2=as.numeric(g_fit_div2$estimate["meanlog"])
g_sigma_div2=as.numeric(g_fit_div2$estimate["sdlog"])


# Simulating data from grande
temp_df=data.frame(minutes=rlnorm(100, meanlog = g_mu_div1, sdlog = g_sigma_div1))
temp_df$type='simulated'
temp_df$division_number=1
temp_df$strain='Grande'

sim_obs_df=rbind(sim_obs_df, temp_df)

temp_df=data.frame(minutes=rlnorm(100, meanlog = g_mu_div2, sdlog = g_sigma_div2))
temp_df$type='simulated'
temp_df$division_number=2
temp_df$strain='Grande'

sim_obs_df=rbind(sim_obs_df, temp_df)



#fitting data for Petite without Delay
temp_df=petite_delay_dt[,c('division_number','minutes')]
temp_df$type='observed'
temp_df$division_number=temp_df$division_number+1
temp_df$strain='Petite w/o Delay'
sim_obs_df=rbind(sim_obs_df, temp_df)

pwd_div1=petite_delay_dt[petite_delay_dt$division_number==0,]
pwd_fit_div1=fitdistr(pwd_div1$minutes, densfun = "lognormal")
pwd_mu_div1=as.numeric(pwd_fit_div1$estimate["meanlog"])
pwd_sigma_div1=as.numeric(pwd_fit_div1$estimate["sdlog"])

pwd_div2=petite_delay_dt[petite_delay_dt$division_number==1,]
pwd_fit_div2=fitdistr(pwd_div2$minutes, densfun = "lognormal")
pwd_mu_div2=as.numeric(pwd_fit_div2$estimate["meanlog"])
pwd_sigma_div2=as.numeric(pwd_fit_div2$estimate["sdlog"])


# Simulating data from grande
temp_df=data.frame(minutes=rlnorm(100, meanlog = pwd_mu_div1, sdlog = pwd_sigma_div1))
temp_df$type='simulated'
temp_df$division_number=1
temp_df$strain='Petite w/o Delay'

sim_obs_df=rbind(sim_obs_df, temp_df)

temp_df=data.frame(minutes=rlnorm(100, meanlog = pwd_mu_div2, sdlog = pwd_sigma_div2))
temp_df$type='simulated'
temp_df$division_number=2
temp_df$strain='Petite w/o Delay'

sim_obs_df=rbind(sim_obs_df, temp_df)

table(sim_obs_df$type, sim_obs_df$division_number, sim_obs_df$strain)

sim_obs_df$division_number=factor(sim_obs_df$division_number)

sim_obs_df %>%
  group_by(strain, division_number, type) %>%
  summarise(min=min(minutes),
            mean=mean(minutes),
            median=median(minutes),
            max=max(minutes),
            sd=sd(minutes),
            var=var(minutes))

ggplot(sim_obs_df, aes(x=division_number, y=minutes, fill=type))+
  facet_wrap(~strain)+
  geom_violin()+
  NULL

ggplot(sim_obs_df, aes(x=division_number, y=minutes/60, fill=type))+
  facet_wrap(~strain)+
  geom_violin()+
  ylab("Hours")+
  xlab("Division number")+
  stat_summary(geom='crossbar', fun='mean')+
  NULL

ggplot(sim_obs_df, aes(x=division_number, y=minutes/60, fill=type))+
  facet_wrap(~strain)+
  geom_violin(position = position_dodge(width = 0.8))+
  ylab("Hours")+
  xlab("Division number")+
  stat_summary(geom='crossbar', fun='mean', 
               position = position_dodge(width = 0.8))+
  NULL

# Data frame with fitted parameters ####

# Create the data frame
df_parameters <- data.frame(
  strain = c("Grande", "Petite", "Petite w/o Delay"),
  first_div_mu = c(g_mu_div1, p_mu_div1, pwd_mu_div1),
  first_div_sigma = c(g_sigma_div1, p_sigma_div1, pwd_sigma_div1),
  second_div_mu = c(g_mu_div2, p_mu_div2, pwd_mu_div2),
  second_div_sigma = c(g_sigma_div2, p_sigma_div2, pwd_sigma_div2)
)




setwd("~/Desktop/results_edge_degree_15/fig4_growth_no_frag_every5_continuous")

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




#### Network Diameter

# Raw data
ggplot(sim_df, aes(x=num_nodes, y=diameter, col=strain, fill=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL

summ_diameter <- sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(
    mean_diameter = mean(diameter),
    sd_diameter = sd(diameter),
    n = n(),
    se_diameter = sd_diameter / sqrt(n),
    ci_lower = mean_diameter - qt(0.975, n-1) * se_diameter,
    ci_upper = mean_diameter + qt(0.975, n-1) * se_diameter,
    .groups = 'drop'
  )

ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  NULL




# Undivided Daughters ####

#raw data
ggplot(sim_df, aes(x=num_nodes, y=cases_mother_with_undivided_cells, col=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL



summ_undivided <- sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(
    mean_undivided = mean(cases_mother_with_undivided_cells),
    sd_undivided = sd(cases_mother_with_undivided_cells),
    n = n(),
    se_undivided = sd_undivided / sqrt(n),
    ci_lower = mean_undivided - qt(0.975, n-1) * se_undivided,
    ci_upper = mean_undivided + qt(0.975, n-1) * se_undivided,
    .groups = 'drop'
  )

ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  NULL





# Max Edge Degree ####

ggplot(sim_df, aes(x=num_nodes, y=max_edge_degree, col=strain))+
  geom_point(alpha=0.1)+
  petite_t200_colors+
  facet_wrap(~strain)+
  guides(col='none', fill='none')+
  NULL


summ_max_edge <- sim_df %>%
  group_by(strain, num_nodes) %>%
  summarise(
    mean_max_edge = mean(max_edge_degree),
    sd_max_edge = sd(max_edge_degree),
    n = n(),
    se_max_edge = sd_max_edge / sqrt(n),
    ci_lower = mean_max_edge - qt(0.975, n-1) * se_max_edge,
    ci_upper = mean_max_edge + qt(0.975, n-1) * se_max_edge,
    .groups = 'drop'
  )

ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  NULL





# Paper Figure ####

# Loading sampled times ####

# doub_t_data=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_300n_1200m/sampled_times.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   doub_t_data=rbind(doub_t_data, temp_df)
# }
# doub_t_data$strain <- ifelse(doub_t_data$strain=='petite', 'Petite', doub_t_data$strain)
# doub_t_data$strain <- ifelse(doub_t_data$strain=='petite-second-only', 'Petite w/o Delay', doub_t_data$strain)
# doub_t_data$strain <- ifelse(doub_t_data$strain=='grande', 'Grande', doub_t_data$strain)
# doub_t_data$strain=factor(doub_t_data$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# doub_t_data$number_divisions=doub_t_data$number_divisions+1
# summary(doub_t_data)
# 
# write.csv(doub_t_data, file="sampled_times_27may2025.csv", row.names = FALSE)


doub_t_data=read.csv("sampled_times_27may2025.csv", header = TRUE)
doub_t_data$strain=factor(doub_t_data$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))

summary(doub_t_data)



dt_distributions_v2=ggplot(doub_t_data[doub_t_data$number_divisions<=2,], 
                           aes(x=as.factor(number_divisions), y=minutes/60, fill=strain))+
  facet_wrap(~strain)+
  geom_violin(adjust=2)+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none')+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  NULL
dt_distributions_v2


p_diam=ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Diameter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_diam


p_undivided=ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Mothers with\nDelayed Daughter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_undivided


p_edge_degree=ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Max\nEdge Degree", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_edge_degree

fig=plot_grid(dt_distributions_v2, p_undivided, p_diam, p_edge_degree,
              labels=c('A', 'B', 'C', 'D'), ncol=2, label_size=12)
fig

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_3_growth_no_frag_28may2025_lognormal.png',
       plot=fig, dpi='retina', width=6.5, height=4.5, bg='white')



