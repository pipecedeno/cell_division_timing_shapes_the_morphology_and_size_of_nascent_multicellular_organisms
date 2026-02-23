
# Date: 15Apr2024
# This code is to test what is the effect of carrying capacity in the variation of the
#selection rate of the settling phase in the settling simulations, as we hypothesize that
#increasing the carrying capacity will decrease the selection rate by making the process
#less stochastic

# Modifications:
# Date: 29apr2024
# Added a section where I also check the effect of doing settling simulation starting the
#clusters at the top of the tube instead of at a random position.

# Date: 20May2024
# Added a section for the simulations using bootstrap to sample the clusters again to increase
# the sample size when performing the settling selections, for this simulations only a 
# growth cycle is performed and then the settling with bootstrap is performed.

# Date: 11June2024
# Added the results of the simulations using bootstrap to sample the cluster and also testing 
# different carrying capacities when performing settling selection with the clusters starting at 
# the top

# Date: 2Feb2026 
# Corrected how the volume was being calculated so now the diameter of a clusters is being obtained 
# by calculating the radius of a sphere from the radius of gyration.


library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(cowplot)

theme_set(theme_classic(base_size = 16))

setwd('~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/settling_selection_simulations_30jan2026/')


get_confidence_interval <- function(data) {
  n <- length(data)
  mean_val <- mean(data)
  se <- sd(data) / sqrt(n)
  
  lower_bound <- mean_val - 1.96 * se
  upper_bound <- mean_val + 1.96 * se
  
  return(c(lower_bound, upper_bound))
}

# Different carrying capacities ####


#### Growth phase ####

# growth_prop=read.csv("proportions_growth_registry_all_sim.csv", header=TRUE)

growth_prop=data.frame()

# , "10mill"
for(j in c("10k", "100k", "1mill")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim/proportions_growth_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  growth_prop=rbind(growth_prop, temp_df)
}
write.csv(growth_prop, file="set_sim_growth_registry_30jan2026.csv", row.names=FALSE)

growth_prop=read.csv("set_sim_growth_registry_30jan2026.csv", header = TRUE)
growth_prop$carr_cap=factor(growth_prop$carr_cap, levels=c("10k", "100k", "1mill", "10mill"))
summary(growth_prop)

#Selection rate

growth_prop$cell_sel_r=(log(growth_prop$cells_pop2_a/growth_prop$cells_pop2_b, base=2)-log(growth_prop$cells_pop1_a/growth_prop$cells_pop1_b, base=2))
growth_prop$cell_sel_r_norm=(log(growth_prop$cells_pop2_a/growth_prop$cells_pop2_b, base=2)-log(growth_prop$cells_pop1_a/growth_prop$cells_pop1_b, base=2))/(log(growth_prop$total_cells_a, base=2)-log(growth_prop$total_cells_b, base=2))


ggplot(growth_prop, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

#### Settling fitness ####

settling_prop=data.frame()

# , "10mill"
for(j in c("10k", "100k", "1mill")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim/proportions_settling_selec_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  #filtering out the row if any population went extinct
  temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
  settling_prop=rbind(settling_prop, temp_df)
}
write.csv(settling_prop, file="set_sim_settling_registry_30jan2026.csv", row.names = FALSE)

settling_prop=read.csv("set_sim_settling_registry_30jan2026.csv", header=TRUE)
settling_prop$carr_cap=factor(settling_prop$carr_cap, levels=c("10k", "100k", "1mill", "10mill"))
summary(settling_prop)

#Selection rate

settling_prop$cell_sel_r=(log(settling_prop$cells_pop2_a/settling_prop$cells_pop2_b, base=2)-log(settling_prop$cells_pop1_a/settling_prop$cells_pop1_b, base=2))

ggplot(settling_prop, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL


#### Do selection rates change over time? ####

ggplot(growth_prop,
       aes(x=transfer, y=cell_sel_r_norm, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL

#Getting rid of the first growth selection rate
ggplot(growth_prop[growth_prop$transfer>1,],
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL


ggplot(settling_prop,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Settling Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL


#### Number of generations per transfer ####

growth_prop$num_gen_1=log2(growth_prop$cells_pop1_a)-log2(growth_prop$cells_pop1_b)
growth_prop$num_gen_2=log2(growth_prop$cells_pop2_a)-log2(growth_prop$cells_pop2_b)
growth_prop$num_generations=log2(growth_prop$total_cells_a)-log2(growth_prop$total_cells_b)

mean_num_gen=growth_prop %>%
  group_by(carr_cap, transfer) %>%
  summarise(
    mean_gen_pop_1 = mean(num_gen_1),
    sd_pop1=sd(num_gen_1),
    mean_gen_pop_2 = mean(num_gen_2),
    sd_pop2=sd(num_gen_2),
    mean_gen_all_pop = mean(num_generations),
    sd_pop_all=sd(num_generations)
  )

# this plot is meant for 10mill

plot_carr_cap="1mill"
ggplot()+
  geom_line(data=growth_prop[growth_prop$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=num_gen_1, col=as.factor(sim_number)))+
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=mean_gen_pop_1), col='black')+
  guides(col='none')+
  facet_wrap(~carr_cap, ncol=1)+
  NULL

ggplot()+
  geom_line(data=growth_prop[growth_prop$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=num_gen_2, col=as.factor(sim_number)))+
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=mean_gen_pop_2), col='black')+
  guides(col='none')+
  facet_wrap(~carr_cap, ncol=1)+
  NULL


mean_sel_r_test=data.frame(transfer=seq(1,max(mean_num_gen$transfer)),
                           gen_petite=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,]$mean_gen_pop_1,
                           gen_second=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,]$mean_gen_pop_2,
                           gen_population=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,]$mean_gen_all_pop)

mean_sel_r_test$gen_petite/mean_sel_r_test$gen_second

mean_sel_r_test$gen_second-mean_sel_r_test$gen_petite

ggplot()+
  geom_line(data=mean_sel_r_test, aes(x=transfer, y=gen_petite, col='Gens Petite'), col='red')+
  geom_line(data=mean_sel_r_test, aes(x=transfer, y=gen_second, col='Gens Petite'), col='blue')+
  geom_line(data=mean_sel_r_test, aes(x=transfer, y=gen_population, col='Difference'), col='black')+
  NULL


ggplot() +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,],
              aes(x=transfer, ymin=mean_gen_pop_1-sd_pop1, ymax=mean_gen_pop_1+sd_pop1),
              alpha=0.2, fill='red') +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,],
              aes(x=transfer, ymin=mean_gen_pop_2-sd_pop2, ymax=mean_gen_pop_2+sd_pop2),
              alpha=0.2, fill='blue') +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,],
              aes(x=transfer, ymin=mean_gen_all_pop-sd_pop_all, ymax=mean_gen_all_pop+sd_pop_all),
              alpha=0.2, fill='black') +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=mean_gen_pop_1, color="Petite")) +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=mean_gen_pop_2, color="Petite w/o Delay")) +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap==plot_carr_cap,], 
            aes(x=transfer, y=mean_gen_all_pop, color="Whole Population")) +
  scale_color_manual(values=c("Petite"="red", 
                              "Petite w/o Delay"="blue", 
                              "Whole Population"="black"),
                     name="Population Type") +
  labs(x="Transfer", y="Mean Number of Generations") +
  theme_classic()+
  NULL


#### Percentage of petite w/o Delay ####

growth_prop$percentage_petite_no_delay = growth_prop$cells_pop2_b/growth_prop$total_cells_b

ggplot(growth_prop[growth_prop$carr_cap==plot_carr_cap,], 
       aes(x=transfer, y=percentage_petite_no_delay, group=as.factor(sim_number)))+
  geom_line(color='blue', alpha=0.1)+
  stat_summary(aes(group=1), fun=mean, geom='line', color='black', linewidth=1)+
  labs(y="Proportion Petite w/o Delay", x='Transfer')+
  NULL

#### Filtered selection rates ####

temp_filt_growth_prop=growth_prop[growth_prop$transfer>1,]
temp_filt_settling_prop=settling_prop #no filtering for the settling phase


filt_combined_sel_r=data.frame(phase=c(rep('Growth', dim(temp_filt_growth_prop)[1]), rep('Settling', dim(temp_filt_settling_prop)[1])),
                               sel_r=c(temp_filt_growth_prop$cell_sel_r, temp_filt_settling_prop$cell_sel_r),
                               carr_cap=c(temp_filt_growth_prop$carr_cap, temp_filt_settling_prop$carr_cap),
                               sim_number=c(temp_filt_growth_prop$sim_number, temp_filt_settling_prop$sim_number))

#[filt_combined_sel_r$carr_cap!='10k',]
ggplot(filt_combined_sel_r, aes(x=phase, y=sel_r, fill=phase))+
  geom_violin()+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('')+
  guides(fill='none')+
  scale_y_continuous(breaks = seq(-0.5, 1.5, by = 0.5))+
  # coord_cartesian(ylim = c(-0.25, 1))+
  stat_summary(fun='mean', geom='crossbar')+
  facet_wrap(~carr_cap)+
  # scale_x_discrete(labels=c("Growth" = paste('Growth\nn=',sum(filt_combined_sel_r$phase=='Growth')), 
  #                           "Settling" = paste('Settling\nn=',sum(filt_combined_sel_r$phase=='Settling'))))+
  NULL

#[filt_combined_sel_r$carr_cap!='10k',]
ggplot(filt_combined_sel_r, aes(x=phase, y=sel_r, col=phase))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  facet_wrap(~carr_cap)+
  NULL

ggplot(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',], 
       aes(x=carr_cap, y=sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10k')),
                            "100k" = paste('100k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',] %>%
  group_by(carr_cap) %>%
  summarize(mean=mean(sel_r),
            var=var(sel_r),
            sd=sd(sel_r),
            CV=sd(sel_r)/mean(sel_r))
#   carr_cap   mean   var    sd    CV
# 1 10k      -0.631 2.87  1.69  -2.68
# 2 100k      0.307 1.61  1.27   4.14
# 3 1mill     0.356 0.656 0.810  2.27
# 4 10mill    
# The variation and standard deviation clearly get reduced, but also the mean
#continues increasing as the carrying capacity increases



# Plot of growth rates
ggplot(filt_combined_sel_r[filt_combined_sel_r$phase=='Growth',], 
       aes(x=carr_cap, y=sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10k')),
                            "100k" = paste('100k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL


filt_combined_sel_r[filt_combined_sel_r$phase=='Growth',] %>%
  group_by(carr_cap) %>%
  summarize(mean=mean(sel_r),
            var=var(sel_r),
            sd=sd(sel_r),
            CV=sd(sel_r)/mean(sel_r))
#   carr_cap  mean      var     sd     CV
# 1 10k      0.558 0.00588  0.0767 0.137 
# 2 100k     0.644 0.000991 0.0315 0.0489
# 3 1mill    0.641 0.000619 0.0249 0.038
# 4 10mill   


# Using mean of the simulations

mean_filt_sel_r=filt_combined_sel_r %>%
  group_by(carr_cap, sim_number, phase) %>%
  summarise(mean_sel_r=mean(sel_r))


ggplot(mean_filt_sel_r[mean_filt_sel_r$phase=='Growth',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='10k')),
                            "100k" = paste('100k\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

mean_filt_sel_r[mean_filt_sel_r$phase=='Growth',] %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r))
# carr_cap  mean
# 1 10k      0.546
# 2 100k     0.643
# 3 1mill    0.641
# 4 10mill   

ggplot(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='10k')),
                            "100k" = paste('100k\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',] %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r))
# carr_cap   mean
# 1 10k      -0.548
# 2 100k      0.683
# 3 1mill     0.431
# 4 10mill    

# Paper plot using mean of the simulations

mean_sim_sel_r=filt_combined_sel_r[filt_combined_sel_r$carr_cap==plot_carr_cap,] %>%
  group_by(phase, sim_number) %>%
  summarise(mean_sel_r=mean(sel_r))

ggplot(mean_sim_sel_r, aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL


ggplot(mean_sim_sel_r, aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  coord_cartesian(ylim = c(0, 1))+
  NULL


# Clusters per time ####

time_top=read.csv('petite_second_10mill_50sim/time_registry_all_sim.csv', header=TRUE)
time_top$total_time=time_top$time+((time_top$transfer-1)*24*60)
time_top$total_hours=time_top$total_time/60
summary(time_top)

table(time_top$sim_number)

ggplot(time_top[time_top$sim_number<=5,], 
       aes(x=total_hours, y=num_clusters, col=strain))+
  geom_line()+
  facet_wrap(~sim_number, ncol=1)+
  ylab('Number of Clusters')+
  xlab('Hours')+
  scale_color_discrete(name = "Strain", labels = c("Ancestor", "Ancestor w/o Delay"))+
  theme(legend.position = "bottom", legend.box = "vertical")+
  NULL


ggplot(time_top[time_top$sim_number<=1,], 
       aes(x=total_hours, y=num_clusters, col=strain))+
  geom_line()+
  # facet_wrap(~sim_number, ncol=1)+
  ylab('Number of Clusters')+
  xlab('Hours')+
  scale_color_discrete(name = "Strain", labels = c("Ancestor", "Ancestor w/o Delay"))+
  theme(legend.position = "bottom", legend.box = "vertical")+
  NULL

# Bootstrap settling phase ####


#### Settling fitness ####

settling_prop_boot=data.frame()

for(j in c("10mill-boot", "100mill-boot", "1bill-boot")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim/proportions_settling_selec_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  #filtering out the row if any population went extinct
  temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
  settling_prop_boot=rbind(settling_prop_boot, temp_df)
}
write.csv(settling_prop_boot, file='set_sim_bootstrap_settling_registry_30jan2026.csv', row.names = FALSE)

settling_prop_boot=read.csv('set_sim_bootstrap_settling_registry_30jan2026.csv', header=TRUE)
settling_prop_boot$carr_cap=factor(settling_prop_boot$carr_cap, levels=c("10mill-boot", "100mill-boot", "1bill-boot"))
summary(settling_prop_boot)

#Selection rate

settling_prop_boot$cell_sel_r=(log(settling_prop_boot$cells_pop2_a/settling_prop_boot$cells_pop2_b, base=2)-log(settling_prop_boot$cells_pop1_a/settling_prop_boot$cells_pop1_b, base=2))

ggplot(settling_prop_boot, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Settling Selection Rate\n(Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

settling_prop_boot %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(cell_sel_r))
# carr_cap      mean
# 1 10mill-boot  0.434
# 2 100mill-boot 0.440
# 3 1bill-boot   0.434

mean_settling_prop_boot=settling_prop_boot %>%
  group_by(sim_number, carr_cap) %>%
  summarise(mean_sel=mean(cell_sel_r))

ggplot(mean_settling_prop_boot, aes(x=carr_cap, y=mean_sel, fill=carr_cap))+
  geom_violin()+
  guides(fill='none')+
  ylab("Settling Selection Rate\n(Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

mean_settling_prop_boot %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel))
# carr_cap      mean
# 1 10mill-boot  0.434
# 2 100mill-boot 0.440
# 3 1bill-boot   0.434


#### Comparison to normal simulations ####

summary(settling_prop)

summary(settling_prop_boot)

temp_settling_prop=settling_prop
colnames(temp_settling_prop)=c("sim_number", "bootstrap_num", "strain1", "strain2", "clusters_pop1_b", "clusters_pop2_b", 
                               "cells_pop1_b", "cells_pop2_b", "total_clusters_b", "total_cells_b", "clusters_pop1_a", 
                               "clusters_pop2_a", "cells_pop1_a", "cells_pop2_a", "total_clusters_a", "total_cells_a", 
                               "carr_cap", "cell_sel_r")

comp_sett_carr_cap=rbind(temp_settling_prop, settling_prop_boot)
summary(comp_sett_carr_cap)

ggplot(comp_sett_carr_cap[comp_sett_carr_cap$carr_cap!='10k',], 
       aes(x=carr_cap, y=cell_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  # scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10k')),
  #                           "100k" = paste('100k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='100k')),
  #                           "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='1mill')),
  #                           "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

comp_sett_carr_cap %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(cell_sel_r))
# carr_cap        mean
# 1 10k          -0.631
# 2 100k          0.307
# 3 1mill         0.356
# 4 10mill        
# 5 10mill-boot   0.434
# 6 100mill-boot  0.440
# 7 1bill-boot    0.434 


mean_comp_sett_carr_cap=comp_sett_carr_cap %>%
  group_by(carr_cap, sim_number) %>%
  summarise(mean_sel_r=mean(cell_sel_r))

ggplot(mean_comp_sett_carr_cap[mean_comp_sett_carr_cap$carr_cap!='10k',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  geom_jitter(alpha=0.5)+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='10mill')),
                            "10mill-boot" = paste('10mill-boot\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='10mill-boot')),
                            "100mill-boot" = paste('100mill-boot\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='100mill-boot')),
                            "1bill-boot" = paste('1bill-boot\nn=',sum(mean_comp_sett_carr_cap$carr_cap=='1bill-boot'))))+
  NULL
#note: there are no 50 data points in 100k because in some simulations the petite strain went extinct in the
#first generation

mean_comp_sett_carr_cap %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r),
            var=var(mean_sel_r),
            sd=sd(mean_sel_r),
            CV=sd(mean_sel_r)/mean(mean_sel_r),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)
# carr_cap       mean
# 1 10k          -0.548
# 2 100k          0.683
# 3 1mill         0.431
# 4 10mill        
# 5 10mill-boot   0.434
# 6 100mill-boot  0.440
# 7 1bill-boot    0.434 


anova_random=aov(mean_sel_r~carr_cap, data=mean_comp_sett_carr_cap)
summary(anova_random)
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# carr_cap      6   6.27  1.0445   4.197 0.000456 ***
# Residuals   294  73.17  0.2489
# the differences are significant
# F(6, 294)=4.197, p=0.000456

pairwise.t.test(mean_comp_sett_carr_cap$mean_sel_r, mean_comp_sett_carr_cap$carr_cap, p.adjust.method='bonferroni')
#              10k     100k    1mill   10mill  10mill-boot 100mill-boot
# 100k         2.3e-05 -       -       -       -           -           
# 1mill        6.5e-05 1.00000 -       -       -           -           
# 10mill       8.6e-05 1.00000 1.00000 -       -           -           
# 10mill-boot  0.00012 1.00000 1.00000 1.00000 -           -           
# 100mill-boot 0.00012 1.00000 1.00000 1.00000 1.00000     -           
# 1bill-boot   0.00013 1.00000 1.00000 1.00000 1.00000     1.00000
# 10k is different from all of the others, and there is no difference between all
# the other concentrations


anova_growth_random=aov(mean_sel_r~carr_cap, data=mean_filt_sel_r[mean_filt_sel_r$phase=='Growth' & mean_filt_sel_r$carr_cap!='10k',])
summary(anova_growth_random)
#              Df   Sum Sq   Mean Sq F value Pr(>F)
# carr_cap      2 0.000493 0.0002464   1.496  0.227
# Residuals   143 0.023548 0.0001647 
# non of the means are different from each other
# F(2, 143)=1.496, p=0.227




# Top Settling selection ####
#This results are doing settling selection starting all the clusters at the top instead of
#them starting at a random position. 


#### Growth phase ####

growth_prop_top=data.frame()

# , "10mill"
for(j in c("100k", "1mill")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim_top/proportions_growth_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  growth_prop_top=rbind(growth_prop_top, temp_df)
}
write.csv(growth_prop_top, file='set_sim_top_growth_registry_30jan2026.csv', row.names = FALSE)

growth_prop_top=read.csv('set_sim_top_growth_registry_30jan2026.csv', header=TRUE)
growth_prop_top$carr_cap=factor(growth_prop_top$carr_cap, levels=c("100k", "1mill", "10mill"))
summary(growth_prop_top)

#Selection rate

growth_prop_top$cell_sel_r=(log(growth_prop_top$cells_pop2_a/growth_prop_top$cells_pop2_b, base=2)-log(growth_prop_top$cells_pop1_a/growth_prop_top$cells_pop1_b, base=2))

ggplot(growth_prop_top, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Growth Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL


#### Settling fitness ####

settling_prop_top=data.frame()

# , "10mill"
for(j in c("100k", "1mill")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim_top/proportions_settling_selec_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  #filtering out the row if any population went extinct
  temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
  settling_prop_top=rbind(settling_prop_top, temp_df)
}
write.csv(settling_prop_top, file='set_sim_top_settling_registry_30jan2026.csv', row.names = FALSE)

settling_prop_top=read.csv("set_sim_top_settling_registry_30jan2026.csv", header=TRUE)
settling_prop_top$carr_cap=factor(settling_prop_top$carr_cap, levels=c("100k", "1mill", "10mill"))
summary(settling_prop_top)

#Selection rate

settling_prop_top$cell_sel_r=(log(settling_prop_top$cells_pop2_a/settling_prop_top$cells_pop2_b, base=2)-log(settling_prop_top$cells_pop1_a/settling_prop_top$cells_pop1_b, base=2))

ggplot(settling_prop_top, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL


#### Do selection rates change over time? ####

ggplot(growth_prop_top,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL

ggplot(growth_prop_top[growth_prop_top$transfer>1,],
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL
#Getting rid of the first selection rate

ggplot(settling_prop_top,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Settling Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~carr_cap, ncol=1)+
  NULL


#### Filtered selection rates ####

temp_filt_growth_prop_top=growth_prop_top[growth_prop_top$transfer>1,]
temp_filt_settling_prop_top=settling_prop_top #no filtering for the settling phase


filt_combined_sel_r_top=data.frame(phase=c(rep('Growth', dim(temp_filt_growth_prop_top)[1]), rep('Settling', dim(temp_filt_settling_prop_top)[1])),
                               sel_r=c(temp_filt_growth_prop_top$cell_sel_r, temp_filt_settling_prop_top$cell_sel_r),
                               carr_cap=c(temp_filt_growth_prop_top$carr_cap, temp_filt_settling_prop_top$carr_cap),
                               sim_number=c(temp_filt_growth_prop_top$sim_number, temp_filt_settling_prop_top$sim_number))


ggplot(filt_combined_sel_r_top, aes(x=phase, y=sel_r, fill=phase))+
  geom_violin()+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('')+
  guides(fill='none')+
  scale_y_continuous(breaks = seq(-0.5, 6, by = 0.5))+
  # coord_cartesian(ylim = c(-0.25, 1))+
  stat_summary(fun='mean', geom='crossbar')+
  facet_wrap(~carr_cap)+
  # scale_x_discrete(labels=c("Growth" = paste('Growth\nn=',sum(filt_combined_sel_r$phase=='Growth')), 
  #                           "Settling" = paste('Settling\nn=',sum(filt_combined_sel_r$phase=='Settling'))))+
  NULL


ggplot(filt_combined_sel_r_top, aes(x=phase, y=sel_r, col=phase))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-0.5, 6, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  facet_wrap(~carr_cap)+
  NULL

ggplot(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',], 
       aes(x=carr_cap, y=sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 6, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',] %>%
  group_by(carr_cap) %>%
  summarize(mean=mean(sel_r),
            var=var(sel_r),
            sd=sd(sel_r),
            CV=sd(sel_r)/mean(sel_r))
#   carr_cap  mean   var    sd    CV
# 1 100k      1.81 0.148 0.384 0.213
# 2 1mill     3.82 1.76  1.33  0.347
# 3 10mill    
# The mean keeps going up so 1 million and 10 million are different values


# Plot of growth rates
ggplot(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Growth',], 
       aes(x=carr_cap, y=sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL


filt_combined_sel_r_top[filt_combined_sel_r_top$phase=='Growth',] %>%
  group_by(carr_cap) %>%
  summarize(mean=mean(sel_r),
            var=var(sel_r),
            sd=sd(sel_r),
            CV=sd(sel_r)/mean(sel_r))
#   carr_cap  mean      var     sd     CV
# 1 100k     0.630 0.00122  0.0349 0.0554
# 2 1mill    0.636 0.000296 0.0172 0.0271
# 3 10mill   


# Using mean of the simulations

mean_filt_sel_r_top=filt_combined_sel_r_top %>%
  group_by(carr_cap, sim_number, phase) %>%
  summarise(mean_sel_r=mean(sel_r))


ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Growth',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Growth',] %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r))
# carr_cap  mean
# 1 100k     0.630
# 2 1mill    0.636
# 3 10mill   

ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 6, by = 0.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',] %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r))
# carr_cap  mean
# 1 100k      1.81
# 2 1mill     3.96
# 3 10mill    



# Bootstrap top settling ####


#### Settling fitness

settling_prop_boot_top=data.frame()

for(j in c("10mill-boot", "100mill-boot", "1bill-boot")){
  temp_df=read.csv(paste('petite_second_', j, "_50sim_top/proportions_settling_selec_registry_all_sim.csv", sep=""), header=TRUE)
  temp_df$carr_cap=j
  #filtering out the row if any population went extinct
  temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
  settling_prop_boot_top=rbind(settling_prop_boot_top, temp_df)
}
write.csv(settling_prop_boot_top,file="set_sim_bootstrap_top_settling_registry_30jan2026.csv", row.names = FALSE)

settling_prop_boot_top=read.csv("set_sim_bootstrap_top_settling_registry_30jan2026.csv", header=TRUE)
settling_prop_boot_top$carr_cap=factor(settling_prop_boot_top$carr_cap, levels=c("10mill-boot", "100mill-boot", "1bill-boot"))
summary(settling_prop_boot_top)

#Selection rate

settling_prop_boot_top$cell_sel_r=(log(settling_prop_boot_top$cells_pop2_a/settling_prop_boot_top$cells_pop2_b, base=2)-log(settling_prop_boot_top$cells_pop1_a/settling_prop_boot_top$cells_pop1_b, base=2))

ggplot(settling_prop_boot_top, 
       aes(x=carr_cap, y=cell_sel_r, fill=carr_cap)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Settling Selection Rate\n(Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

settling_prop_boot_top %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(cell_sel_r))
# carr_cap      mean
# 1 10mill-boot   4.50
# 2 100mill-boot  4.62
# 3 1bill-boot    4.53

mean_settling_prop_boot_top=settling_prop_boot_top %>%
  group_by(sim_number, carr_cap) %>%
  summarise(mean_sel=mean(cell_sel_r))

ggplot(mean_settling_prop_boot_top, aes(x=carr_cap, y=mean_sel, fill=carr_cap))+
  geom_violin()+
  guides(fill='none')+
  ylab("Settling Selection Rate\n(Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

mean_settling_prop_boot_top %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel))
# carr_cap      mean
# 1 10mill-boot   4.50
# 2 100mill-boot  4.62
# 3 1bill-boot    4.53


#### Comparison to normal simulations ####

summary(settling_prop_top)

summary(settling_prop_boot_top)

temp_settling_prop_top=settling_prop_top
colnames(temp_settling_prop_top)=c("sim_number", "bootstrap_num", "strain1", "strain2", "clusters_pop1_b", "clusters_pop2_b", 
                               "cells_pop1_b", "cells_pop2_b", "total_clusters_b", "total_cells_b", "clusters_pop1_a", 
                               "clusters_pop2_a", "cells_pop1_a", "cells_pop2_a", "total_clusters_a", "total_cells_a", 
                               "carr_cap", "cell_sel_r")

comp_sett_carr_cap_top=rbind(temp_settling_prop_top, settling_prop_boot_top)
summary(comp_sett_carr_cap_top)

ggplot(comp_sett_carr_cap_top[comp_sett_carr_cap_top$carr_cap!='10k',], 
       aes(x=carr_cap, y=cell_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  # scale_x_discrete(labels=c("10k" = paste('10k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10k')),
  #                           "100k" = paste('100k\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='100k')),
  #                           "1mill" = paste('1mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='1mill')),
  #                           "10mill" = paste('10mill\nn=',sum(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',]$carr_cap=='10mill'))))+
  NULL

comp_sett_carr_cap_top %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(cell_sel_r))
# 1 100k          1.81
# 2 1mill         3.82
# 3 10mill
# 4 10mill-boot   4.50
# 5 100mill-boot  4.62
# 6 1bill-boot    4.53

mean_comp_sett_carr_cap_top=comp_sett_carr_cap_top %>%
  group_by(carr_cap, sim_number) %>%
  summarise(mean_sel_r=mean(cell_sel_r))

ggplot(mean_comp_sett_carr_cap_top, 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  # geom_beeswarm()+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  geom_jitter(alpha=0.5)+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  scale_x_discrete(labels=c("100k" = paste('100k\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='100k')),
                            "1mill" = paste('1mill\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='1mill')),
                            "10mill" = paste('10mill\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='10mill')),
                            "10mill-boot" = paste('10mill-boot\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='10mill-boot')),
                            "100mill-boot" = paste('100mill-boot\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='100mill-boot')),
                            "1bill-boot" = paste('1bill-boot\nn=',sum(mean_comp_sett_carr_cap_top$carr_cap=='1bill-boot'))))+
  NULL


mean_comp_sett_carr_cap_top %>%
  group_by(carr_cap) %>%
  summarise(mean=mean(mean_sel_r),
            var=var(mean_sel_r),
            sd=sd(mean_sel_r),
            CV=sd(mean_sel_r)/mean(mean_sel_r),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)
# carr_cap      mean
# 1 100k          1.81
# 2 1mill         3.96
# 3 10mill        
# 4 10mill-boot   4.50
# 5 100mill-boot  4.62
# 6 1bill-boot    4.53 

anova_growth_top=aov(mean_sel_r~carr_cap, data=mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Growth',])
summary(anova_growth_top)
#             Df   Sum Sq   Mean Sq F value Pr(>F)
# carr_cap     2 0.000185 9.265e-05   0.392  0.677
# Residuals   77 0.018182 2.361e-04 
# no difference is significant
# F(2, 77)=0.392, p=0.677

anova_top=aov(mean_sel_r~carr_cap, data=mean_comp_sett_carr_cap_top)
summary(anova_top)
#              Df Sum Sq Mean Sq F value Pr(>F)    
# carr_cap      5   99.1  19.821   37.34 <2e-16 ***
# Residuals   224  118.9   0.531   
# There is a group that is different from the others
# F(5, 224)=37.34, p<2e-16

pairwise.t.test(mean_comp_sett_carr_cap_top$mean_sel_r, mean_comp_sett_carr_cap_top$carr_cap, p.adjust.method='bonferroni')
#              100k    1mill   10mill 10mill-boot 100mill-boot
# 1mill        1.5e-11 -       -      -           -           
# 10mill       < 2e-16 1.4e-08 -      -           -           
# 10mill-boot  < 2e-16 8.7e-09 1      -           -           
# 100mill-boot < 2e-16 2.4e-09 1      1           -           
# 1bill-boot   < 2e-16 2.6e-09 1      1           1
# 100k and 1mill has a different mean that the others, so at least 10 millions is needed to estimate the mean








# Supplementary Figure 8 ####


supp_mean_growth=ggplot(mean_filt_sel_r[mean_filt_sel_r$phase=='Growth' & mean_filt_sel_r$carr_cap!='10k',], 
                        aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed', linewidth=0.3)+
  geom_jitter(alpha=0.5)+
  ylab('Growth Selection Rate \n(Petite w/o delay - Petite)\nper selection phase')+
  # xlab('Carrying Capacity')+
  xlab('')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  coord_cartesian(ylim = c(0, 1))+
  theme_classic(base_size = 11)+
  scale_x_discrete(labels=c("100k" = "1e5",
                            "1mill" = "1e6",
                            "10mill" = "1e7"))+
  NULL
supp_mean_growth


supp_mean_settling=ggplot(mean_comp_sett_carr_cap[mean_comp_sett_carr_cap$carr_cap!='10k',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed', linewidth=0.3)+
  geom_jitter(alpha=0.5)+
  ylab('Settling Selection Rate \n(Petite w/o delay - Petite)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-0, 6, by = 1))+
  coord_cartesian(ylim = c(-0.5, 6))+
  theme_classic(base_size = 11)+
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(labels=c("100k" = expression(10^5),
                            "1mill" = expression(10^6),
                            "10mill" = expression(10^7),
                            "10mill-boot" = expression(paste(10^7, " boot")),
                            "100mill-boot" = expression(paste(10^8, " boot")),
                            "1bill-boot" = expression(paste(10^9, " boot"))))+
  NULL
supp_mean_settling


supp_mean_growth_top=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Growth',], 
       aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed', linewidth=0.3)+
  geom_jitter(alpha=0.5)+
  # ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  ylab('')+
  # xlab('Carrying Capacity')+
  xlab('')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  coord_cartesian(ylim = c(0, 1))+
  theme_classic(base_size = 11)+
  scale_x_discrete(labels=c("100k" = "1e5",
                            "1mill" = "1e6",
                            "10mill" = "1e7"))+
  NULL
supp_mean_growth_top


supp_mean_settling_top=ggplot(mean_comp_sett_carr_cap_top, 
                              aes(x=carr_cap, y=mean_sel_r, col=carr_cap))+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed', linewidth=0.3)+
  geom_jitter(alpha=0.5)+
  ylab('')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-0, 6, by = 1))+
  coord_cartesian(ylim = c(-0.5, 6))+
  theme_classic(base_size = 11)+
  scale_x_discrete(labels=c("100k" = expression(10^5),
                            "1mill" = expression(10^6),
                            "10mill" = expression(10^7),
                            "10mill-boot" = expression(paste(10^7, " boot")),
                            "100mill-boot" = expression(paste(10^8, " boot")),
                            "1bill-boot" = expression(paste(10^9, " boot"))))+
  NULL
supp_mean_settling_top


# supp_carr_cap=plot_grid(supp_mean_growth, supp_mean_growth_top, supp_mean_settling, supp_mean_settling_top, 
#                         labels=c('A', 'B', 'C', 'D'), ncol=2, align='hv')
# supp_carr_cap
# 
# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig8_carrying_capacity_test_23apr2025.png',
#        plot=supp_carr_cap, dpi='retina', width=9, height=7)

supp_carr_cap=plot_grid(supp_mean_settling, supp_mean_settling_top, 
                        labels=c('A', 'B'), ncol=2, align='hv')
supp_carr_cap

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig9_carrying_capacity_test_6june2025.png',
#        plot=supp_carr_cap, dpi='retina', width=9, height=3.5)


# Supplementary Figure 6 ####

growth_s_rate=ggplot(growth_prop[growth_prop$carr_cap=='10mill',],
                     aes(x=transfer, y=cell_sel_r, group=sim_number,))+
  geom_line(color='#657060FF')+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Petite w/o delay / Petite)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  theme_classic(base_size = 12)+
  NULL
growth_s_rate

settling_s_rate=ggplot(settling_prop[growth_prop$carr_cap=='10mill',],
                       aes(x=transfer, y=cell_sel_r, group=sim_number))+
  geom_line(color='#657060FF')+
  guides(col='none')+
  ylab("Settling Selection Rate \n(Petite w/o delay / Petite)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  theme_classic(base_size = 12)+
  NULL
settling_s_rate

p_num_gen=ggplot() +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',],
              aes(x=transfer, ymin=mean_gen_pop_1-sd_pop1, ymax=mean_gen_pop_1+sd_pop1),
              alpha=0.5, fill='#F0D77BFF') +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',],
              aes(x=transfer, ymin=mean_gen_pop_2-sd_pop2, ymax=mean_gen_pop_2+sd_pop2),
              alpha=0.5, fill='#B4DAE5FF') +
  geom_ribbon(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',],
              aes(x=transfer, ymin=mean_gen_all_pop-sd_pop_all, ymax=mean_gen_all_pop+sd_pop_all),
              alpha=0.5, fill='black') +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',], 
            aes(x=transfer, y=mean_gen_pop_1, color="Petite")) +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',], 
            aes(x=transfer, y=mean_gen_pop_2, color="Petite w/o Delay")) +
  geom_line(data=mean_num_gen[mean_num_gen$carr_cap=='10mill',], 
            aes(x=transfer, y=mean_gen_all_pop, color="Whole Population")) +
  scale_color_manual(values=c("Petite"="#F0D77BFF", 
                              "Petite w/o Delay"="#B4DAE5FF", 
                              "Whole Population"="black"),
                     name="Population Type") +
  labs(x="Transfer", y="Mean Number of\nGenerations") +
  theme_classic(base_size = 12)+
  NULL
p_num_gen


supp_s_rate=plot_grid(growth_s_rate,settling_s_rate, p_num_gen,
                      labels=c('A', 'B', 'C'), ncol=1)
supp_s_rate

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig7_selection_rate_per_transfer_3oct2025.png',
#        plot=supp_s_rate, dpi='retina', width=6, height=8)



# Figure 6 ####

fig6a_v2=ggplot(mean_sim_sel_r, 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay vs Petite)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 6, by = 1))+
  coord_cartesian(ylim = c(0, 6))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle('Random Starting Position')+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig6a_v2


fig6b_v2=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$carr_cap=='10mill',], 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay - Petite)\nper selection phase')+
  # ylab('')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 6, by = 1))+
  coord_cartesian(ylim = c(0, 6))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle("Top Starting Position")+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig6b_v2


figure_6_v3=plot_grid(fig6a_v2, fig6b_v2, labels=c('A', 'B'), ncol=1, align='hv')
figure_6_v3

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_selection_rates_23apr2025.png',
#        plot=figure_6_v3, dpi='retina', width=4, height=6)

# Merged figure 7 and supp Fig 9 ####

fig_a=ggplot(mean_sim_sel_r, 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay vs Petite)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 6, by = 1))+
  coord_cartesian(ylim = c(0, 6))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle('Random Starting Position')+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig_a


fig_c=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$carr_cap=='10mill',], 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay - Petite)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 6, by = 1))+
  coord_cartesian(ylim = c(0, 6))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle("Top Starting Position")+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig_c
