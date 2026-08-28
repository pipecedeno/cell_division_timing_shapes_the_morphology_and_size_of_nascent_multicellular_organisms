
# Date: 2024Feb2

library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(ghibli)
library(cowplot)
library(effsize)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/settling_selection_simulations_30jan2026/")

# Random Position ####

#### Growth phase ####

# growth_prop=data.frame()
# temp_df=read.csv("petite_second-15min_10mill_50sim/proportions_growth_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Equal Growth Rates'
# growth_prop=rbind(growth_prop, temp_df)
# temp_df=read.csv("petite_second_10mill_50sim/proportions_growth_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Normal'
# growth_prop=rbind(growth_prop, temp_df)
# growth_prop$experiment=factor(growth_prop$experiment, levels=c('Normal', 'Equal Growth Rates'))
# summary(growth_prop)
# write.csv(growth_prop, file="set_sim_growth_registry_equal_growth_rate_12mar2026.csv", row.names = FALSE)

growth_prop=read.csv("set_sim_growth_registry_equal_growth_rate_12mar2026.csv", header=TRUE)

growth_prop$experiment=factor(growth_prop$experiment, levels=c('Normal', 'Equal Growth Rates'))
summary(growth_prop)

#Selection rate

growth_prop$cell_sel_r=(log(growth_prop$cells_pop2_a/growth_prop$cells_pop2_b, base=2)-log(growth_prop$cells_pop1_a/growth_prop$cells_pop1_b, base=2))

ggplot(growth_prop, 
       aes(x=experiment, y=cell_sel_r, fill=experiment)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

#### Transfers to extinction ####

extinction_inf=growth_prop %>%
  group_by(sim_number, experiment) %>%
  summarize(total_transfers=max(transfer))

extinction_inf %>%
  group_by(experiment) %>%
  summarize(mean_transfers=mean(total_transfers))
#   experiment         mean_transfers
# 1 Normal                       9.82
# 2 Equal Growth Rates          22.0

ggplot(extinction_inf[extinction_inf$experiment=='Equal Growth Rates',], 
       aes(x=total_transfers, col=experiment, fill=experiment))+
  geom_histogram(position="identity", alpha=0.5, bins=10)+
  theme_classic(base_size=16)+
  guides(fill='none', col='none')+
  labs(x="Transfers until Fixation", y="Count")+
  # facet_wrap(~experiment)+
  geom_vline(xintercept=mean(extinction_inf[extinction_inf$experiment=='Equal Growth Rates',]$total_transfers), linetype='dashed', size=1.2)+
  NULL


ggplot(extinction_inf, 
       aes(x=total_transfers, col=experiment, fill=experiment))+
  geom_histogram(position="identity", alpha=0.5, bins=10)+
  theme_classic(base_size=16)+
  # guides(fill='none', col='none')+
  labs(x="Transfers until Fixation", y="Count")+
  # facet_wrap(~experiment)+
  # geom_vline(xintercept=mean(extinction_inf[extinction_inf$experiment=='Equal Growth Rates',]$total_transfers), linetype='dashed', size=1.2)+
  NULL


#### Settling fitness ####

# settling_prop=data.frame()
# 
# temp_df=read.csv("petite_second-15min_10mill_50sim/proportions_settling_selec_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Equal Growth Rates'
# temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
# settling_prop=rbind(settling_prop, temp_df)
# temp_df=read.csv("petite_second_10mill_50sim/proportions_settling_selec_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Normal'
# temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
# settling_prop=rbind(settling_prop, temp_df)
# settling_prop$experiment=factor(settling_prop$experiment, levels=c('Normal', 'Equal Growth Rates'))
# summary(settling_prop)
# 
# write.csv(settling_prop, file="set_sim_settling_registry_equal_growth_rate_12mar2026.csv", row.names=FALSE)


settling_prop=read.csv("set_sim_settling_registry_equal_growth_rate_12mar2026.csv", header=TRUE)

settling_prop$experiment=factor(settling_prop$experiment, levels=c('Normal', 'Equal Growth Rates'))
summary(settling_prop)

#Selection rate

settling_prop$cell_sel_r=(log(settling_prop$cells_pop2_a/settling_prop$cells_pop2_b, base=2)-log(settling_prop$cells_pop1_a/settling_prop$cells_pop1_b, base=2))

ggplot(settling_prop, 
       aes(x=experiment, y=cell_sel_r, fill=experiment)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL


#### Do selection rates change over time? ####

ggplot(growth_prop,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  # scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~experiment, ncol=1)+
  NULL

ggplot(growth_prop[growth_prop$transfer>1,],
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~experiment, ncol=1)+
  NULL
#Getting rid of the first selection rate

ggplot(settling_prop,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  guides(col='none')+
  ylab("Settling Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  # scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~experiment, ncol=1)+
  NULL


#### Filtered selection rates ####

temp_filt_growth_prop=rbind(growth_prop[growth_prop$experiment=='Normal' & growth_prop$transfer>1,],
                            growth_prop[growth_prop$experiment=='Equal Growth Rates',])
temp_filt_settling_prop=settling_prop #no filtering for the settling phase


filt_combined_sel_r=data.frame(phase=c(rep('Growth', dim(temp_filt_growth_prop)[1]), rep('Settling', dim(temp_filt_settling_prop)[1])),
                               sel_r=c(temp_filt_growth_prop$cell_sel_r, temp_filt_settling_prop$cell_sel_r),
                               experiment=c(temp_filt_growth_prop$experiment, temp_filt_settling_prop$experiment),
                               sim_number=c(temp_filt_growth_prop$sim_number, temp_filt_settling_prop$sim_number))


ggplot(filt_combined_sel_r, aes(x=phase, y=sel_r, fill=phase))+
  geom_violin()+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('')+
  guides(fill='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  stat_summary(fun='mean', geom='crossbar')+
  facet_wrap(~experiment)+
  # scale_x_discrete(labels=c("Growth" = paste('Growth\nn=',sum(filt_combined_sel_r$phase=='Growth')), 
  #                           "Settling" = paste('Settling\nn=',sum(filt_combined_sel_r$phase=='Settling'))))+
  NULL


ggplot(filt_combined_sel_r, aes(x=phase, y=sel_r, col=phase))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  facet_wrap(~experiment)+
  NULL

ggplot(filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',], 
       aes(x=experiment, y=sel_r, col=experiment))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Carrying Capacity')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL

filt_combined_sel_r[filt_combined_sel_r$phase=='Settling',] %>%
  group_by(experiment) %>%
  summarize(mean=mean(sel_r),
            var=var(sel_r),
            sd=sd(sel_r),
            CV=sd(sel_r)/mean(sel_r))
#   experiment          mean   var    sd    CV
# 1 Normal             0.401 0.335 0.579  1.44
# 2 Equal Growth Rates 0.468 0.266 0.516  1.10

# Using mean of the simulations

mean_filt_sel_r=filt_combined_sel_r %>%
  group_by(experiment, sim_number, phase) %>%
  summarise(mean_sel_r=mean(sel_r))


ggplot(mean_filt_sel_r[mean_filt_sel_r$phase=='Growth',], 
       aes(x=experiment, y=mean_sel_r, col=experiment))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Test Condition')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL

mean_filt_sel_r[mean_filt_sel_r$phase=='Growth',] %>%
  group_by(experiment) %>%
  summarise(mean=mean(mean_sel_r))
#   experiment            mean
# 1 Normal              0.642 
# 2 Equal Growth Rates -0.0669

ggplot(mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',], 
       aes(x=experiment, y=mean_sel_r, col=experiment))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Test Condition')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL

mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',] %>%
  group_by(experiment) %>%
  summarise(mean=mean(mean_sel_r))
#   experiment          mean
# 1 Normal             0.435
# 2 Equal Growth Rates 0.484


ggplot(mean_filt_sel_r, 
       aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection Phase')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  facet_wrap(~experiment)+
  NULL

mean_filt_sel_r %>%
  group_by(experiment, phase) %>%
  summarise(mean=mean(mean_sel_r))
#   experiment         phase       mean
# 1 Normal             Growth    0.642 
# 2 Normal             Settling  0.435 
# 3 Equal Growth Rates Growth   -0.0669
# 4 Equal Growth Rates Settling  0.484

mean_filt_sel_r %>%
  group_by(experiment, phase) %>%
  summarise(mean=mean(mean_sel_r),
            var=var(mean_sel_r),
            sd=sd(mean_sel_r),
            CV=sd(mean_sel_r)/mean(mean_sel_r),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)





# Top Settling selection ####
#This results are doing settling selection starting all the clusters at the top instead of
#them starting at a random position. 


#### Growth phase ####

# growth_prop_top=data.frame()
# 
# temp_df=read.csv("petite_second-15min_10mill_50sim_top/proportions_growth_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Equal Growth Rates'
# growth_prop_top=rbind(growth_prop_top, temp_df)
# temp_df=read.csv("petite_second_10mill_50sim_top/proportions_growth_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Normal'
# growth_prop_top=rbind(growth_prop_top, temp_df)
# growth_prop_top$experiment=factor(growth_prop_top$experiment, levels=c('Normal', 'Equal Growth Rates'))
# summary(growth_prop_top)
# 
# write.csv(growth_prop_top, file="set_sim_top_growth_registry_equal_growth_rate_12mar2026.csv", row.names = FALSE)


growth_prop_top=read.csv("set_sim_top_growth_registry_equal_growth_rate_12mar2026.csv", header=TRUE)

growth_prop_top$experiment=factor(growth_prop_top$experiment, levels=c('Normal', 'Equal Growth Rates'))
summary(growth_prop_top)

growth_prop_top$cell_sel_r=(log(growth_prop_top$cells_pop2_a/growth_prop_top$cells_pop2_b, base=2)-log(growth_prop_top$cells_pop1_a/growth_prop_top$cells_pop1_b, base=2))

ggplot(growth_prop_top, 
       aes(x=experiment, y=cell_sel_r, fill=experiment)) +
  geom_violin()+
  guides(fill='none')+
  ylab("Selection Rate (Ancestor w/o delay / Ancestor)")+
  xlab('')+
  stat_summary(fun='mean', geom='crossbar')+
  NULL


#### Transfers to extinction ####

extinction_inf_top=growth_prop_top %>%
  group_by(sim_number, experiment) %>%
  summarize(total_transfers=max(transfer))

extinction_inf_top %>%
  group_by(experiment) %>%
  summarize(mean_transfers=mean(total_transfers))
#   experiment         mean_transfers
# 1 Normal                       2.4 
# 2 Equal Growth Rates           2.62

ggplot(extinction_inf_top[extinction_inf_top$experiment=='Equal Growth Rates',], 
       aes(x=total_transfers, col=experiment, fill=experiment))+
  geom_histogram(position="identity", alpha=0.5, bins=10)+
  theme_classic(base_size=16)+
  guides(fill='none', col='none')+
  labs(x="Transfers until Fixation", y="Count")+
  # facet_wrap(~experiment)+
  geom_vline(xintercept=mean(extinction_inf_top[extinction_inf_top$experiment=='Equal Growth Rates',]$total_transfers), linetype='dashed', size=1.2)+
  NULL


ggplot(extinction_inf_top, 
       aes(x=total_transfers, col=experiment, fill=experiment))+
  geom_histogram(position="identity", alpha=0.5, bins=10)+
  theme_classic(base_size=16)+
  # guides(fill='none', col='none')+
  labs(x="Transfers until Fixation", y="Count")+
  # facet_wrap(~experiment)+
  # geom_vline(xintercept=mean(extinction_inf_top[extinction_inf_top$experiment=='Equal Growth Rates',]$total_transfers), linetype='dashed', size=1.2)+
  NULL

#### Settling fitness ####

# settling_prop_top=data.frame()
# 
# temp_df=read.csv("petite_second-15min_10mill_50sim_top/proportions_settling_selec_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Equal Growth Rates'
# temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
# settling_prop_top=rbind(settling_prop_top, temp_df)
# temp_df=read.csv("petite_second_10mill_50sim_top/proportions_settling_selec_registry_all_sim.csv", header=TRUE)
# temp_df$experiment='Normal'
# temp_df=temp_df[temp_df$cells_pop1_a!=0 & temp_df$cells_pop2_a!=0,]
# settling_prop_top=rbind(settling_prop_top, temp_df)
# settling_prop_top$experiment=factor(settling_prop_top$experiment, levels=c('Normal', 'Equal Growth Rates'))
# summary(settling_prop_top)
# 
# write.csv(settling_prop_top, file="set_sim_top_settling_registry_equal_growth_rate_12mar2026.csv", row.names = FALSE)


settling_prop_top=read.csv("set_sim_top_settling_registry_equal_growth_rate_12mar2026.csv", header=TRUE)

settling_prop_top$experiment=factor(settling_prop_top$experiment, levels=c('Normal', 'Equal Growth Rates'))
summary(settling_prop_top)

#Selection rate

settling_prop_top$cell_sel_r=(log(settling_prop_top$cells_pop2_a/settling_prop_top$cells_pop2_b, base=2)-log(settling_prop_top$cells_pop1_a/settling_prop_top$cells_pop1_b, base=2))

ggplot(settling_prop_top, 
       aes(x=experiment, y=cell_sel_r, fill=experiment)) +
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
  guides(col='none')+
  ylab("Growth Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~experiment)+
  NULL


ggplot(settling_prop_top,
       aes(x=transfer, y=cell_sel_r, col=as.factor(sim_number)))+
  geom_line()+
  geom_point()+
  guides(col='none')+
  ylab("Settling Selection Rate \n(Ancestor w/o delay / Ancestor)")+
  xlab('Transfer')+
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
  facet_wrap(~experiment, ncol=1)+
  NULL

#### Combined selection rates ####

growth_prop_top$phase='Growth'
settling_prop_top$phase='Settling'

temp_growth_prop_top=rbind(growth_prop_top[growth_prop_top$experiment=='Equal Growth Rates', -3],
                           growth_prop_top[growth_prop_top$experiment=='Normal' & growth_prop_top$transfer>1, -3])

combined_sel_r_top=rbind(temp_growth_prop_top, settling_prop_top) 
summary(combined_sel_r_top)

# Using mean of the simulations

mean_filt_sel_r_top=combined_sel_r_top %>%
  group_by(experiment, sim_number, phase) %>%
  summarise(mean_sel_r=mean(cell_sel_r))


ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Growth',], 
       aes(x=experiment, y=mean_sel_r, col=experiment))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black', linetype='dashed')+
  ylab('Growth Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Test Condition')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL


ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',], 
       aes(x=experiment, y=mean_sel_r, col=experiment))+
  # geom_beeswarm()+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Settling Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Test Condition')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  NULL


ggplot(mean_filt_sel_r_top, 
       aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  ylab('Selection Rate \n(Ancestor w/o delay - Ancestor)\nper selection phase')+
  xlab('Selection Phase')+
  guides(col='none')+
  # scale_y_continuous(breaks = seq(-1.5, 2.5, by = 0.5))+
  # coord_cartesian(ylim = c(-1.5, 2.5))+
  facet_wrap(~experiment)+
  NULL

mean_filt_sel_r_top %>%
  group_by(experiment, phase) %>%
  summarise(mean=mean(mean_sel_r))
#   experiment         phase       mean
# 1 Normal             Growth    0.640 
# 2 Normal             Settling  4.10  
# 3 Equal Growth Rates Growth   -0.0682
# 4 Equal Growth Rates Settling  4.06

mean_filt_sel_r_top %>%
  group_by(experiment, phase) %>%
  summarise(mean=mean(mean_sel_r),
            var=var(mean_sel_r),
            sd=sd(mean_sel_r),
            CV=sd(mean_sel_r)/mean(mean_sel_r),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

t.test(mean_sel_r~experiment, mean_filt_sel_r_top[mean_filt_sel_r_top$phase=='Settling',])
# t = 0.41021, df = 88.534, p-value = 0.6826




# Supplementary Figure 6 ####

mean_sim_sel_r=filt_combined_sel_r %>%
  group_by(phase, sim_number, experiment) %>%
  summarise(mean_sel_r=mean(sel_r))

mean_filt_sel_r %>%
  group_by(experiment, phase) %>%
  summarise(mean=mean(mean_sel_r))
#   experiment         phase       mean
# 1 Normal             Growth    0.642 
# 2 Normal             Settling  0.435 
# 3 Equal Growth Rates Growth   -0.0669
# 4 Equal Growth Rates Settling  0.484

fig5c_v2=ggplot(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay +15 min vs \nPetite) per selection phase')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 5.5, by = 1))+
  coord_cartesian(ylim = c(0, 5.5))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle('Random Starting Position')+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig5c_v2


fig5d_v2=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay +15 min vs \nPetite) per selection phase')+
  # ylab('')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_y_continuous(breaks = seq(0, 5.5, by = 1))+
  coord_cartesian(ylim = c(0, 5.5))+
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  ggtitle("Top Starting Position")+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=11))+
  NULL
fig5d_v2


supp_fig=plot_grid(fig5c_v2, fig5d_v2, labels=c('A', 'B'), ncol=1, align='hv')
supp_fig

# ggsave(filename="~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig8_selection_rates_equal_growth_23apr2025.png",
#        plot=supp_fig, dpi='retina', width=4, height=6)


# Figure 7 updated merged with supplement ####

mean_sim_sel_r %>%
  group_by(experiment, phase) %>%
  summarise(median=median(mean_sel_r))
#   experiment         phase     median
# 1 Normal             Growth    0.641 
# 2 Normal             Settling  0.470 
# 3 Equal Growth Rates Growth   -0.0671
# 4 Equal Growth Rates Settling  0.496

fig_a=ggplot(mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_color_manual(values = rev(ghibli_palette("MarnieMedium2"))) +
  scale_y_continuous(breaks = seq(0, 7, by = 1)) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  theme_classic(base_size = 10)+
  theme(
    plot.title = element_text(hjust = 0.5, size=10),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )+
  NULL
fig_a



wilcox.test(mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal',]$mean_sel_r~mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal',]$phase)
# p-value = 1.075e-10

cliff.delta(mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal',]$mean_sel_r~mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal',]$phase)
# delta estimate: 0.7496 (large)

fig_b=ggplot(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',], 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_color_manual(values = rev(ghibli_palette("MarnieMedium2"))) +
  scale_y_continuous(breaks = seq(0, 7, by = 1)) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  theme_classic(base_size = 10)+
  theme(
    plot.title = element_text(hjust = 0.5, size=10),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )+
  NULL
fig_b


wilcox.test(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',]$mean_sel_r~mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',]$phase)
# p-value < 2.2e-16

cliff.delta(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',]$mean_sel_r~mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',]$phase)
# delta estimate: -1 (large)


mean_filt_sel_r_top %>%
  group_by(experiment, phase) %>%
  summarise(median=median(mean_sel_r))
#   experiment         phase     median
# 1 Normal             Growth    0.640 
# 2 Normal             Settling  4.05  
# 3 Equal Growth Rates Growth   -0.0681
# 4 Equal Growth Rates Settling  4.14

fig_c=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate')+
  # ylab('')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_color_manual(values = rev(ghibli_palette("MarnieMedium2"))) +
  scale_y_continuous(breaks = seq(0, 7, by = 1)) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5, size=10))+
  NULL
fig_c


wilcox.test(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal',]$mean_sel_r~mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal',]$phase)
# p-value < 2.2e-16

cliff.delta(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal',]$mean_sel_r~mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal',]$phase)
# delta estimate: -1 (large)


fig_d=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',], 
                aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate')+
  # ylab('')+
  xlab('Selection phase')+
  guides(col='none')+
  scale_color_manual(values = rev(ghibli_palette("MarnieMedium2"))) +
  scale_y_continuous(breaks = seq(0, 7, by = 1)) +
  coord_cartesian(ylim = c(0, 7)) +
  scale_x_discrete(labels=c('Growth'='Growth',
                            'Settling'='Settling'))+
  theme_classic(base_size = 10)+
  theme(
    plot.title = element_text(hjust = 0.5, size=10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )+
  NULL
fig_d


wilcox.test(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',]$mean_sel_r~mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',]$phase)
# p-value < 2.2e-16

cliff.delta(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',]$mean_sel_r~mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',]$phase)
# delta estimate: -1 (large)

#Comparison settling top starting standard and equal growth rates
## Top settling selection
wilcox.test(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates' & mean_filt_sel_r_top$phase=='Settling',]$mean_sel_r, mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal' & mean_filt_sel_r_top$phase=='Settling',]$mean_sel_r)
# p-value = 0.9862
cliff.delta(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates' & mean_filt_sel_r_top$phase=='Settling',]$mean_sel_r, mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Normal' & mean_filt_sel_r_top$phase=='Settling',]$mean_sel_r)
# delta estimate: -0.0024 (negligible)

## Random position settling
wilcox.test(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates' & mean_sim_sel_r$phase=='Settling',]$mean_sel_r, mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal' & mean_sim_sel_r$phase=='Settling',]$mean_sel_r)
# p-value = 0.1845
cliff.delta(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates' & mean_sim_sel_r$phase=='Settling',]$mean_sel_r, mean_sim_sel_r[mean_sim_sel_r$experiment=='Normal' & mean_sim_sel_r$phase=='Settling',]$mean_sel_r)
# delta estimate: 0.1544 (small) 


figure_7=plot_grid(fig_a, fig_b, fig_c, fig_d, labels=c('A', 'B', 'C', 'D'), ncol=2, align='hv')
figure_7


# Add column labels at the top
column_labels <- plot_grid(
  NULL,
  ggdraw() + draw_label("Standard Growth Rates", size = 10, fontface="bold"),
  ggdraw() + draw_label("Equal Growth Rates", size = 10, fontface="bold"),
  NULL,
  ncol = 4,
  rel_widths = c(0.05, 1, 1, 0.05)
)

# Add row labels on the right
row_labels <- plot_grid(
  ggdraw() + draw_label("Random Starting Position", angle = 270, size = 10, fontface="bold"),
  ggdraw() + draw_label("Top Starting Position", angle = 270, size = 10, fontface="bold"),
  ncol = 1
)

# Combine everything
figure_7_annotated <- ggdraw(
  plot_grid(
    column_labels,
    plot_grid(
      figure_7,
      row_labels,
      ncol = 2,
      rel_widths = c(1, 0.05)
    ),
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

figure_7_annotated


ggsave(filename="~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_selection_rates_10mill_12mar2026.png",
       plot=figure_7_annotated, dpi='retina', width=6.5, height=5)





