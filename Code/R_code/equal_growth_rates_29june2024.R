
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

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))


#### Growth phase ####

growth_prop=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/settling_simulation_data/set_sim_growth_registry_equal_growth_rate_15apr2024.csv", header=TRUE)

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

extinction_inf=growth_prop %>%
  group_by(sim_number, experiment) %>%
  summarize(total_transfers=max(transfer))

extinction_inf %>%
  group_by(experiment) %>%
  summarize(mean_transfers=mean(total_transfers))

ggplot(extinction_inf[extinction_inf$experiment=="Normal",], 
       aes(x=total_transfers, col=experiment, fill=experiment))+
  geom_histogram(position="identity", alpha=0.5, bins=10)+
  theme_classic(base_size=16)+
  guides(fill='none', col='none')+
  labs(x="Transfers until Fixation", y="Count")+
  # facet_wrap(~experiment)+
  geom_vline(xintercept=6.62, linetype='dashed', size=1.2)+
  NULL

#### Settling fitness ####
settling_prop=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/settling_simulation_data/set_sim_settling_registry_equal_growth_rate_15apr2024.csv", header=TRUE)

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
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
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
  scale_x_continuous(breaks = seq(1, 20, by = 1))+
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
# experiment          mean   var    sd    CV
# 1 Normal             0.382 0.594 0.771  2.02
# 2 Equal Growth Rates 0.460 0.429 0.655  1.42

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
# experiment            mean
# 1 Normal              0.644 
# 2 Equal Growth Rates -0.0670

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
# experiment          mean
# 1 Normal             0.450
# 2 Equal Growth Rates 0.493


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
# experiment         phase       mean
# 1 Normal             Growth    0.644 
# 2 Normal             Settling  0.450 
# 3 Equal Growth Rates Growth   -0.0670
# 4 Equal Growth Rates Settling  0.493 

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

t.test(mean_sel_r~experiment, mean_filt_sel_r[mean_filt_sel_r$phase=='Settling',])


# Top Settling selection ####
#This results are doing settling selection starting all the clusters at the top instead of
#them starting at a random position. 


#### Growth phase ####

growth_prop_top=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/settling_simulation_data/set_sim_top_growth_registry_equal_growth_rate_15apr2024.csv", header=TRUE)

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

#### Settling fitness ####

settling_prop_top=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/settling_simulation_data/set_sim_top_settling_registry_equal_growth_rate_15apr2024.csv", header=TRUE)

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
# experiment         phase       mean
# 1 Normal             Growth    0.639 
# 2 Normal             Settling  4.28  
# 3 Equal Growth Rates Growth   -0.0664
# 4 Equal Growth Rates Settling  3.80 

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
# p-value = 0.01578
# The mean between using the unmodified doubling time and increasing the mean for the Ancestor w/o
# Delay make the settling selection rate have significantly different means (p-value = 0.01578)



# Supplementary Figure 6 ####

mean_sim_sel_r=filt_combined_sel_r %>%
  group_by(phase, sim_number, experiment) %>%
  summarise(mean_sel_r=mean(sel_r))

fig5c_v2=ggplot(mean_sim_sel_r[mean_sim_sel_r$experiment=='Equal Growth Rates',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay +15 min vs \nPetite) per selection phase')+
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
fig5c_v2


fig5d_v2=ggplot(mean_filt_sel_r_top[mean_filt_sel_r_top$experiment=='Equal Growth Rates',], 
             aes(x=phase, y=mean_sel_r, col=phase))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun = 'mean', geom = 'crossbar', linetype = 'dashed', alpha = 0.5, col = 'black')+
  ylab('Mean Selection Rate \n(Petite w/o delay +15 min vs \nPetite) per selection phase')+
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
fig5d_v2


supp_fig=plot_grid(fig5c_v2, fig5d_v2, labels=c('A', 'B'), ncol=1, align='hv')
supp_fig

ggsave(filename="~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig8_selection_rates_equal_growth_23apr2025.png",
       plot=supp_fig, dpi='retina', width=4, height=6)






