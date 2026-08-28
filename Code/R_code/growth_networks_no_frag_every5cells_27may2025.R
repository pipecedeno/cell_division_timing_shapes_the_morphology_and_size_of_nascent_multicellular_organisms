


#Date: 30may2025
# This code is used to analyze the results of growing the networks without fragmentation
# so it is only growing them until 1200 nodes in size, now the properties are being tracked
# every 5 cell divisions

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

library(slider)
library(effectsize)
library(ggbeeswarm)
library(colorspace)


theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig4_growth_no_frag_18june2025/")

# Loading the data ####

# sim_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_300n_1300m/network_information.csv", sep=""), header=TRUE)
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




# Network Diameter ####

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

# Saving data for supplementary figure 5 of network normalization
# write.csv(summ_diameter, "~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_4_network_growth_with_fragmentation/mean_diameter_values_3oct2025.csv", row.names = FALSE)


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




### Functional ANOVA


# Paper Figure ####


# img_petite <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_diameter_undivided_23apr2025.png")
# img_plot_petite_net <- rasterGrob(img_petite, interpolate = TRUE)
# 
# img_grande <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_diameter_undivided_23apr2025.png")
# img_plot_grande_net <- rasterGrob(img_grande, interpolate = TRUE)
# 
# # Create text annotations
# text_petite <- textGrob("Petite", gp = gpar(fontsize = 14, fontface = "bold"))
# text_grande <- textGrob("Grande", gp = gpar(fontsize = 14, fontface = "bold"))
# 
# # Create ggplot objects for the images with annotations
# p_petite <- ggplot() + 
#   annotation_custom(img_plot_petite_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   annotation_custom(text_petite, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
#   theme_void() +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p_petite
# 
# p_grande <- ggplot() + 
#   annotation_custom(img_plot_grande_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
#   annotation_custom(text_grande, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
#   theme_void() +
#   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p_grande



doubled=read.csv('~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/first_second_div_doubling_data_2024dec2.csv', header=TRUE)
doubled$division_number=factor(doubled$division_number)
doubled$timepoint=factor(doubled$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

petite_dt=doubled[doubled$strain=='petite' & as.numeric(doubled$division_number)<=2,]
petite_dt$name='Petite'

only_second_doubling=petite_dt[petite_dt$division_number==2,]
temp_first=only_second_doubling
temp_first$division_number=1
petite_second_doubling=rbind(only_second_doubling, temp_first)
petite_second_doubling$strain='second_doub'
petite_second_doubling$name='Petite w/o Delay'

petite_synthetic_data=rbind(petite_dt, petite_second_doubling)
summary(petite_synthetic_data)

grande_dt_temp=doubled[doubled$strain=='grande' & as.numeric(doubled$division_number)<=2,]
grande_dt_temp$name='Grande'

colnames(grande_dt_temp)

doub_t_data=rbind(petite_synthetic_data, grande_dt_temp)
doub_t_data$strain=factor(doub_t_data$strain, levels=c("petite", "second_doub", "grande"))
doub_t_data$name=factor(doub_t_data$name, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
table(doub_t_data$strain)

facet.labs=c('Ancestor', 'Ancestor \nw/o Delay', 'Evolved')
names(facet.labs)=c('Ancestor', 'Ancestor w/o Delay', 't200')




your_colors <- c("#F0D77B", "#B4DAE5", "#AE93BE")
darker_colors <- darken(your_colors, amount = 0.3)

dt_distributions_v2=ggplot(doub_t_data, aes(x=division_number, y=hours, fill=name))+
  facet_wrap(~name)+
  geom_violin(adjust=2)+
  # geom_beeswarm(cex=0.5)+
  # geom_jitter(alpha=0.5, size=0.25)+
  stat_summary(fun='median', geom='crossbar')+
  geom_quasirandom(aes(col=name), method="tukeyDense", size=0.1)+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none', color='none')+
  petite_t200_colors+
  scale_color_manual(values=darker_colors)+
  theme_classic(base_size = 10)+
  NULL
dt_distributions_v2


p_diam = ggplot(summ_diameter, aes(x=num_nodes, y=mean_diameter, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  geom_point(data = data.frame(
    num_nodes = 150,
    mean_diameter = summ_diameter[summ_diameter$strain=='Grande' & summ_diameter$num_nodes==150,]$mean_diameter
  ), 
    aes(x = num_nodes, y = mean_diameter), color = "#AE93BEFF", shape = 15, inherit.aes = FALSE, size = 1) +
  geom_point(data = data.frame(
    num_nodes = 150,
    mean_diameter = summ_diameter[summ_diameter$strain=='Petite' & summ_diameter$num_nodes==150,]$mean_diameter
  ), 
  aes(x = num_nodes, y = mean_diameter), color = "#F0D77BFF", shape = 15, inherit.aes = FALSE, size = 1) +
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Diameter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_diam
# Grande, mean: 13.33, observed: 13
# Petite, mean: 12.22667, observed: 12

# Petite vs Petite w/o Delay
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$diameter ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# p-value < 2.2e-16

cliff.delta(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$diameter ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# delta estimate: -0.6945111 (large)

# Petite w/o Delay vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$diameter ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$strain)
# p-value = 0.09466

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite', ])
cliff.delta(sub_df$diameter ~ sub_df$strain)
# delta estimate: 0.07146667 (negligible)

# Petite vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$diameter ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$strain)
# p-value < 2.2e-16

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite w/o Delay', ])
cliff.delta(sub_df$diameter ~ sub_df$strain)
# delta estimate: -0.6767111 (large)


p_undivided=ggplot(summ_undivided, aes(x=num_nodes, y=mean_undivided, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  geom_point(data = data.frame(
    num_nodes = 100,
    mean_diameter = summ_undivided[summ_undivided$strain=='Grande' & summ_undivided$num_nodes==100,]$mean_undivided
  ), 
  aes(x = num_nodes, y = mean_diameter), color = "#AE93BEFF", shape = 15, inherit.aes = FALSE, size = 1) +
  geom_point(data = data.frame(
    num_nodes = 100,
    mean_diameter = summ_undivided[summ_undivided$strain=='Petite' & summ_undivided$num_nodes==100,]$mean_undivided
  ), 
  aes(x = num_nodes, y = mean_diameter), color = "#F0D77BFF", shape = 15, inherit.aes = FALSE, size = 1) +
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  scale_y_continuous(breaks = seq(0, 110, 30)) +
  labs(y = "Mean Mothers with\nDelayed Daughter", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_undivided
# Grande, mean: 3.106667, observed: 3
# Petite, mean 8.806667, observed: 12


# Petite vs Petite w/o Delay
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$cases_mother_with_undivided_cells ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# p-value < 2.2e-16

cliff.delta(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$cases_mother_with_undivided_cells ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# delta estimate: 0.7811778 (large)

# Petite w/o Delay vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$cases_mother_with_undivided_cells ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$strain)
# p-value < 2.2e-16

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite', ])
cliff.delta(sub_df$cases_mother_with_undivided_cells ~ sub_df$strain)
# delta estimate: 0.4687222 (medium)

# Petite vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$cases_mother_with_undivided_cells ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$strain)
# p-value < 2.2e-16

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite w/o Delay', ])
cliff.delta(sub_df$cases_mother_with_undivided_cells ~ sub_df$strain)
# delta estimate: 0.9339778 (large)

p_edge_degree=ggplot(summ_max_edge, aes(x=num_nodes, y=mean_max_edge, col=strain, fill=strain))+
  geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), alpha=0.5, color=NA)+
  geom_line()+
  geom_point(data = data.frame(
    num_nodes = 150,
    mean_diameter = summ_max_edge[summ_max_edge$strain=='Grande' & summ_max_edge$num_nodes==150,]$mean_max_edge
  ), 
  aes(x = num_nodes, y = mean_diameter), color = "#AE93BEFF", shape = 15, inherit.aes = FALSE, size = 1) +
  geom_point(data = data.frame(
    num_nodes = 150,
    mean_diameter = summ_max_edge[summ_max_edge$strain=='Petite' & summ_max_edge$num_nodes==150,]$mean_max_edge
  ), 
  aes(x = num_nodes, y = mean_diameter), color = "#F0D77BFF", shape = 15, inherit.aes = FALSE, size = 1) +
  petite_t200_colors+
  theme_classic(base_size = 10)+
  scale_x_continuous(breaks = seq(0, 1200, 200)) +
  labs(y = "Mean Max\nEdge Degree", x = "Number of Nodes") +
  guides(col='none', fill='none')+
  NULL
p_edge_degree
# Grande, mean: 12.38, observed: 12
# Petite, mean: 13.51667, observed: 14


# Petite vs Petite w/o Delay
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$max_edge_degree ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# p-value < 2.2e-16

cliff.delta(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$max_edge_degree ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Grande',]$strain)
# delta estimate: 0.6696 (large)

# Petite w/o Delay vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$max_edge_degree ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite',]$strain)
# p-value = 0.2401

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite', ])
cliff.delta(sub_df$max_edge_degree ~ sub_df$strain)
# delta estimate: -0.04982222 (negligible)

# Petite vs Grande
wilcox.test(sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$max_edge_degree ~ 
              sim_df[sim_df$num_nodes==100 & sim_df$strain!='Petite w/o Delay',]$strain)
# p-value < 2.2e-16

sub_df <- droplevels(sim_df[sim_df$num_nodes == 100 & sim_df$strain != 'Petite w/o Delay', ])
cliff.delta(sub_df$max_edge_degree ~ sub_df$strain)
# delta estimate: 0.6709778 (large)

fig=plot_grid(dt_distributions_v2, p_undivided, p_diam, p_edge_degree,
              labels=c('A', 'B', 'C', 'D'), ncol=2, label_size=12)
fig

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_3_growth_no_frag_1apr2025_every5cells_jitter.svg',
       plot=fig, dpi='retina', width=6.5, height=4.5, bg='white')


