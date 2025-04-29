

#Date: 18nov2024
# This code is used to analyze the results of growing the networks without fragmentation
# so it is only growing them until 200 nodes in size

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


theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

setwd("~/work_dir/observed_synchrony/paper_results_2dec2024/growth_without_fragmentation_3dec2024/")

# Network Diameter ####

diameter_df=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/fig_3_network_growth_no_fragmentation/growth_no_frag_diameter_3dec2024.csv", header = TRUE)
diameter_df$strain=factor(diameter_df$strain, levels=c('Petite', 'Petite w/o delay', 'Grande'))
summary(diameter_df)


ggplot(diameter_df, aes(x=strain, y=diameter, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  NULL

pairwise.t.test(diameter_df$diameter, diameter_df$strain, p.adjust.method = "bonferroni")
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 0.059           
# P value adjustment method: bonferroni 

tapply(diameter_df$diameter, diameter_df$strain, summary)
# $Petite
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.0    13.0    13.0    13.1    14.0    15.0 
# $`Petite w/o delay`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12.00   14.00   14.00   14.19   15.00   16.00 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   14.00   14.00   14.05   14.00   15.00 


ggplot(diameter_df, aes(x=strain, y=cases_mother_with_undivided_cells, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  NULL

tapply(diameter_df$cases_mother_with_undivided_cells, diameter_df$strain, summary)
# $Petite
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.0    15.0    18.0    17.5    20.0    27.0 
# $`Petite w/o delay`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.000   7.000   9.000   9.047  11.000  17.000 
# $Grande
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    4.00    5.00    5.32    7.00   13.00

pairwise.t.test(diameter_df$cases_mother_with_undivided_cells, diameter_df$strain, p.adjust.method="bonferroni")
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16          
# P value adjustment method: bonferroni 

# Paper Figure ####


img_petite <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_diameter_undivided_23apr2025.png")
img_plot_petite_net <- rasterGrob(img_petite, interpolate = TRUE)

img_grande <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_diameter_undivided_23apr2025.png")
img_plot_grande_net <- rasterGrob(img_grande, interpolate = TRUE)

# Create text annotations
text_petite <- textGrob("Petite", gp = gpar(fontsize = 14, fontface = "bold"))
text_grande <- textGrob("Grande", gp = gpar(fontsize = 14, fontface = "bold"))

# Create ggplot objects for the images with annotations
p_petite <- ggplot() + 
  annotation_custom(img_plot_petite_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_petite, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_petite

p_grande <- ggplot() + 
  annotation_custom(img_plot_grande_net, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotation_custom(text_grande, xmin = 0.5, xmax = 0.5, ymin = 1, ymax = 1) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p_grande



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
petite_second_doubling$name='Petite w/o delay'

petite_synthetic_data=rbind(petite_dt, petite_second_doubling)
summary(petite_synthetic_data)

grande_dt_temp=doubled[doubled$strain=='grande' & as.numeric(doubled$division_number)<=2,]
grande_dt_temp$name='Grande'

colnames(grande_dt_temp)

doub_t_data=rbind(petite_synthetic_data, grande_dt_temp)
doub_t_data$strain=factor(doub_t_data$strain, levels=c("petite", "second_doub", "grande"))
doub_t_data$name=factor(doub_t_data$name, levels=c('Petite', 'Petite w/o delay', 'Grande'))
table(doub_t_data$strain)

facet.labs=c('Ancestor', 'Ancestor \nw/o delay', 'Evolved')
names(facet.labs)=c('Ancestor', 'Ancestor w/o delay', 't200')



dt_distributions_v2=ggplot(doub_t_data, aes(x=division_number, y=hours, fill=name))+
  facet_wrap(~name)+
  geom_violin(adjust=2)+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none')+
  petite_t200_colors+
  NULL
dt_distributions_v2


p_diam_v2=ggplot(diameter_df, aes(x=strain, y=diameter, fill=strain))+
  geom_violin()+
  # stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Network Diameter (# edges)')+
  NULL
p_diam_v2

p_undivided_v2=ggplot(diameter_df, aes(x=strain, y=cases_mother_with_undivided_cells, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  labs(x='Strain', y='Mothers with\nundivided nodes')+
  NULL
p_undivided_v2

fig_restructured_alt_v2=plot_grid(dt_distributions_v2, p_petite, p_grande, p_undivided_v2, p_diam_v2,
                               labels=c('A', 'B', 'C', 'D', 'E'), ncol=3, label_size=16)
fig_restructured_alt_v2

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_3_growth_no_frag_23apr2025_ancestors.png',
       plot=fig_restructured_alt_v2, dpi='retina', width=15, height=8, bg='white')
