

#Date: 16oct2024
# This code is used to analyze the results using al clusters from an exponentially growing tree
#using edge degree fragmentation, updated to use Petite, Petite w/o Delay and Grande strains

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
library(effsize)


theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))


# setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig5_frag_all-clusters_edge_degree")


# Fragmentation summary ####

# frag_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_15e_100g_all-clusters/fragmentation_inf.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   frag_df=rbind(frag_df, temp_df)
# }
# frag_df$strain <- ifelse(frag_df$strain=='petite', 'Petite', frag_df$strain)
# frag_df$strain <- ifelse(frag_df$strain=='petite-second-only', 'Petite w/o Delay', frag_df$strain)
# frag_df$strain <- ifelse(frag_df$strain=='grande', 'Grande', frag_df$strain)
# frag_df$strain=factor(frag_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# summary(frag_df)
# 
# table(frag_df$generation)
# 
# write.csv(frag_df, "growth_frag_all_frag_df_29may2025.csv", row.names=FALSE)


frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_4_network_growth_with_fragmentation/growth_frag_all_frag_df_29may2025.csv", header = TRUE)
frag_df$strain=factor(frag_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
summary(frag_df)

table(frag_df$generation)


#### Size at fracture ####

ggplot(frag_df,
       aes(x=as.factor(generation), y=cluster_size, fill=strain))+
  geom_violin()+
  ylab("Cluster size at fracture")+xlab("Generations")+
  petite_t200_colors+
  theme(legend.position="bottom")+
  NULL

ggplot(frag_df, aes(x=strain, y=cluster_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('Strain')+ylab('Mean Fracture Size')+
  petite_t200_colors+
  guides(fill='none')+
  NULL


tapply(frag_df$cluster_size, frag_df$strain, mean)
# edge degree 14:
#   Petite Petite w/o delay           Grande 
# 158.7926         248.3067         224.1929 
# edge degree 15:
#   Petite Petite w/o Delay           Grande 
# 218.1666         353.9265         340.0217

pairwise.t.test(frag_df$cluster_size, frag_df$strain, method="bonferroni")
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16          
# P value adjustment method: holm 

mean_clust_size=frag_df %>%
  group_by(strain, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))


ggplot(mean_clust_size, aes(x=generation, y=mean_size, col=strain))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=strain), 
              alpha=0.3, linetype='blank')+
  # ggtitle("Mean Fracture Size")+
  ylab("Mean fracure size")+xlab('Generation')+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL


#### undivided daughter cells ####


mean_undivided_cells=frag_df %>%
  group_by(strain, sim_number) %>%
  summarise(mean_undivided=mean(cases_mother_with_undivided_cells))

mean_undivided_cells %>%
  group_by(strain) %>%
  summarise(mean_mean_undivided=mean(mean_undivided))
# edge degree: 14
# strain           mean_mean_undivided
# 1 Petite                          5.93
# 2 Petite w/o delay                4.78
# 3 Grande                          3.00
# edge degree: 15
# strain           mean_mean_undivided
# 1 Petite                          8.02
# 2 Petite w/o Delay                6.78
# 3 Grande                          4.37


#### Fracture proportion ####

summary(frag_df)

mean_frag_prop=frag_df %>%
  group_by(strain, generation) %>%
  summarise(mean_frac=mean(proportion_propagule),
            sd_frac=sd(proportion_propagule),
            min_frac=min(proportion_propagule),
            max_frac=max(proportion_propagule))

ggplot(mean_frag_prop)+
  geom_line(aes(x=generation, y=mean_frac, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_frac, ymin=mean_frac-sd_frac, ymax=mean_frac+sd_frac, fill=strain), 
              alpha=0.2, linetype='blank')+
  # ggtitle('Mean Propagule Proportion')+
  ylab("Mean propagule proportion")+xlab('Generation')+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL


ggplot(frag_df, 
       aes(x=as.factor(generation), y=proportion_propagule, fill=strain))+
  geom_violin()+
  xlab("Generation")+ylab("Proportion propagule size")+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL

ggplot(frag_df, aes(x=strain, y=proportion_propagule, fill=strain))+
  geom_violin()+
  xlab('Generation')+
  ylab('Proportion Propagule Size')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

frag_df %>%
  group_by(strain) %>%
  summarise(mean=mean(proportion_propagule))
# edge degree: 14
# strain            mean
# 1 Petite           0.395
# 2 Petite w/o delay 0.392
# 3 Grande           0.453
# edge degree: 15
# strain            mean
# 1 Petite           0.393
# 2 Petite w/o Delay 0.394
# 3 Grande           0.454

pairwise.t.test(frag_df$proportion_propagule, frag_df$strain)
# edge degree: 14
#                  Petite  Petite w/o delay
# Petite w/o delay 7.9e-08 -               
# Grande           <2e-16  <2e-16          
# P value adjustment method: holm 
# edge degree: 15
#                Petite Petite w/o Delay
# Petite w/o Delay 0.73   -               
# Grande           <2e-16 <2e-16

#### Propagule size ####

mean_prop_size=frag_df %>%
  group_by(strain, generation) %>%
  summarise(mean_prop_size=mean(size_propagule))

ggplot(mean_prop_size, aes(x=strain, y=mean_prop_size, fill=strain))+
  geom_violin()+
  xlab('Generation')+
  ylab('Proportion Propagule Size')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

aggregate(mean_prop_size ~ strain, data = mean_prop_size, FUN = summary)

aggregate(size_propagule ~ strain, data = frag_df, FUN = summary)

frag_df %>% 
  group_by(strain) %>%
  summarise(mean = mean(size_propagule),
            lower_ci = mean - qt(0.975, n()-1) * (sd(size_propagule)/sqrt(n())),
            upper_ci = mean + qt(0.975, n()-1) * (sd(size_propagule)/sqrt(n())))
# edge degree: 15
# strain            mean lower_ci upper_ci
# 1 Petite            85.1     84.8     85.3
# 2 Petite w/o Delay 138.     137.     138. 
# 3 Grande           154.     154.     155. 

# Network Diameter ####

# diameter_df=data.frame()
# 
# for(j in c("petite", "grande", "petite-second-only")){
#   temp_df=read.csv(paste(j, "_15e_100g_all-clusters/networks_diameter.csv", sep=""), header=TRUE)
#   temp_df$strain=j
#   diameter_df=rbind(diameter_df, temp_df)
# }
# diameter_df$strain <- ifelse(diameter_df$strain=='petite', 'Petite', diameter_df$strain)
# diameter_df$strain <- ifelse(diameter_df$strain=='petite-second-only', 'Petite w/o Delay', diameter_df$strain)
# diameter_df$strain <- ifelse(diameter_df$strain=='grande', 'Grande', diameter_df$strain)
# diameter_df$strain=factor(diameter_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
# summary(diameter_df)
# 
# table(diameter_df$generation)
# 
# write.csv(diameter_df, "growth_frag_all_diameter_df_29may2025.csv", row.names=FALSE)

diameter_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_4_network_growth_with_fragmentation/growth_frag_all_diameter_df_29may2025.csv", header = TRUE)

diameter_df$strain=factor(diameter_df$strain, levels=c('Petite', 'Petite w/o Delay', 'Grande'))
summary(diameter_df)
summary(frag_df)

#To normalize the network diameter metric I need to first join the diameter and size at fracture information
#to later divide network diameter by cluster size

diam_size <- left_join(diameter_df, frag_df, 
                       by = c("sim_number", "generation", "strain", "cluster_id"))

summary(diam_size)

# diam_size$norm_diameter=diam_size$diameter/diam_size$cluster_size
diam_size$norm_diameter=diam_size$diameter/(6.64385619*log10(diam_size$cluster_size)-1)
summary(diam_size)

#normalizing diameter
ggplot(diam_size, aes(x=strain, y=norm_diameter, fill=strain))+
  geom_violin()+
  xlab('Strain')+ylab('Normalized Diameter')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  guides(fill='none')+
  NULL

diam_size %>%
  group_by(strain) %>%
  summarise(mean=mean(norm_diameter, na.rm=TRUE))
# edge degree: 14
# strain              mean
# <fct>              <dbl>
# 1 Petite           0.910
# 2 Petite w/o delay 0.995
# 3 Grande           0.988
# edge degree: 15
# strain            mean
# 1 Petite           0.913
# 2 Petite w/o Delay 0.997
# 3 Grande           0.990

pairwise.t.test(diam_size$norm_diameter, diam_size$strain)
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16


# Plots of diameter by generation
summ_diam=diam_size %>%
  group_by(strain, generation) %>%
  summarise(mean_diam=mean(norm_diameter),
            sd_diam=sd(norm_diameter),
            median_diam=median(norm_diameter))

mean_summ_diam=summ_diam %>%
  group_by(strain, generation) %>%
  summarise(across(mean_diam, list(mean = mean, sd = sd), .names = "{.fn}_diam"))

ggplot(summ_diam)+
  geom_line(aes(x=generation, y=mean_diam, col=strain))+
  geom_ribbon(aes(x=generation, y=mean_diam, ymin=mean_diam-sd_diam, ymax=mean_diam+sd_diam, fill=strain), alpha=0.2)+
  # ggtitle("Mean of Normalized Network Diameter")+
  xlab("Generation")+
  ylab("Normalized Diameter")+
  petite_t200_colors+
  theme(legend.position='bottom')+
  NULL



# Paper Figure 4 ####


# Mean of means


frag_df %>% 
  group_by(strain) %>%
  summarise(mean_frac=mean(cluster_size),
            mean_prop=mean(proportion_propagule))
# edge degree: 14
# strain                 mean_frac mean_prop
# 1 Petite                158.7926     0.395
# 2 Petite w/o delay      248.3067     0.392
# 3 Grande                224.1929     0.453
# edge degree: 15
#   strain           mean_frac mean_prop
# 1 Petite                218.     0.393
# 2 Petite w/o Delay      354.     0.394
# 3 Grande                340.     0.454

mean_frag_sim=frag_df %>% 
  group_by(strain, sim_number) %>%
  summarise(mean_clust_size=mean(cluster_size),
            mean_propagule_prop=mean(proportion_propagule))
summary(mean_frag_sim)

mean_frag_sim %>% 
  group_by(strain) %>%
  summarise(mean_frac=mean(mean_clust_size),
            mean_prop=mean(mean_propagule_prop))
# edge degree: 14
# strain                mean_frac mean_prop
# 1 Petite            158.7926    0.395
# 2 Petite w/o delay  248.3067    0.392
# 3 Grande            224.1929    0.453
# The means are the same as if all the data was used at the same time
# edge degree: 15
#   strain              mean_frac  mean_prop
# 1 Petite                218.1666     0.393
# 2 Petite w/o Delay      353.9265     0.394
# 3 Grande                340.0217     0.454

anova_frac_sim=aov(mean_clust_size~strain, data=mean_frag_sim)
summary(anova_frac_sim)
# edge degree 14
#              Df Sum Sq Mean Sq F value Pr(>F)    
# strain        2 429048  214524   19828 <2e-16 ***
# Residuals   297   3213      11  
# edge degree 15
#              Df  Sum Sq Mean Sq F value Pr(>F)    
# strain        2 1115758  557879   21935 <2e-16 ***
# Residuals   297    7554      25 

pairwise.t.test(mean_frag_sim$mean_clust_size, mean_frag_sim$strain, p.adjust.method='bonferroni')
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16


anova_propagule_sim=aov(mean_propagule_prop~strain, data=mean_frag_sim)
summary(anova_propagule_sim)
# edge degree 14
#              Df  Sum Sq Mean Sq F value Pr(>F)    
# strain        2 0.23223 0.11612    8858 <2e-16 ***
# Residuals   297 0.00389 0.00001  
# edge degree 15
#              Df  Sum Sq Mean Sq F value Pr(>F)    
# strain        2 0.24496 0.12248    8992 <2e-16 ***
# Residuals   297 0.00405 0.00001     

pairwise.t.test(mean_frag_sim$mean_propagule_prop, mean_frag_sim$strain, p.adjust.method='bonferroni')
# edge degree 14
#                  Petite  Petite w/o delay
# Petite w/o delay 5.3e-06 -               
# Grande           < 2e-16 < 2e-16  
# edge degree 15
#                  Petite Petite w/o Delay
# Petite w/o Delay 1      -               
# Grande           <2e-16 <2e-16 



mean_frag_sim %>% 
  group_by(strain) %>%
  summarise(median_frac=median(mean_clust_size),
            median_prop=median(mean_propagule_prop))
#   strain           median_frac median_prop
# 1 Petite                  218.       0.393
# 2 Petite w/o Delay        354.       0.394
# 3 Grande                  339.       0.454


pairwise.wilcox.test(mean_frag_sim$mean_clust_size, mean_frag_sim$strain, p.adjust.method='bonferroni')
#                  Petite Petite w/o Delay
# Petite w/o Delay <2e-16 -               
# Grande           <2e-16 <2e-16 

pairwise.wilcox.test(mean_frag_sim$mean_propagule_prop, mean_frag_sim$strain, p.adjust.method='bonferroni')
#                  Petite Petite w/o Delay
# Petite w/o Delay 1      -               
# Grande           <2e-16 <2e-16 


p1_temp=ggplot(mean_frag_sim, aes(x=strain, y=mean_clust_size, fill=strain))+
  geom_violin()+
  stat_summary(fun='mean', geom='crossbar')+
  xlab('')+ylab('Mean Fracture\nSize')+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  guides(fill='none')+
  NULL
p1_temp

p2_temp=ggplot(mean_frag_sim, aes(x=strain, y=mean_propagule_prop, fill=strain))+
  geom_violin()+
  xlab('')+
  ylab('Mean Propagule\nProportion')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  guides(fill='none')+
  NULL
p2_temp


mean_diam_sim=diam_size %>% 
  group_by(strain, sim_number) %>%
  summarise(mean_norm_diam=mean(norm_diameter))
summary(mean_diam_sim)

mean_diam_sim %>% 
  group_by(strain) %>%
  summarise(mean_diam=mean(mean_norm_diam))
# edge degree 14
# strain             mean_diam
# 1 Ancestor               0.910
# 2 Ancestor w/o delay     0.995
# 3 PA2_t200               0.988
# edge degree 15
# strain           mean_diam
# 1 Petite               0.913
# 2 Petite w/o Delay     0.997
# 3 Grande               0.990

anova_diam_sim=aov(mean_norm_diam~strain, data=mean_diam_sim)
summary(anova_diam_sim)
# edge degree 14
#              Df Sum Sq Mean Sq F value Pr(>F)    
# strain        2 0.4422 0.22112   22264 <2e-16 ***
# Residuals   297 0.0029 0.00001 
# edge degree 15
#              Df Sum Sq Mean Sq F value Pr(>F)    
# strain        2 0.4394 0.21969   23051 <2e-16 ***
# Residuals   297 0.0028 0.00001    

pairwise.t.test(mean_diam_sim$mean_norm_diam, mean_diam_sim$strain, p.adjust.method='bonferroni')
#                  Petite Petite w/o delay
# Petite w/o delay <2e-16 -               
# Grande           <2e-16 <2e-16 


mean_diam_sim %>% 
  group_by(strain) %>%
  summarise(median_diam=median(mean_norm_diam))
#   strain           median_diam
# 1 Petite                 0.912
# 2 Petite w/o Delay       0.997
# 3 Grande                 0.990

pairwise.wilcox.test(mean_diam_sim$mean_norm_diam, mean_diam_sim$strain, p.adjust.method='bonferroni')
#                  Petite Petite w/o Delay
# Petite w/o Delay <2e-16 -               
# Grande           <2e-16 <2e-16  

p3_temp=ggplot(mean_diam_sim, aes(x=strain, y=mean_norm_diam, fill=strain))+
  geom_violin()+
  xlab('')+ylab('Mean Normalized\nDiameter')+
  stat_summary(fun='mean', geom='crossbar')+
  petite_t200_colors+
  theme_classic(base_size = 10)+
  guides(fill='none')+
  NULL
p3_temp





fig_restructured_updated=plot_grid(p1_temp, p3_temp, p2_temp,
                               labels=c('A', 'B', 'C'), ncol=1, label_size=11)
fig_restructured_updated

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/fig_4_network_properties_all-clusters_28aug2026.png',
       plot=fig_restructured_updated, dpi='retina', width=3.5, height=6.3, bg='white')
