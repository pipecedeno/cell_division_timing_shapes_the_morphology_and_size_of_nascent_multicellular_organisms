

# Date 25apr2023
# This code was made to analyze the results of using synthetic data to analyze the effects
# of variation and delay in the doubling time distributions regarding the cluster properties
# The results were made with the python script sim_frag_clust_diff_mean_dists_24apr2024.py


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
library(ghibli)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$YesterdayMedium)
syn_data_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("YesterdayMedium", direction = -1))

default_mean=105
std_vector=c(0, 17.33, 34.66, 52)
delay_vector=seq(0, 105, length.out=4)
delay_percentages=c(0, 33, 66, 100)


color_diff_mean=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
                                          breaks = delay_vector, limits = c(min(delay_vector), max(delay_vector))))
color_diff_mean_perc=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
                                          breaks = delay_percentages, limits = c(min(delay_percentages), max(delay_percentages))))
color_variation=list(scale_color_gradient( name = "Variation", low = "pink", high = "darkred", 
                                           breaks = std_vector, limits = c(min(std_vector), max(std_vector))))

setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/corrections_supreme_version_sep2026/fig5_effects_delay_variation_14sep2026")


# Fragmentation summary ####

# frag_df=data.frame()
# 
# for (temp_var in std_vector){
#   for (temp_diff in delay_vector){
#     temp_df=read.csv(paste("test_", temp_var,"_var_", temp_diff,"_diff/fragmentation_inf.csv", sep=""), header=TRUE)
#     temp_df$variation=temp_var
#     temp_df$diff_mean=temp_diff
#     frag_df=rbind(frag_df, temp_df)
#   }
# }
# 
# summary(frag_df)
# 
# write.csv(frag_df, "frag_syn_frag_df_15sep2026.csv", row.names = FALSE)

frag_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_5_delay_variation_in_cluster_properties/frag_syn_frag_df_15sep2026.csv",
                 header=TRUE)

summary(frag_df)


table(frag_df$generation)


#### Size at fracture ####

mean_clust_size=frag_df %>%
  group_by(variation, diff_mean, generation) %>%
  summarise(mean_size=mean(cluster_size),
            sd_size=sd(cluster_size),
            min_size=min(cluster_size),
            max_size=max(cluster_size))


ggplot(mean_clust_size, aes(x=generation, y=mean_size, col=variation))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_size, ymin=mean_size-sd_size, ymax=mean_size+sd_size, fill=variation), 
              alpha=0.2, linetype='blank')+
  ggtitle("Mean Fracture Size")+
  ylab("Mean fracure size")+xlab('Generation')+
  facet_grid(variation~diff_mean)+
  # syn_data_colors+
  NULL


mean_mean_clust_size <- mean_clust_size %>%
  group_by(variation, diff_mean) %>%
  summarize(mean = mean(mean_size),
            sd = sd(mean_size),
            n = n(),
            se = sd / sqrt(n),
            lower_ci = mean - qt(0.975, df = n - 1) * se,
            upper_ci = mean + qt(0.975, df = n - 1) * se)

ggplot(mean_mean_clust_size, aes(x=variation, y=mean, col=diff_mean, group=diff_mean)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  xlab('Variation') +
  ylab('Mean fracture size') +
  color_diff_mean+
  NULL


ggplot(mean_mean_clust_size, aes(x=diff_mean, y=mean, col=variation, group=variation)) +
  geom_point() +
  geom_line() +
  xlab('Delay First Division') +
  ylab('Mean fracture size') +
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1) +
  color_variation+
  NULL

#### Fracture proportion ####

summary(frag_df)

mean_frag_prop=frag_df %>%
  group_by(variation, diff_mean, generation) %>%
  summarise(mean_propagule=mean(proportion_propagule),
            sd_propagule=sd(proportion_propagule))


ggplot(mean_frag_prop)+
  theme_bw()+
  geom_line(aes(x=generation, y=mean_propagule, col=variation))+
  geom_ribbon(aes(x=generation, y=mean_propagule, ymin=mean_propagule-sd_propagule, ymax=mean_propagule+sd_propagule, fill=variation), 
              alpha=0.2, linetype='blank')+
  facet_grid(variation~diff_mean)+
  ggtitle('Propagule proportion after fracture')+
  ylab("Mean fracture proportion")+xlab('Generation')+
  NULL


ggplot(frag_df[frag_df$generation%%10==0,], 
       aes(x=as.factor(generation), y=proportion_propagule, fill=variation))+
  geom_violin()+
  theme_bw()+
  xlab("Generation")+ylab("Proportion offspring size")+
  facet_grid(variation~diff_mean)+
  # syn_data_colors+
  NULL


mean_mean_frag_prop=mean_frag_prop %>%
  group_by(variation, diff_mean) %>%
  summarize(mean_prop = mean(mean_propagule),
            sd_prop = sd(mean_propagule),
            n = n(),
            se = sd_prop / sqrt(n),
            lower_ci = mean_prop - qt(0.975, df = n - 1) * se,
            upper_ci = mean_prop + qt(0.975, df = n - 1) * se)

ggplot(mean_mean_frag_prop, aes(x=diff_mean, y=mean_prop, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay First Division')+
  ylab('Mean propagule proportion')+
  guides(col = guide_legend(title = "Variation"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  NULL

ggplot(mean_mean_frag_prop, aes(x=variation, y=mean_prop, col=diff_mean, group=diff_mean))+
  geom_point()+
  geom_line()+
  xlab('Variation')+
  ylab('Mean propagule proportion')+
  guides(col = guide_legend(title = "Delay"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean+
  NULL


# Network Diameter ####

# diameter_df=data.frame()
# 
# for (temp_var in std_vector){
#   for (temp_diff in delay_vector){
#     temp_df=read.csv(paste("test_", temp_var,"_var_", temp_diff,"_diff/networks_diameter.csv", sep=""), header=TRUE)
#     temp_df$variation=temp_var
#     temp_df$diff_mean=temp_diff
#     diameter_df=rbind(diameter_df, temp_df)
#   }
# }
# 
# summary(diameter_df)
# 
# write.csv(diameter_df, "frag_syn_diameter_df_15sep2026.csv", row.names = FALSE)


diameter_df=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/fig_5_delay_variation_in_cluster_properties/frag_syn_diameter_df_15sep2026.csv", header=TRUE)

summary(diameter_df)


summary(frag_df)

#To normalize the network diameter metric I need to first join the diameter and size at fracture information
#to later divide network diameter by cluster size

diam_size <- left_join(frag_df, diameter_df, 
                       by = c("sim_number", "generation", "variation", 'diff_mean'))

summary(diam_size)

# diam_size$norm_diameter=diam_size$diameter/diam_size$cluster_size
diam_size$norm_diameter=diam_size$diameter/(6.64385619*log10(diam_size$cluster_size)-1)


# Plots of diameter by generation
summ_diam=diam_size %>%
  group_by(variation, diff_mean, generation) %>%
  summarise(mean_diam=mean(norm_diameter),
            sd_diam=sd(norm_diameter),
            median_diam=median(norm_diameter))



ggplot(summ_diam, aes(x=generation, y=mean_diam, col=variation))+
  geom_line()+
  geom_ribbon(aes(x=generation, y=mean_diam, ymin=mean_diam-sd_diam, ymax=mean_diam+sd_diam, fill=variation), 
              alpha=0.2, linetype='blank')+
  ggtitle("Mean of Normalized Network Diameter")+
  xlab("Generation")+
  ylab("Normalized network diameter")+
  facet_grid(variation~diff_mean)+
  NULL

mean_summ_diam=summ_diam %>%
  group_by(variation, diff_mean) %>%
  summarize(mean_d = mean(mean_diam),
            sd_d = sd(mean_diam),
            n = n(),
            se = sd_d / sqrt(n),
            lower_ci = mean_d - qt(0.975, df = n - 1) * se,
            upper_ci = mean_d + qt(0.975, df = n - 1) * se)

ggplot(mean_summ_diam, aes(x=diff_mean, y=mean_d, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay First Division')+
  ylab('Mean Diameter')+
  guides(col = guide_legend(title = "Variation"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  NULL

ggplot(mean_summ_diam, aes(x=variation, y=mean_d, col=diff_mean, group=diff_mean))+
  geom_point()+
  geom_line()+
  xlab('Variation')+
  ylab('Mean Diameter')+
  guides(col = guide_legend(title = "Delay"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean+
  NULL




#### Figure 5 ####


# Fracture size
clust_var=ggplot(mean_mean_clust_size, aes(x=variation, y=mean, col=floor(diff_mean/default_mean*100), group=floor(diff_mean/default_mean*100)))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  xlab('Variation')+
  ylab('Mean Fracture\nSize (# nodes)')+
  # guides(col = guide_legend(title = "Delay"))+
  guides(col='none')+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
clust_var

clust_delay=ggplot(mean_mean_clust_size, aes(x=diff_mean/default_mean*100, y=mean, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean Fracture\nSize (# nodes)')+
  # guides(col = guide_legend(title = "Variation"))+
  guides(col='none')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
clust_delay


# Fracture Proportion
prop_delay=ggplot(mean_mean_frag_prop, aes(x=diff_mean/default_mean*100, y=mean_prop, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean propagule\nproportion')+
  # guides(col = guide_legend(title = "Variation"))+
  guides(col='none')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
prop_delay

prop_var=ggplot(mean_mean_frag_prop, aes(x=variation, y=mean_prop, col=floor(diff_mean/default_mean*100), group=floor(diff_mean/default_mean*100)))+
  geom_point()+
  geom_line()+
  xlab('Variation (SD)')+
  ylab('Mean propagule\nproportion')+
  # guides(col = guide_legend(title = "Delay"))+
  guides(col='none')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
prop_var


#Network diameter
diam_delay=ggplot(mean_summ_diam, aes(x=diff_mean/default_mean*100, y=mean_d, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean Normalized\nNetwork Diameter')+
  # guides(col = guide_legend(title = "Variation"))+
  guides(col='none')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
diam_delay

diam_var=ggplot(mean_summ_diam, aes(x=variation, y=mean_d, col=floor(diff_mean/default_mean*100), group=floor(diff_mean/default_mean*100)))+
  geom_point()+
  geom_line()+
  xlab('Variation')+
  ylab('Mean Network\nDiameter')+
  # guides(col = guide_legend(title = "Delay"))+
  guides(col='none')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )+
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 15))+
  NULL
diam_var



common_legend_delay <- cowplot::get_plot_component(
  ggplot(mean_summ_diam, aes(x = variation, y = mean_d, 
                             col = floor(diff_mean/default_mean*100), 
                             group = floor(diff_mean/default_mean*100))) +
    geom_point() +
    geom_line() +
    xlab('Variation') +
    ylab('Mean Network\nDiameter') +
    guides(col = guide_legend(title = "Delay", direction = "horizontal")) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
    color_diff_mean_perc +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom"),
  pattern = "guide-box-bottom"
)

common_legend_variation <- cowplot::get_plot_component(
  ggplot(mean_summ_diam, aes(x = diff_mean/default_mean*100, y = mean_d, 
                             col = variation, group = variation)) +
    geom_point() +
    geom_line() +
    xlab('Delay (% Second Division)') +
    ylab('Mean Network\nDiameter') +
    guides(col = guide_legend(title = "Variation", direction = "horizontal")) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.1) +
    color_variation +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom"),
  pattern = "guide-box-bottom"
)




supp_net_prop_v4=plot_grid(clust_delay, clust_var, 
                           diam_delay, diam_var, 
                           prop_delay, prop_var, 
                           common_legend_variation, common_legend_delay,
                           labels=c('A', 'B', 'C', 'D', 'E', 'F'), ncol=2, 
                          label_size=10, rel_heights = c(1,1,1,0.1), rel_widths = c(1, 0.82))
supp_net_prop_v4

supp_net_prop_v4_test=plot_grid(common_legend_variation, common_legend_delay,
                           clust_delay, clust_var, 
                           diam_delay, diam_var, 
                           prop_delay, prop_var, 
                           labels=c('','','A', 'B', 'C', 'D', 'E', 'F'), ncol=2, 
                           label_size=10, rel_heights = c(0.1,1,1,1.17), rel_widths = c(1, 0.82))
supp_net_prop_v4_test

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/fig_5_effects_delay_variation_26oct2026.png',
       plot=supp_net_prop_v4_test, dpi='retina', width=5.6, height=5.6, bg='white')



