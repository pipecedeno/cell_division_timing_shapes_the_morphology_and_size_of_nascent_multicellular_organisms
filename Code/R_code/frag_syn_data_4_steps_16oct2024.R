

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


std_vector=seq(0, 30, 10)
delay_vector=seq(0, 60, 20)
delay_percentages=c(0, 33, 66, 100)


color_diff_mean=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
                                          breaks = delay_vector, limits = c(min(delay_vector), max(delay_vector))))
color_diff_mean_perc=list(scale_color_gradient(name = "Delay", low = "lightblue", high = "darkblue", 
                                          breaks = delay_percentages, limits = c(min(delay_percentages), max(delay_percentages))))
color_variation=list(scale_color_gradient( name = "Variation", low = "pink", high = "darkred", 
                                           breaks = std_vector, limits = c(min(std_vector), max(std_vector))))

# setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig6_effects_delay_and_variation")
setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/results_edge_degree_15/fig6_effects_delay_and_variation_more_variation")


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
# write.csv(frag_df, "frag_syn_frag_df_24sep2025.csv", row.names = FALSE)

# frag_df=read.csv("frag_syn_frag_df_2june2025.csv", header=TRUE)
frag_df=read.csv("frag_syn_frag_df_24sep2025.csv", header=TRUE)

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
# write.csv(diameter_df, "frag_syn_diameter_df_24sep2025.csv", row.names = FALSE)


# diameter_df=read.csv("frag_syn_diameter_df_2june2025.csv", header=TRUE)
diameter_df=read.csv("frag_syn_diameter_df_24sep2025.csv", header=TRUE)

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




#### Effect on network Properties ####


# Fracture size
clust_var=ggplot(mean_mean_clust_size, aes(x=variation, y=mean, col=floor(diff_mean/60*100), group=floor(diff_mean/60*100)))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  xlab('Variation')+
  ylab('Mean fracture size')+
  guides(col = guide_legend(title = "Delay"))+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  NULL
clust_var

clust_delay=ggplot(mean_mean_clust_size, aes(x=diff_mean/60*100, y=mean, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean fracture size')+
  guides(col = guide_legend(title = "Variation"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  NULL
clust_delay


# Fracture Proportion
prop_delay=ggplot(mean_mean_frag_prop, aes(x=diff_mean/60*100, y=mean_prop, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean propagule\nproportion')+
  guides(col = guide_legend(title = "Variation"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  NULL
prop_delay

prop_var=ggplot(mean_mean_frag_prop, aes(x=variation, y=mean_prop, col=floor(diff_mean/60*100), group=floor(diff_mean/60*100)))+
  geom_point()+
  geom_line()+
  xlab('Variation')+
  ylab('Mean propagule\nproportion')+
  guides(col = guide_legend(title = "Delay"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  NULL
prop_var


#Network diameter
diam_delay=ggplot(mean_summ_diam, aes(x=diff_mean/60*100, y=mean_d, col=variation, group=variation))+
  geom_point()+
  geom_line()+
  xlab('Delay (% Second Division)')+
  ylab('Mean Network\nDiameter')+
  guides(col = guide_legend(title = "Variation"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_variation+
  theme_classic(base_size = 10)+
  NULL
diam_delay

diam_var=ggplot(mean_summ_diam, aes(x=variation, y=mean_d, col=floor(diff_mean/60*100), group=floor(diff_mean/60*100)))+
  geom_point()+
  geom_line()+
  xlab('Variation')+
  ylab('Mean Network\nDiameter')+
  guides(col = guide_legend(title = "Delay"))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1)+
  color_diff_mean_perc+
  theme_classic(base_size = 10)+
  NULL
diam_var



#### Figure 5 ####

# supp_net_prop_v3=plot_grid(clust_delay,clust_var, prop_delay,prop_var, 
#                            diam_delay,diam_var,
#                            labels=c('A', 'B', 'C', 'D', 'E', 'F'), ncol=2, 
#                            align='hv', label_size=12)
# supp_net_prop_v3
# 
# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_6_effects_delay_variation_25sep2025.png',
#        plot=supp_net_prop_v3, dpi='retina', width=7, height=6, bg='white')


supp_net_prop_v4=plot_grid(clust_delay,clust_var, diam_delay,diam_var, prop_delay,prop_var,
                           labels=c('A', 'B', 'C', 'D', 'E', 'F'), ncol=2, 
                           align='hv', label_size=11)
supp_net_prop_v4

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_6_effects_delay_variation_29mar2026.png',
       plot=supp_net_prop_v4, dpi='retina', width=6.5, height=5.6, bg='white')


# OLD Figure 5 ####

mean_mean_clust_size$var_cat=NA
mean_mean_clust_size$var_cat=ifelse(mean_mean_clust_size$variation==0, '0',
                                    ifelse(mean_mean_clust_size$variation==5, 'l',
                                           ifelse(mean_mean_clust_size$variation==10, 'm', 'h')))
mean_mean_clust_size$var_cat=factor(mean_mean_clust_size$var_cat, levels=c('0', 'l', 'm', 'h'))

mean_mean_clust_size$diff_cat=NA
mean_mean_clust_size$diff_cat=ifelse(mean_mean_clust_size$diff_mean==0, '0',
                                     ifelse(mean_mean_clust_size$diff_mean==20, 'l',
                                            ifelse(mean_mean_clust_size$diff_mean==40, 'm', 'h')))
mean_mean_clust_size$diff_cat=factor(mean_mean_clust_size$diff_cat, levels=c('0', 'l', 'm', 'h'))


clust_size_diff=mean_mean_clust_size[mean_mean_clust_size$var_cat=='l',]
clust_size_diff$category=clust_size_diff$diff_cat
clust_size_diff$Variable='Delay'

clust_size_var=mean_mean_clust_size[mean_mean_clust_size$diff_cat=='l',]
clust_size_var$category=clust_size_var$var_cat
clust_size_var$Variable='Variation'

clust_size_min=min(c(clust_size_var$mean, clust_size_diff$mean))-10
clust_size_max=max(c(clust_size_var$mean, clust_size_diff$mean))+10

clust_p2=ggplot(clust_size_var, aes(x=variation, y=mean))+
  geom_line(col='#DCCA2CFF')+
  # geom_point(col='#DCCA2CFF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(clust_size_min, clust_size_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Fracture\nSize', x='Variation')+
  NULL
clust_p2

clust_p1=ggplot(clust_size_diff, aes(x=diff_mean/60*100, y=mean))+
  geom_line(col='#92BBD9FF')+
  # geom_point(col='#92BBD9FF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(clust_size_min, clust_size_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Fracture\nSize', x='Delay (% Second Division)')+
  NULL
clust_p1


mean_mean_frag_prop$var_cat=NA
mean_mean_frag_prop$var_cat=ifelse(mean_mean_frag_prop$variation==0, '0',
                                   ifelse(mean_mean_frag_prop$variation==5, 'l',
                                          ifelse(mean_mean_frag_prop$variation==10, 'm', 'h')))
mean_mean_frag_prop$var_cat=factor(mean_mean_frag_prop$var_cat, levels=c('0', 'l', 'm', 'h'))

mean_mean_frag_prop$diff_cat=NA
mean_mean_frag_prop$diff_cat=ifelse(mean_mean_frag_prop$diff_mean==0, '0',
                                    ifelse(mean_mean_frag_prop$diff_mean==20, 'l',
                                           ifelse(mean_mean_frag_prop$diff_mean==40, 'm', 'h')))
mean_mean_frag_prop$diff_cat=factor(mean_mean_frag_prop$diff_cat, levels=c('0', 'l', 'm', 'h'))


propag_prop_diff=mean_mean_frag_prop[mean_mean_frag_prop$var_cat=='l',]
propag_prop_diff$category=propag_prop_diff$diff_cat
propag_prop_diff$Variable='Delay'

propag_prop_var=mean_mean_frag_prop[mean_mean_frag_prop$diff_cat=='l',]
propag_prop_var$category=propag_prop_var$var_cat
propag_prop_var$Variable='Variation'

prop_propag_min=min(c(propag_prop_var$mean_prop, propag_prop_diff$mean_prop))-0.005
prop_propag_max=max(c(propag_prop_var$mean_prop, propag_prop_diff$mean_prop))+0.005

propagule_p2=ggplot(propag_prop_var, aes(x=variation, y=mean_prop))+
  geom_line(col='#DCCA2CFF')+
  # geom_point(col='#DCCA2CFF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(prop_propag_min, prop_propag_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Propagule\nProportion', x='Variation')+
  NULL
propagule_p2

propagule_p1=ggplot(propag_prop_diff, aes(x=diff_mean/60*100, y=mean_prop))+
  geom_line(col='#92BBD9FF')+
  # geom_point(col='#92BBD9FF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(prop_propag_min, prop_propag_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Propagule\nProportion', x='Delay (% Second Division)')+
  NULL
propagule_p1


mean_summ_diam$var_cat=NA
mean_summ_diam$var_cat=ifelse(mean_summ_diam$variation==0, '0', 
                              ifelse(mean_summ_diam$variation==5, 'l',
                                     ifelse(mean_summ_diam$variation==10, 'm', 'h')))
mean_summ_diam$var_cat=factor(mean_summ_diam$var_cat, levels=c('0', 'l', 'm', 'h'))

mean_summ_diam$diff_cat=NA
mean_summ_diam$diff_cat=ifelse(mean_summ_diam$diff_mean==0, '0', 
                               ifelse(mean_summ_diam$diff_mean==20, 'l',
                                      ifelse(mean_summ_diam$diff_mean==40, 'm', 'h')))
mean_summ_diam$diff_cat=factor(mean_summ_diam$diff_cat, levels=c('0', 'l', 'm', 'h'))


net_diam_diff=mean_summ_diam[mean_summ_diam$var_cat=='l',]
net_diam_diff$category=net_diam_diff$diff_cat
net_diam_diff$Variable='Delay'

net_diam_var=mean_summ_diam[mean_summ_diam$diff_cat=='l',]
net_diam_var$category=net_diam_var$var_cat
net_diam_var$Variable='Variation'

net_diam_min=min(c(net_diam_var$mean_d, net_diam_diff$mean_d))-0.02
net_diam_max=max(c(net_diam_var$mean_d, net_diam_diff$mean_d))+0.02

diameter_p2=ggplot(net_diam_var, aes(x=variation, y=mean_d))+
  geom_line(col='#DCCA2CFF')+
  # geom_point(col='#DCCA2CFF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(net_diam_min, net_diam_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Diameter', x='Variation')+
  NULL
diameter_p2

diameter_p1=ggplot(net_diam_diff, aes(x=diff_mean/60*100, y=mean_d))+
  geom_line(col='#92BBD9FF')+
  # geom_point(col='#92BBD9FF')+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci), width=0.1, col='black')+
  ylim(c(net_diam_min, net_diam_max))+
  theme_classic(base_size = 10)+
  labs(y='Mean Diameter', x='Delay (% Second Division)')+
  NULL
diameter_p1


legend <- get_legend(
  ggplot() + 
    geom_line(data = data.frame(x = 1, y = 1, group = factor(c("Delay", "Variation"))),
              aes(x = x, y = y, color = group)) +
    guides(color = guide_legend(title = "")) +
    scale_color_manual(values = c('#92BBD9FF', '#DCCA2CFF'), labels = c("Delay", "Variation")) +
    theme(legend.position = "bottom",
          legend.box.just = "center",
          legend.justification = "center",
          legend.text = element_text(size = 8))
)

plot_summary_effects_temp <- plot_grid(clust_p1, clust_p2, diameter_p1, diameter_p2,
                                       propagule_p1, propagule_p2,
                                  labels = c("A", "B", "C", "D", "E", "F"), 
                                  ncol = 2, align = 'v', axis = 'l',
                                  rel_heights = c(1, 1, 1),
                                  rel_widths = c(1, 1),label_size=12)
plot_summary_effects_temp


# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_6_syn_summ_doub_t_effects_30july2025.png',
#        plot=plot_summary_effects_temp, dpi='retina', width=6, height=6, bg='white')



