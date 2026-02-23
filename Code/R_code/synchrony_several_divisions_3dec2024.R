#Date 2 December 2024
#The following code is used to generate the plot of the doubling time distribution by the
# number of divisions which is Supplementary Figure 1 in the paper

library(ggplot2)
library(ggpubr)
# library(ggridges)
library(stringr)
library(plyr)  # load before dplyr
library(dplyr)
library(stringi)
library(purrr)
library(ghibli)
library(tidyr)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

#loading data set
doubled=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/first_second_div_doubling_data_2024dec2.csv", header = TRUE)
doubled$timepoint=factor(doubled$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))
doubled$division_number=factor(doubled$division_number)


# Doubling time ####

ggplot(doubled[doubled$timepoint!='t0',], aes(x=timepoint,y=hours,fill=timepoint))+
  geom_violin() +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  facet_wrap(~strain)+
  ylim(c(1,4))+
  NULL


dt_by_div_num=doubled %>%
  group_by(strain_timepoint, division_number) %>%
  summarise(mean_dt = round(mean(hours), 2), .groups = "drop") %>%
  pivot_wider(
    names_from = division_number,
    values_from = mean_dt,
    names_prefix = "Division_"
  )

# write.csv(dt_by_div_num, "~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/summary_doubling_time_per_division_number.csv")





# Doubling time vs Times a cell has divided ####

table(doubled$division_number, doubled$strain_timepoint)

for(temp_strain in sort(unique(doubled$strain))){
  div_num_vs_time=ggplot(doubled[doubled$strain==temp_strain,]
                         , aes(x=division_number, y=hours, fill=division_number)) +
    geom_violin()+
    theme_bw()+
    #geom_jitter(color="black", size=0.4, alpha=0.9)+
    facet_wrap(~timepoint)+
    ggtitle(temp_strain)+
    xlab('Division Number')+ylab('Doubling Time (Hours)')+
    guides(fill = guide_legend(title = "Strain"))+
    NULL
  print(div_num_vs_time)
}


# Selecting t600 because it is the time point with more data for all replicate populations

tapply(doubled[doubled$timepoint=='t600',]$division_number, doubled[doubled$timepoint=='t600',]$strain, table)

# Keeping distributions that have more than 20 data points
t600_doubled=doubled[(doubled$strain_timepoint=='PA1_t600' & as.numeric(doubled$division_number)<=4) |
                       (doubled$strain_timepoint=='PA2_t600' & as.numeric(doubled$division_number)<=4) |
                       (doubled$strain_timepoint=='PA3_t600' & as.numeric(doubled$division_number)<=3) | 
                       (doubled$strain_timepoint=='PA4_t600' & as.numeric(doubled$division_number)<=3) |
                       (doubled$strain_timepoint=='PA5_t600' & as.numeric(doubled$division_number)<=4),]

ggplot(t600_doubled, 
       aes(x=division_number, y=hours, fill=division_number)) +
  geom_violin(adjust=2)+
  facet_wrap(~strain_timepoint)+
  xlab('Division Number')+ylab('Doubling Time (Hours)')+
  guides(fill='none')+
  theme_classic()+
  stat_summary(fun='mean', geom='crossbar')+
  NULL

table(t600_doubled$strain_timepoint, t600_doubled$division_number)
#             1    2    3    4    5    6    7
# PA1_t600 1115  519  131   28    0    0    0
# PA2_t600  634  285   95   27    0    0    0
# PA3_t600  504  221   61    0    0    0    0
# PA4_t600  379  199   64    0    0    0    0
# PA5_t600 1221  599  219   55    0    0    0



# Supplementary Figure ####

div_num_t600=ggplot(t600_doubled, 
       aes(x=division_number, y=hours, fill=division_number)) +
  geom_violin(adjust=2)+
  facet_wrap(~strain_timepoint)+
  xlab('Division Number')+ylab('Doubling Time (Hours)')+
  guides(fill='none')+
  theme_classic()+
  stat_summary(fun='mean', geom='crossbar')+
  # petite_t200_colors+
  scale_fill_manual(values=c("#BBA78CFF", "#333544FF", "#B50A2AFF", "#44A57CFF"))+
  NULL
div_num_t600

# ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/supp_fig1_similar_mean_several_div_3oct2025.png',
#        plot=div_num_t600, dpi='retina', width=6, height=4, bg='white')

