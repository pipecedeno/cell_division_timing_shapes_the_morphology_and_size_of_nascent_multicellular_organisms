
library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
# library(glue)
library(purrr)
library(stringr)
library(ggridges)
library(ghibli)
# library(ggpubr)
# library(tidyverse)
library(ggnewscale)
library(cowplot)
library(png)
library(grid)

theme_set(theme_classic(base_size = 16))
mycolors <- rev(ghibli_palettes$LaputaMedium)
petite_t200_colors <- list(scale_color_manual(values = mycolors), scale_fill_ghibli_d("LaputaMedium", direction = -1))

cell_sync=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/cell_sync_mother_daughter_2024dec2.csv", header=TRUE)

cell_sync$ordered_timepoint=cell_sync$timepoint
cell_sync$ordered_timepoint=ifelse(cell_sync$strain=='petite', 'petite',
                                   ifelse(cell_sync$strain=='grande', 'grande',
                                          cell_sync$timepoint))
cell_sync$timepoint=factor(cell_sync$timepoint, levels=c('t0', 't200', 't400', 't600', 't1000'))
cell_sync$ordered_timepoint=factor(cell_sync$ordered_timepoint, levels=c('grande', 'petite', 't200', 't400', 't600', 't1000'))
summary(cell_sync)



ggplot(cell_sync, aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.1)+
  facet_wrap(~ordered_timepoint)+
  geom_abline()+
  xlab('Mother doubling time (min)')+
  ylab('Daughter doubling time (min)')+
  NULL

cor_res=cell_sync %>%
  group_by(ordered_timepoint) %>%
  summarise(correlation=cor(mother_minutes, daughter_minutes))
cor_res
# ordered_timepoint correlation
# 1 grande                  0.914
# 2 petite                  0.508
# 3 t200                    0.844
# 4 t400                    0.789
# 5 t600                    0.891
# 6 t1000                   0.749

table(cell_sync$ordered_timepoint)

# Branch synchrony ####

df_all=read.csv("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/timelapse_doubling_time_inf_2024dec2_clean.csv", header = TRUE)
df_all$division_number=as.factor(df_all$division_number)
df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

summary(df_all)


doubled=df_all[df_all$N_PREDECESSORS==1 & df_all$N_SUCCESSORS==2,]


# normalizing the frames ####
# The idea is to make the doubling time of the cell 
normal_frames=data.frame(matrix(ncol = dim(doubled)[2], nrow = 0))
colnames(normal_frames)=colnames(doubled)
for (i in unique(doubled$unique_id_branch)){
  temp_df=doubled[doubled$unique_id_branch==i & doubled$FRAME!=180,]
  
  if(dim(temp_df)[1]>2){
    #normalizing the FRAME value by making it start at 0
    first_frame=min(temp_df$FRAME)
    temp_df$FRAME=temp_df$FRAME-first_frame
    normal_frames=rbind(normal_frames,temp_df)
  }
}


# Selected time lapse plots ####

summary(normal_frames)

table(normal_frames[normal_frames$strain_timepoint=='petite_t0',]$id_file)

ggplot(normal_frames[normal_frames$id_file=='2022nov22_gob440_2',], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram()+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL


ggplot(normal_frames[normal_frames$id_file=='2022sep30_gob21_6',], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram()+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL




# Figure 1 ####


img_petite_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/petite_2hr_5f_scale_masks.png")
img_plot_petite_v2 <- rasterGrob(img_petite_v2, interpolate = TRUE)

img_grande_v2 <- readPNG("~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/Images/grande_1hr_5f_scale_0.5x_masks.png")
img_plot_grande_v2 <- rasterGrob(img_grande_v2, interpolate = TRUE)


ancestor_tl_dt=normal_frames[normal_frames$strain=='petite' | normal_frames$strain=='grande',]
table(ancestor_tl_dt$strain)
ancestor_tl_dt$strain=factor(ancestor_tl_dt$strain, levels=c("petite", "grande"))


p2_ancestors_v2=ggplot(ancestor_tl_dt, 
                    aes(x = 5*FRAME, fill = strain)) +
  # geom_histogram(position='dodge', binwidth = 20)+
  xlim(c(-15, 350))+
  # scale_y_continuous(breaks = c(0, 3, 6, 9), limits=c(0,11))+
  geom_histogram(binwidth = 10, position = 'identity')+
  facet_wrap(~strain, ncol=1, scales='free')+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(fill = "Strain")+
  scale_fill_manual(labels = c("Petite", "Grande"),
                    values = c("#add8e6", "#ff767e"))+
  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  NULL
p2_ancestors_v2


cell_sync_ancestors=cell_sync[cell_sync$strain=='grande' | cell_sync$strain=='petite',]
cell_sync_ancestors$strain=ifelse(cell_sync_ancestors$strain=='petite', 'Petite', 'Grande')

p3_ancestors=ggplot(cell_sync_ancestors, 
                    aes(x=mother_minutes, y=daughter_minutes, col=strain))+
  geom_point(alpha=0.5)+
  geom_abline(linetype='dashed')+
  facet_wrap(~strain, ncol=5)+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter Doubling\nTime (minutes)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  # scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p3_ancestors

petite_dt=doubled[doubled$strain=='petite' & as.numeric(doubled$division_number)<=2,]
petite_dt$name='Petite'

grande_dt=doubled[doubled$strain=='grande' & as.numeric(doubled$division_number)<=2,]
grande_dt$name='Grande'

p4_ancestors_dt=rbind(petite_dt, grande_dt)
table(p4_ancestors_dt$strain)
summary(p4_ancestors_dt)


p4_doubling_times_violin_ancestors=ggplot(p4_ancestors_dt, 
                                          aes(x=division_number, y=hours, fill=timepoint, col=name))+
  facet_wrap(~name, ncol=5)+
  # geom_boxplot()+
  geom_jitter(alpha=0.5)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  xlab('Division Number')+
  ylab('\nHours')+
  guides(fill='none', col='none')+
  petite_t200_colors+
  scale_fill_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times_violin_ancestors


# Create the first column with two plots
column1 <- plot_grid(img_plot_petite_v2, img_plot_grande_v2, 
                     labels = c('A)', 'B)'), 
                     ncol = 1, 
                     label_size = 16, 
                     rel_heights = c(0.8, 1))

column2=plot_grid(p3_ancestors, p4_doubling_times_violin_ancestors, 
                  labels=c('D)', 'E)'),
                  ncol=1,
                  label_size = 16,
                  rel_heights = c(1,1))


row_v2=plot_grid(p2_ancestors_v2, column2,
                 labels = c('C)', ''),
                 ncol = 2,
                 label_size = 16,
                 rel_widths = c(1.4, 1.4))

figure_1_paper_ancestors_dt_v2 <- plot_grid(column1, row_v2, 
                                         ncol = 1, 
                                         rel_heights = c(0.8, 1))
figure_1_paper_ancestors_dt_v2


ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_1_synchronization_ancestors_dt_4nov2024.png',
       plot=figure_1_paper_ancestors_dt_v2, dpi='retina', width=10, height=12, bg='white')


# Figure 7 ####

p3_paper=ggplot(cell_sync[cell_sync$strain!='grande',], 
                aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.5)+
  geom_abline(linetype='dashed')+
  facet_wrap(~timepoint, ncol=5)+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter Doubling\nTime (minutes)')+
  scale_x_continuous(breaks = seq(100, 300, by = 100))+
  scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p3_paper

p4_doubling_times_violin=ggplot(doubled[as.numeric(doubled$division_number)<=2 & doubled$strain!='grande',], 
                                aes(x=division_number, y=hours, fill=timepoint, col=timepoint))+
  facet_wrap(~timepoint, ncol=5)+
  # geom_boxplot()+
  geom_jitter(alpha=0.2)+
  # geom_violin(adjust=2, color='black')+
  stat_summary(fun='mean', geom='crossbar', col='black')+
  xlab('Division Number')+
  ylab('\nHours')+
  guides(fill='none', col='none')+
  petite_t200_colors+
  scale_fill_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  scale_color_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times_violin

figure_evol_sync=plot_grid(p3_paper, p4_doubling_times_violin, 
                           labels=c('A)', 'B)'), ncol=1, label_size=16)
figure_evol_sync

ggsave(filename='~/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Paper_figures/fig_7_evolution_synchrony_PA_line_16oct2024.png',
       plot=figure_evol_sync, dpi='retina', width=10, height=6, bg='white')

