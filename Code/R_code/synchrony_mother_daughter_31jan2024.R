
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

# cell_sync=read.csv("~/work_dir/timelapses/cell_sync_info_all_tl_2023apr5_v3.csv", header=TRUE)
cell_sync=read.csv("~/work_dir/timelapses/cell_sync_info_2024dec2.csv", header=TRUE)

cell_sync$ordered_timepoint=cell_sync$timepoint
cell_sync$ordered_timepoint=ifelse(cell_sync$strain=='petite', 'petite',
                                   ifelse(cell_sync$strain=='grande', 'grande',
                                          cell_sync$timepoint))
cell_sync$timepoint=factor(cell_sync$timepoint, levels=c('t0', 't200', 't400', 't600', 't1000'))
cell_sync$ordered_timepoint=factor(cell_sync$ordered_timepoint, levels=c('grande', 'petite', 't200', 't400', 't600', 't1000'))
summary(cell_sync)


#Here the petite and the grande are combined in the same table
ggplot(cell_sync, aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.1)+
  facet_wrap(~timepoint)+
  geom_abline()+
  NULL

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

ggplot(cor_res, aes(x=ordered_timepoint, y=correlation, group(ordered_timepoint)))+
  geom_point()+
  xlab('Timepoint')+
  ylab('Correlation')+
  NULL

table(cell_sync$ordered_timepoint)

# Branch synchrony ####

# df_all=read.csv("~/work_dir/timelapses/tl_merged_all_tl_2023apr5.csv", header = TRUE)
df_all=read.csv("~/work_dir/timelapses/tl_merged_2024dec2.csv", header = TRUE)

df_all$minutes=df_all$DELTA_T/60
df_all$hours=df_all$minutes/60
df_all$mother_delta_hours=df_all$mother_delta_t/3600

df_all$division_number=as.factor(df_all$division_number+1)

df_all$transfer=as.numeric(str_replace(df_all$timepoint,'t',''))

df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

#adding variable to make a unique id for each branch
df_all$unique_id_branch=paste(df_all$complete_id,df_all$TRACK_ID,sep='_')

summary(df_all)

doubled=df_all[df_all$N_PREDECESSORS==1 & df_all$N_SUCCESSORS==2,]

#filtering out slow dividing tl ####
for (i in unique(doubled$strain)){
  print(ggplot(doubled[doubled$strain==i,], aes(x=complete_id, y=hours, fill=timepoint))+
          geom_violin()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          ggtitle(paste("Individual time lapses for",i))+
          NULL)
}

#Filtering out the time lapses where the cells take longer to divide
doubled_filter <- doubled[!(doubled$complete_id %in% c("PA1_t1000_3_2022sep27","PA4_t200_2_2022nov23",
                                                       "PA3_t1000_7_2022dec1")), ]

# normalizing the frames ####
normal_frames=data.frame(matrix(ncol = dim(doubled_filter)[2], nrow = 0))
colnames(normal_frames)=colnames(doubled_filter)
for (i in unique(doubled_filter$unique_id_branch)){
  temp_df=doubled_filter[doubled_filter$unique_id_branch==i & doubled_filter$FRAME!=180,]
  
  if(dim(temp_df)[1]>2){
    #normalizing the FRAME value by making it start at 0
    first_frame=min(temp_df$FRAME)
    temp_df$FRAME=temp_df$FRAME-first_frame
    normal_frames=rbind(normal_frames,temp_df)
  }
}

ids_order=c("grande_t0","petite_t0","PA1_t200","PA1_t400","PA1_t600","PA1_t1000","PA2_t200","PA2_t400","PA2_t600","PA2_t1000",
            "PA3_t200","PA3_t400" ,"PA3_t600" ,"PA3_t1000","PA4_t200","PA4_t400","PA4_t600","PA4_t1000","PA5_t200","PA5_t400",
            "PA5_t600","PA5_t1000")  

normal_frames$strain_timepoint=factor(normal_frames$strain_timepoint, levels=ids_order)

frames_by_strain=ggplot(normal_frames, aes(x = FRAME, y=strain_timepoint, fill = strain)) +
  geom_density_ridges(stat="binline",binwidth=1)+
  theme_ridges()+ 
  theme(legend.position = "none")+
  NULL
frames_by_strain

frames_by_strain_200=ggplot(normal_frames[normal_frames$transfer<=200,], 
                            aes(x = FRAME, y=strain_timepoint, fill = strain)) +
  geom_density_ridges(stat="binline",binwidth=1)+
  theme_ridges()+
  theme(legend.position = "none")+
  NULL
frames_by_strain_200

for(i in c(200,400,600,1000)){
  frames_by_strain_timepoint=ggplot(normal_frames[normal_frames$transfer==i,], 
                                    aes(x = FRAME, y=strain_timepoint, fill = strain)) +
    geom_density_ridges(stat="binline",binwidth=1)+
    theme_ridges()+
    theme(legend.position = "none")+
    ggtitle(paste("Doubling distribution by frame t",as.character(i),sep=""))+
    NULL
  print(frames_by_strain_timepoint)
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

ggplot(df_all[df_all$id_file=='2022nov22_gob440_2' & df_all$N_SUCCESSORS==2,], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram()+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL

# Looking at which set of grande strain data points is better

ggplot(normal_frames[normal_frames$strain=='grande',],
       aes(x=5*FRAME, fill=strain))+
  geom_histogram(binwidth = 10, position = 'identity')+
  facet_wrap(~id_file)+
  xlim(c(-15, 350))+
  NULL

old_data_points=c('2022nov22_gob8_6', '2022nov22_gob8_7', '2022nov23_gob8_7')

old_grande_strain=normal_frames[normal_frames$id_file %in% old_data_points,]
new_grande_strain=normal_frames[normal_frames$strain=='grande' & ! normal_frames$id_file %in% old_data_points,]

ggplot(old_grande_strain, aes(x = FRAME*5, fill = strain))+
  geom_histogram(binwidth = 10, position = 'identity')+
  xlim(c(-15, 350))+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL

ggplot(new_grande_strain, aes(x = FRAME*5, fill = strain))+
  geom_histogram(binwidth = 10, position = 'identity')+
  xlim(c(-15, 350))+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  NULL

# Grande data comparison ####
# Checking number of doublings for each batch of data

old_doubled_grande=doubled_filter[doubled_filter$strain=='grande' & doubled_filter$id_file %in% old_data_points,]
new_doubled_grande=doubled_filter[doubled_filter$strain=='grande' & ! doubled_filter$id_file %in% old_data_points,]

table(old_doubled_grande$division_number)
#   1   2   3   4   5   6   7   8 
# 144  36   3   0   0   0   0   0

dim(old_doubled_grande) # 183

table(new_doubled_grande$division_number)
#   1   2   3   4   5   6   7   8 
# 147  53   3   0   0   0   0   0

dim(new_doubled_grande) # 203

# Paper figure 1 plot ####

img <- readPNG("~/work_dir/observed_synchrony/evolution_results/images_timelapse/t200_montage/15min_montage/gob440_montage_outline_scale_time_15min_4px_line.png")
img_plot <- rasterGrob(img, interpolate = TRUE)


p2=ggplot(df_all[df_all$id_file=='2022nov22_gob440_2' & df_all$N_SUCCESSORS==2,], 
       aes(x = FRAME*5, fill = strain))+
  geom_histogram(fill='lightblue')+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Frequency')+
  theme_classic(base_size = 11)+
  NULL
p2

p2_norm=ggplot(normal_frames[normal_frames$id_file=='2022nov22_gob440_2',], 
               aes(x = FRAME*5, fill = strain))+
  geom_histogram(fill='lightblue')+
  theme(legend.position = "none")+
  xlab('Minutes')+
  ylab('Number of\nDivisions')+
  theme_classic(base_size = 11)+
  NULL
p2_norm

p3=ggplot(cell_sync[cell_sync$strain!='grande',], 
          aes(x=mother_minutes, y=daughter_minutes, col=timepoint))+
  geom_point(alpha=0.5)+
  facet_wrap(~timepoint, ncol=5)+
  geom_abline()+
  guides(col='none')+
  xlab('Mother Doubling Time (minutes)')+
  ylab('Daughter Doubling\nTime (minutes)')+
  theme_classic(base_size = 11)+
  NULL
p3


figure_1=plot_grid(img_plot,p2_norm,p3, labels=c('A)', 'B)', 'C)', 'D)'), ncol=1, align='hv', label_size=11, rel_heights=c(1,1.5,1.5))
figure_1

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/fig_1_synchronization_strains_26feb2024_4px.png',
#        plot=figure_1, dpi='retina', width=10, height=12, bg='white')



#### Updated Figure 1 ####


ggplot(normal_frames[normal_frames$id_file=='2022nov22_gob440_2' | normal_frames$id_file=='2022sep30_gob21_6',], 
       aes(x = FRAME, y=strain_timepoint, fill = strain)) +
  geom_density_ridges(stat="binline", binwidth=1)+
  theme_ridges()+ 
  theme(legend.position = "none")+
  NULL


p2_ancestor_data=normal_frames[normal_frames$id_file=='2022sep30_gob21_6' | normal_frames$id_file=='2022nov22_gob440_2',]
table(p2_ancestor_data$strain)
p2_ancestor_data$strain=factor(p2_ancestor_data$strain, levels=c("petite", "PA3"))


p2_ancestor_paper_v1=ggplot(p2_ancestor_data, 
                         aes(x = 5*FRAME, fill = strain)) +
  # geom_histogram(position='dodge', binwidth = 20)+
  geom_histogram(binwidth = 10, alpha=0.75, position = 'identity')+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(fill = "Strain")+
  scale_fill_manual(labels = c("t0 (Ancestor)", "t200"),
                    values = c("#ff767e", "#add8e6"))+
  xlim(c(-15, 350))+
  NULL
p2_ancestor_paper_v1

p2_ancestor_paper_v2=ggplot(p2_ancestor_data, 
                   aes(x = 5*FRAME, fill = strain)) +
  # geom_histogram(position='dodge', binwidth = 20)+
  xlim(c(-15, 350))+
  scale_y_continuous(breaks = c(0, 3, 6, 9), limits=c(0,11))+
  geom_histogram(binwidth = 10, position = 'identity')+
  facet_wrap(~strain, ncol=1, scales='free')+
  xlab('Minutes')+
  ylab('Number of Divisions')+
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))+
  labs(fill = "Strain")+
  scale_fill_manual(labels = c("t0 (Ancestor)", "t200"),
                    values = c("#ff767e", "#add8e6"))+
  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  NULL
p2_ancestor_paper_v2



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

petite_dt=doubled_filter[doubled_filter$strain=='petite' & as.numeric(doubled_filter$division_number)<=2,]
petite_dt$name='t0'

t200_dt_conf=doubled_filter[doubled_filter$strain=='PA2' & doubled_filter$timepoint=='t200' & as.numeric(doubled_filter$division_number)<=2,]
t200_dt_conf$name='t200'

t400_dt_conf=doubled_filter[doubled_filter$strain=='PA2' & doubled_filter$timepoint=='t400' & as.numeric(doubled_filter$division_number)<=2,]
t400_dt_conf$name='t400'

t600_dt_conf=doubled_filter[doubled_filter$strain=='PA2' & doubled_filter$timepoint=='t600' & as.numeric(doubled_filter$division_number)<=2,]
t600_dt_conf$name='t600'

t1000_dt_conf=doubled_filter[doubled_filter$strain=='PA2' & doubled_filter$timepoint=='t1000' & as.numeric(doubled_filter$division_number)<=2,]
t1000_dt_conf$name='t1000'

ancestor_PA2_dt=rbind(petite_dt, t200_dt_conf, t400_dt_conf, t600_dt_conf, t1000_dt_conf)
summary(ancestor_PA2_dt)


p4_doubling_times=ggplot(doubled_filter[as.numeric(doubled_filter$division_number)<=2 & doubled_filter$strain!='grande',], 
                         aes(x=division_number, y=hours, fill=timepoint))+
  facet_wrap(~timepoint, ncol=5)+
  geom_boxplot()+
  # geom_violin(adjust=2)+
  # stat_summary(fun='mean', geom='crossbar')+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none')+
  petite_t200_colors+
  scale_fill_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_doubling_times


p4_doubling_times_violin=ggplot(doubled_filter[as.numeric(doubled_filter$division_number)<=2 & doubled_filter$strain!='grande',], 
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


p4_PA2_dt=ggplot(ancestor_PA2_dt, aes(x=division_number, y=hours, fill=name))+
  facet_wrap(~name, ncol=5)+
  geom_boxplot()+
  # stat_summary(fun='mean', geom='crossbar')+
  xlab('Division Number')+
  ylab('Hours')+
  guides(fill='none')+
  petite_t200_colors+
  scale_fill_manual(values=c("#ff767e","#add8e6","#00BF7D","#00B0F6","#E76BF3"))+
  # theme_classic(base_size = 11)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )+
  NULL
p4_PA2_dt


figure_1_paper=plot_grid(img_plot,p2_ancestor_paper_v1,p3_paper, labels=c('A)', 'B)', 'C)'), ncol=1, align='hv', label_size=16, rel_heights=c(1,1.5,1.5))
figure_1_paper


# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/fig_1_synchronization_strains_v2_12jun2024.png',
#        plot=figure_1_paper, dpi='retina', width=10, height=12, bg='white')


figure_1_paper_v2=plot_grid(img_plot,p2_ancestor_paper_v2,p3_paper, p4_doubling_times_violin, 
                            labels=c('A)', 'B)', 'C)', 'D)'), ncol=1, label_size=16, 
                            rel_heights=c(1.1,2,1.5,1.5))
figure_1_paper_v2

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_1_synchronization_strains_27sep2024.png',
#        plot=figure_1_paper_v2, dpi='retina', width=10, height=12, bg='white')


#### Evolution synchrony plot ####
# Plot evolution of synchronous cell divisions in the Anaerobic line


figure_evol_sync=plot_grid(p3_paper, p4_doubling_times_violin, 
                            labels=c('A)', 'B)'), ncol=1, label_size=16)
figure_evol_sync

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_evolution_synchrony_PA_line_16oct2024.png',
#        plot=figure_evol_sync, dpi='retina', width=10, height=6, bg='white')

# Ancestors plot ####


img_petite <- readPNG("~/work_dir/observed_synchrony/evolution_results/images_timelapse/petite_montage/petite_1hr_montage_time_v1.png")
img_plot_petite <- rasterGrob(img_petite, interpolate = TRUE)

img_grande <- readPNG("~/work_dir/observed_synchrony/evolution_results/images_timelapse/grande_montage/grande_montage_time_4f.png")
img_plot_grande <- rasterGrob(img_grande, interpolate = TRUE)



ancestor_tl_dt=normal_frames[normal_frames$id_file=='2022nov22_gob8_6' | normal_frames$id_file=='2022sep27_gob21_4',]
table(ancestor_tl_dt$strain)
ancestor_tl_dt$strain=factor(ancestor_tl_dt$strain, levels=c("petite", "grande"))


p2_ancestors=ggplot(ancestor_tl_dt, 
                            aes(x = 5*FRAME, fill = strain)) +
  # geom_histogram(position='dodge', binwidth = 20)+
  xlim(c(-15, 350))+
  scale_y_continuous(breaks = c(0, 3, 6, 9), limits=c(0,11))+
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
p2_ancestors

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


petite_dt=doubled_filter[doubled_filter$strain=='petite' & as.numeric(doubled_filter$division_number)<=2,]
petite_dt$name='Petite'

grande_dt=doubled_filter[doubled_filter$strain=='grande' & as.numeric(doubled_filter$division_number)<=2,]
grande_dt$name='Grande'

p4_ancestors_dt=rbind(petite_dt, grande_dt)
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


figure_1_paper_ancestors=plot_grid(img_plot_petite, img_plot_grande, 
                                   p2_ancestors, p4_doubling_times_violin_ancestors, 
                            labels=c('A)', 'B)', 'C)', 'D)'), ncol=1, label_size=16, 
                            rel_heights=c(1.6,1,1.4,1.4))
figure_1_paper_ancestors

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_1_synchronization_ancestors_11oct2024.png',
#        plot=figure_1_paper_ancestors, dpi='retina', width=10, height=14, bg='white')


# Create the first column with two plots
column1 <- plot_grid(img_plot_petite, img_plot_grande, 
                     labels = c('A)', 'B)'), 
                     ncol = 1, 
                     label_size = 16, 
                     rel_heights = c(1.6, 1))

column2=plot_grid(p3_ancestors, p4_doubling_times_violin_ancestors, 
                  labels=c('D)', 'E)'),
                  ncol=1,
                  label_size = 16,
                  rel_heights = c(1,1))

# Create the second row with two plots
row2 <- plot_grid(p2_ancestors, p4_doubling_times_violin_ancestors,
                  labels = c('C)', 'D)'),
                  ncol = 2,
                  label_size = 16,
                  rel_widths = c(1.4, 1.4))


row_v2=plot_grid(p2_ancestors, column2,
                    labels = c('C)', ''),
                    ncol = 2,
                    label_size = 16,
                    rel_widths = c(1.4, 1.4))

# Figure 1 without doubling time distributions
figure_1_paper_ancestors <- plot_grid(column1, row2, 
                                      ncol = 1, 
                                      rel_heights = c(2.6, 1.6))
figure_1_paper_ancestors

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_1_synchronization_ancestors_11oct2024.png',
#        plot=figure_1_paper_ancestors, dpi='retina', width=10, height=14, bg='white')


# Figure 1 with doubling time distributions
figure_1_paper_ancestors_dt <- plot_grid(column1, row_v2, 
                                      ncol = 1, 
                                      rel_heights = c(2.6, 1.6))
figure_1_paper_ancestors_dt

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_1_synchronization_ancestors_dt_16oct2024.png',
#        plot=figure_1_paper_ancestors_dt, dpi='retina', width=10, height=14, bg='white')


# Ancestors plot v2 ####


img_petite_v2 <- readPNG("~/work_dir/observed_synchrony/evolution_results/images_timelapse/petite_montage_v2/gob21_2hr_5f_scale_masks.png")
img_plot_petite_v2 <- rasterGrob(img_petite_v2, interpolate = TRUE)

img_grande_v2 <- readPNG("~/work_dir/observed_synchrony/evolution_results/images_timelapse/grande_montage_v2/grande_1hr_5f_scale_0.5x_masks.png")
img_plot_grande_v2 <- rasterGrob(img_grande_v2, interpolate = TRUE)


old_data_points=c('2022nov22_gob8_6', '2022nov22_gob8_7', '2022nov23_gob8_7')
ancestor_tl_dt=normal_frames[normal_frames$strain=='petite' | (normal_frames$strain=='grande' & ! normal_frames$id_file %in% old_data_points),]
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


petite_dt=doubled_filter[doubled_filter$strain=='petite' & as.numeric(doubled_filter$division_number)<=2,]
petite_dt$name='Petite'

grande_dt=doubled_filter[doubled_filter$strain=='grande' & as.numeric(doubled_filter$division_number)<=2
                         & !doubled_filter$id_file %in% old_data_points,]
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

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/restructured_paper/fig_1_synchronization_ancestors_dt_4nov2024.png',
#        plot=figure_1_paper_ancestors_dt_v2, dpi='retina', width=10, height=12, bg='white')

# ggsave(filename='~/work_dir/observed_synchrony/evolution_results/images_paper/updated_grande/fig_1_synchronization_ancestors_dt_4nov2024.png',
#        plot=figure_1_paper_ancestors_dt_v2, dpi='retina', width=10, height=12, bg='white')
