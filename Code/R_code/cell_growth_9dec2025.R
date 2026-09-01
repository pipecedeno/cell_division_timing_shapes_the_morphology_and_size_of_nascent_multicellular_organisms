


library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(effsize)

# Cell properties over time ####

ancestor_spots=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_1_cell_growth/cell_growth_data_ancestors_9dec2025.csv", 
                        header=TRUE)


summary(ancestor_spots)


ancestor_spots$Strain=factor(ancestor_spots$Strain,
                                 levels=c("Grande","Petite"))

colnames(ancestor_spots)


p_dt_diam=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
       aes(x=dt_minutes, y=ELLIPSE_MINOR*2, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Doubling Time (min)", y=expression("Cell Diameter"~(µ*m)))+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  NULL
p_dt_diam

p_cell_cycle_diam=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
       aes(x=percentage_cell_cycle, y=ELLIPSE_MINOR*2, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  geom_smooth(aes(group=1), method="loess", color="black", linewidth=0.75, se=FALSE)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Proportion of Cell Cycle", y=expression("Cell Diameter"~(µ*m)))+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  scale_x_continuous(labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x)))+
  NULL
p_cell_cycle_diam


p_dt_asp_r=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
            aes(x=dt_minutes, y=ELLIPSE_ASPECTRATIO, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Doubling Time (min)", y="Cell Aspect Ratio")+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  NULL
p_dt_asp_r

p_cell_cycle_asp_r=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
                    aes(x=percentage_cell_cycle, y=ELLIPSE_ASPECTRATIO, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  geom_smooth(aes(group=1), method="loess", color="black", linewidth=0.75, se=FALSE)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Proportion of Cell Cycle", y="Cell Aspect Ratio")+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  scale_x_continuous(labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x)))+
  NULL
p_cell_cycle_asp_r

p_dt_area=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
               aes(x=dt_minutes, y=AREA, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Doubling Time (min)", y=expression("Cell Area"~(µ*m^2)))+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  NULL
p_dt_area

p_cell_cycle_area=ggplot(ancestor_spots[ancestor_spots$division_number<=2,], 
                       aes(x=percentage_cell_cycle, y=AREA, col=Strain))+
  geom_point(alpha=0.3, size=0.5)+
  geom_smooth(aes(group=1), method="loess", color="black", linewidth=0.75, se=FALSE)+
  theme_classic(base_size = 8)+
  facet_grid(Strain~division_number)+
  labs(x="Proportion of Cell Cycle", y=expression("Cell Area"~(µ*m^2)))+
  scale_color_manual(labels = c("Grande", "Petite"),
                     values = c("#AE93BEFF", "#F0D77BFF")) +
  guides(col='none')+
  scale_x_continuous(labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sprintf("%.2f", x)))+
  NULL
p_cell_cycle_area

# plot_grid(p_dt_diam, p_dt_asp_r, p_dt_area, p_cell_cycle_diam, p_cell_cycle_asp_r, p_cell_cycle_area,
#           ncol=3)

p_cell_props=plot_grid(p_dt_diam, p_cell_cycle_diam, p_dt_asp_r, p_cell_cycle_asp_r, p_dt_area, p_cell_cycle_area,
                       ncol=2, labels=LETTERS[1:6], label_size=10)
p_cell_props

# ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_cell_growth_properties_9dec2025.png',
#        plot=p_cell_props, dpi='retina', width=6, height=6, bg='white')


#### Supplementary Figure 1 ####

p_cell_props_v2=plot_grid(p_dt_asp_r, p_dt_diam,
                       ncol=2, labels=LETTERS[1:6], label_size=10)
p_cell_props_v2

# ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_1_cell_growth_properties_28aug2026.png',
#        plot=p_cell_props_v2, dpi='retina', width=6.5, height=3, bg='white')

ancestor_spots[ancestor_spots$division_number==2 &
                 ancestor_spots$ELLIPSE_ASPECTRATIO!=1,] %>%
  group_by(strain, division_id) %>%
  summarise(median_area=median(AREA)) %>%
  group_by(strain) %>%
  summarise(median_median_area=median(median_area))
#   strain median_median_area
# 1 grande               23.3
# 2 petite               16.8
# From static measurements:
# petite: 14.95149
# grande: 16.9404





#### Cell properties statistics ####

colnames(ancestor_spots)

median_cell_properties=ancestor_spots[ancestor_spots$division_number==2 &
                                      ancestor_spots$ELLIPSE_ASPECTRATIO!=1,] %>%
  group_by(strain, division_id) %>%
  summarise(median_ar=median(ELLIPSE_ASPECTRATIO),
            median_diameter=median(ELLIPSE_MAJOR*2))

ggplot(median_cell_properties, aes(x=strain, y=median_ar, fill=strain))+
  geom_violin()+
  theme_classic()+
  labs(x="Strain", y="Median Cell Aspect Ratio")+
  scale_fill_manual(values = c("grande" = "#AE93BEFF", "petite" = "#F0D77BFF"))+
  stat_summary(fun='median', geom='crossbar')+
  NULL

ggplot(median_cell_properties, aes(x=strain, y=median_diameter, fill=strain))+
  geom_violin()+
  theme_classic()+
  labs(x="Strain", y="Median Cell Diameter")+
  scale_fill_manual(values = c("grande" = "#AE93BEFF", "petite" = "#F0D77BFF"))+
  stat_summary(fun='median', geom='crossbar')+
  NULL

table(median_cell_properties$strain)

median_cell_properties %>%
  group_by(strain) %>%
  summarise(median_asp_ratio=median(median_ar),
            median_cell_diam=median(median_diameter))
#   strain median_asp_ratio median_cell_diam
# 1 grande             1.16             5.97
# 2 petite             1.18             4.97

t.test(median_cell_properties$median_diameter~median_cell_properties$strain)
# p-value = 2.77e-12

cohen.d(median_cell_properties$median_diameter~median_cell_properties$strain)
# d estimate: 1.83729 (large)

wilcox.test(median_cell_properties$median_diameter~median_cell_properties$strain)
# p-value = 1.256e-11

cliff.delta(median_cell_properties$median_diameter~median_cell_properties$strain)
# delta estimate: 0.8049689 (large)

t.test(median_cell_properties$median_ar~median_cell_properties$strain)
# p-value = 0.7321

cohen.d(median_cell_properties$median_ar~median_cell_properties$strain)
# d estimate: 0.07928122 (negligible)

wilcox.test(median_cell_properties$median_ar~median_cell_properties$strain)
# p-value = 0.8907

cliff.delta(median_cell_properties$median_ar~median_cell_properties$strain)
# delta estimate: 0.01863354 (negligible)


# Budding angle ####


new_cell_position=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/supp_fig_1_cell_growth/cell_budding_angles_all_divisions_30july2025.csv", 
                           header = TRUE)



ggplot(new_cell_position[(new_cell_position$strain=='grande' | 
                            new_cell_position$strain=='petite') &
                           !is.na(new_cell_position$angle_to_bud_scar),],
       aes(x=angle_to_bud_scar, color=strain))+
  stat_ecdf(geom = "step", linewidth = 1)+
  scale_color_manual(values = c("grande" = "#AE93BEFF", 
                                "petite" = "#F0D77BFF"))+
  labs(y = "Cumulative Proportion")+
  theme_classic()+
  NULL

ggplot(new_cell_position[(new_cell_position$strain=='grande' | 
                            new_cell_position$strain=='petite') &
                           !is.na(new_cell_position$angle_to_bud_scar),],
       aes(x=angle_to_bud_scar, color=strain))+
  stat_ecdf(geom = "step", linewidth = 1)+
  scale_color_manual(values = c("grande" = "#AE93BEFF", 
                                "petite" = "#F0D77BFF"))+
  labs(y = "Cumulative Proportion")+
  theme_classic()+
  facet_wrap(~division_number)+
  NULL

# Extract the angle_to_bud_scar values for each strain
grande_angles <- new_cell_position$angle_to_bud_scar[
  new_cell_position$strain == 'grande' &
    !is.na(new_cell_position$angle_to_bud_scar)
]

petite_angles <- new_cell_position$angle_to_bud_scar[
  new_cell_position$strain == 'petite' &
    !is.na(new_cell_position$angle_to_bud_scar)
]

# Perform the two-sample KS test
ks_result <- ks.test(grande_angles, petite_angles)

# View the results
print(ks_result)

# adjusted p-value 0.05/4 = 0.0125
# all divisions: p-value = 0.0008095 (significant)
# first division: p-value = 0.003044 (significant)
# second division: p-value = 0.4739
# third division: p-value = 0.4288

table(new_cell_position[(new_cell_position$strain=='grande' | 
                           new_cell_position$strain=='petite') &
                          !is.na(new_cell_position$angle_to_bud_scar),]$strain,
      new_cell_position[(new_cell_position$strain=='grande' | 
                           new_cell_position$strain=='petite') &
                          !is.na(new_cell_position$angle_to_bud_scar),]$division_number)
#          1   2   3 total
# grande 155 120  21   296
# petite 125  99  19   243


# Create ECDF functions
ecdf_grande <- ecdf(grande_angles)
ecdf_petite <- ecdf(petite_angles)

# Create a fine grid of points to evaluate
all_points <- sort(unique(c(grande_angles, petite_angles)))

# Calculate the difference at each point
differences <- abs(ecdf_grande(all_points) - ecdf_petite(all_points))

# Find the maximum
max_diff_idx <- which.max(differences)
max_diff_value <- differences[max_diff_idx]
max_diff_location <- all_points[max_diff_idx]

cat("Maximum difference:", max_diff_value, "\n")
cat("Location (angle):", max_diff_location, "\n")



ggplot(new_cell_position[(new_cell_position$strain=='grande' | 
                            new_cell_position$strain=='petite') &
                           !is.na(new_cell_position$angle_to_bud_scar),],
       aes(x=strain, y=angle_to_bud_scar, fill=strain))+
  geom_violin()+
  scale_fill_manual(values = c("grande" = "#AE93BEFF", 
                               "petite" = "#F0D77BFF"))+
  stat_summary(fun="median", geom="crossbar")+
  labs(y = "Cumulative Proportion")+
  theme_classic()+
  # geom_hline(yintercept = 35.27)+
  NULL

wrs_result = wilcox.test(grande_angles, petite_angles)
print(wrs_result)

median(petite_angles) # 47.48933
median(grande_angles) # 35.27551

cliff.delta(grande_angles, petite_angles)
# delta estimate: -0.1106384 (negligible)

ggplot(new_cell_position[(new_cell_position$strain=='grande' | 
                            new_cell_position$strain=='petite') &
                           !is.na(new_cell_position$angle_to_bud_scar),],
       aes(x=angle_to_bud_scar, color=strain))+
  stat_ecdf(geom = "step", linewidth = 1)+
  labs(y = "Cumulative Proportion")+
  theme_classic()+
  facet_wrap(~division_number)+
  NULL



