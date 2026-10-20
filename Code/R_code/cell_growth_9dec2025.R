

library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(effsize)
library(tidyr)
library(purrr)
library(broom)
library(ggh4x)

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

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_1_cell_growth_properties_28aug2026.png',
       plot=p_cell_props_v2, dpi='retina', width=6.5, height=3, bg='white')

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



# Correlation time lapse duration vs dt ####

df_all=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/timelapse_doubling_time_inf_2025may5.csv", header = TRUE)

unique(df_all[df_all$strain=='grande',]$id_file)
# "2024nov5_gob8_4"   "2022nov22_gob8_6"  "2024oct31_gob8_2"  "2024oct30_gob8_11" "2022nov22_gob8_7"  "2024nov5_gob8_12" "2022nov23_gob8_7"  "2024oct30_gob8_6"

df_all=df_all[! df_all$id_file %in% c("2022nov22_gob8_6", "2022nov22_gob8_7", "2022nov23_gob8_7"),]

df_all$division_number=as.factor(df_all$division_number)
df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

summary(df_all)

doubled=df_all[df_all$N_PREDECESSORS==1 & df_all$N_SUCCESSORS==2,]

ggplot(doubled[doubled$strain=='petite',],
       aes(x=FRAME, y=hours))+
  geom_point()+
  facet_wrap(~complete_id)+
  # geom_smooth(method = "lm", formula = y ~ x)+
  theme_classic()+
  NULL

doubled_correlation=doubled[doubled$id_file!='2022nov23_gob385_2',] %>% #
  mutate(FRAME_hours=(FRAME*5/60)) %>%
  group_by(complete_id, timepoint, strain, strain_timepoint) %>%
    summarise(correlation=cor(FRAME, hours),
            slope=unname(coef(glm(hours~FRAME_hours))[2]),
            p_val=summary(glm(hours~FRAME_hours))$coefficients[2, 4])

sum(doubled_correlation$p_val<(0.05/125))
doubled_correlation[doubled_correlation$p_val<(0.05/125),]$complete_id
# [1] "PA1_t600_14_2023jan24" "PA1_t600_5_2023jan24"  "PA2_t1000_2_2022nov22" "PA2_t1000_5_2022nov22" "PA2_t400_2_2023jan13" 
# [6] "PA2_t600_10_2023jan26" "PA4_t1000_3_2022dec1"  "PA4_t1000_8_2022dec1"  "PA4_t600_10_2023jan26" "PA4_t600_2_2023jan24" 
# [11] "PA4_t600_3_2023jan26"  "PA4_t600_5_2023jan26"  "PA5_t400_8_2023jan13"  "PA5_t600_15_2023jan26" "PA5_t600_7_2023jan26"
sum(doubled_correlation[doubled_correlation$p_val<(0.05/125),]$slope<0)

ggplot(doubled[doubled$complete_id %in% doubled_correlation[doubled_correlation$p_val<(0.05/125),]$complete_id,],
       aes(x=FRAME, y=hours, col=timepoint))+
  geom_point()+
  facet_wrap(~complete_id)+
  # geom_smooth(method = "lm", formula = y ~ x)+
  theme_classic()+
  scale_color_manual(values= c("#D8AF39FF", "#DE7862FF",  "#E75B64FF"))+
  NULL


ggplot(doubled_correlation, aes(x=strain, y=correlation, fill=timepoint))+
  # geom_violin()+
  geom_boxplot(alpha=0.75)+
  geom_jitter(width = 0, height = 0, col='grey')+
  facet_wrap(~timepoint, ncol=5, scales="free_x") +
  theme_classic()+
  guides(fill='none')+
  force_panelsizes(cols = c(2, 5, 5, 5, 5)) +
  scale_fill_manual(values= c("#AE93BEFF", "#E8C4A2FF", "#D8AF39FF", "#DE7862FF",  "#E75B64FF"))+
  NULL


ggplot(doubled_correlation, 
       aes(x=strain, y=slope, fill=timepoint))+
  # geom_violin()+
  geom_boxplot()+
  geom_jitter(width = 0, height = 0, col='grey')+
  facet_wrap(~timepoint, ncol=5, scales="free_x") +
  force_panelsizes(cols = c(2, 5, 5, 5, 5)) +
  theme_classic()+
  guides(fill='none')+
  geom_hline(yintercept = 0, linetype='dashed')+
  scale_fill_manual(values= c("#AE93BEFF", "#E8C4A2FF", "#D8AF39FF", "#DE7862FF",  "#E75B64FF"))+
  NULL


doubled_test <- doubled_correlation %>%
  group_by(strain_timepoint) %>%
  summarise(
    slope_test = list(tidy(wilcox.test(slope, mu = 0))),
    corr_test  = list(tidy(wilcox.test(correlation, mu = 0))),
    .groups = "drop"
  ) %>%
  unnest(slope_test, names_sep = "_") %>%
  unnest(corr_test, names_sep = "_")


sum(doubled_test$slope_test_p.value<(0.05/22))
sum(doubled_test$corr_test_p.value<(0.05/22))



for (i in unique(doubled$strain)){
  print(ggplot(doubled[doubled$strain==i,], aes(x=complete_id, y=hours, fill=timepoint))+
          geom_violin()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          ggtitle(paste("Individual time lapses for",i))+
          NULL)
}

ggplot(doubled[doubled$strain=='PA4' & doubled$timepoint=='t200',], aes(x=complete_id, y=hours, fill=timepoint))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(paste("Individual time lapses for",i))+
  NULL



ggplot(doubled[doubled$strain=='grande',],
       aes(x=FRAME, y=hours))+
  geom_point()+
  facet_wrap(~complete_id)+
  # geom_smooth(method = "lm", formula = y ~ x)+
  theme_classic()+
  NULL

time_vs_dt_data=doubled[doubled$complete_id %in% c("grande_t0_4_2024nov5", "grande_t0_6_2024oct30",
                                                   "petite_t0_2_2022sep30", "petite_t0_4_2022sep27"),] %>%
  group_by(strain) %>%
  mutate(pair_number = dense_rank(complete_id)) %>%
  ungroup()


ggplot(time_vs_dt_data,
       aes(x=FRAME, y=hours, col=strain))+
  geom_point(alpha=0.75)+
  facet_grid(pair_number~strain)+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col='black', linetype='dashed')+
  theme_classic()+
  scale_color_manual(values = c("#AE93BEFF", "#F0D77BFF"))+
  guides(col='none')+
  NULL

summary(glm(hours~FRAME, data=time_vs_dt_data[time_vs_dt_data$strain=='petite' & time_vs_dt_data$pair_number==2,]))
# p-value of 0.378


#### Permutation test ####



set.seed(123)          # reproducibility
n_perm <- 10000

df <- doubled[doubled$id_file != '2022nov23_gob385_2', ] %>%
  mutate(FRAME_hours = FRAME * 5 / 60)

## 1. Observed statistics per group
observed <- df %>%
  group_by(complete_id, timepoint, strain, strain_timepoint) %>%
  summarise(
    slope       = unname(coef(glm(hours ~ FRAME_hours))[2]),
    correlation = cor(FRAME_hours, hours),
    .groups = "drop"
  )

## 2. Permutations per group
perm_results <- df %>%
  group_by(complete_id, timepoint, strain, strain_timepoint) %>%
  group_modify(~ {
    x <- .x$FRAME_hours
    y <- .x$hours
    sx <- sd(x)
    
    map_dfr(seq_len(n_perm), function(i) {
      y_perm <- sample(y)                    # shuffle hours
      r <- cor(x, y_perm)
      tibble(
        perm        = i,
        correlation = r,
        slope       = r * sd(y_perm) / sx
      )
    })
  }) %>%
  ungroup()

## 3. (Optional) merge with observed to get empirical p-values
perm_pvals <- perm_results %>%
  left_join(observed, by = c("complete_id", "timepoint", "strain", "strain_timepoint"),
            suffix = c("_perm", "_obs")) %>%
  group_by(complete_id, timepoint, strain, strain_timepoint) %>%
  summarise(
    slope_obs        = first(slope_obs),
    correlation_obs  = first(correlation_obs),
    p_slope          = mean(abs(slope_perm) >= abs(slope_obs)),
    p_correlation    = mean(abs(correlation_perm) >= abs(correlation_obs)),
    .groups = "drop"
  )

# write.csv(perm_results, "~/perm_results_9sep2026_10k_permutations.csv", row.names = FALSE)
# write.csv(perm_pvals, "~/perm_pvals_9sep2026_10k_permutations.csv", row.names = FALSE)

sum(perm_pvals$p_slope<(0.05/125)) # 15
sum(perm_pvals$p_correlation<(0.05/125)) # 15

perm_pvals[perm_pvals$p_slope<(0.05/125),]$complete_id
# [1] "PA1_t600_14_2023jan24" "PA1_t600_5_2023jan24"  "PA2_t1000_2_2022nov22" "PA2_t1000_5_2022nov22" "PA2_t600_10_2023jan26"
# [6] "PA4_t1000_3_2022dec1"  "PA4_t1000_8_2022dec1"  "PA4_t400_14_2023jan12" "PA4_t600_10_2023jan26" "PA4_t600_2_2023jan24" 
# [11] "PA4_t600_3_2023jan26"  "PA4_t600_5_2023jan26"  "PA5_t400_8_2023jan13"  "PA5_t600_15_2023jan26" "PA5_t600_7_2023jan26" 

perm_pvals[perm_pvals$p_correlation<(0.05/125),]$complete_id


# Plot distribution of permutations for significant values
sig_groups <- perm_pvals %>%
  filter(p_slope < (0.05 / 125)) %>%
  select(complete_id, timepoint, strain, strain_timepoint, slope_obs, p_slope)

sig_perm_slopes <- perm_results %>%
  semi_join(sig_groups, by = c("complete_id", "timepoint", "strain", "strain_timepoint")) %>%
  left_join(sig_groups, by = c("complete_id", "timepoint", "strain", "strain_timepoint"))

ggplot(sig_perm_slopes, aes(x = slope)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(
    data = sig_groups,
    aes(xintercept = slope_obs),
    color = "red", linewidth = 1
  ) +
  facet_wrap(~ complete_id, scales = "free") +
  labs(
    x = "Permuted slope",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 8))


# Plot petite permutation results
sig_groups <- perm_pvals %>%
  filter(strain == 'petite') %>%
  select(complete_id, timepoint, strain, strain_timepoint, slope_obs, p_slope)

sig_perm_slopes <- perm_results %>%
  semi_join(sig_groups, by = c("complete_id", "timepoint", "strain", "strain_timepoint")) %>%
  left_join(sig_groups, by = c("complete_id", "timepoint", "strain", "strain_timepoint"))

ggplot(sig_perm_slopes, aes(x = slope)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(
    data = sig_groups,
    aes(xintercept = slope_obs),
    color = "red", linewidth = 1
  ) +
  facet_wrap(~ complete_id, scales = "free") +
  labs(
    x = "Permuted slope",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 8)) +
  NULL




#### Updated supplementary figure 1 ####

p_t_vs_dt=ggplot(time_vs_dt_data,
       aes(x=FRAME*5/60, y=hours, col=strain))+
  geom_point(alpha=0.75)+
  facet_grid(pair_number~strain, labeller = labeller(strain = c("grande" = "Grande", "petite" = "Petite")))+
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, col='black', linetype='dashed')+
  theme_classic(base_size = 8)+
  scale_color_manual(values = c("#AE93BEFF", "#F0D77BFF"))+
  guides(col='none')+
  labs(y="Cell Doubling time (hours)", x="Time lapse duration (hours)")+
  NULL
p_t_vs_dt

doubled_correlation_plot <- doubled_correlation %>%
  mutate(fill_group = case_when(
    timepoint == "t0" & strain == "petite" ~ "t0_petite",
    timepoint == "t0" & strain == "grande" ~ "t0_grande",
    TRUE ~ as.character(timepoint)
  ))

p_correlation <- ggplot(doubled_correlation_plot, aes(x = strain, y = correlation, fill = fill_group)) +
  # y=correlation or slope
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0, height = 0, size = 1, col = 'grey', alpha = 0.75) + # 
  facet_wrap(~timepoint, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 8) +
  guides(fill = 'none') +
  force_panelsizes(cols = c(2, 5, 5, 5, 5)) +
  scale_fill_manual(values = c(
    "t0_petite" = "#F0D77BFF",   # your custom color for petite at t0
    "t0_grande" = "#AE93BEFF",   # your custom color for grande at t0
    "t200" = "#E8C4A2FF",
    "t400" = "#D8AF39FF",
    "t600" = "#DE7862FF",
    "t1000" = "#E75B64FF"
  )) +
  scale_x_discrete(labels = function(x) {
    x[x == "petite"] <- "P"
    x[x == "grande"] <- "G"
    x
  }) +
  labs(x = "Strain", y = "Correlation") +
  NULL
p_correlation


perm_pvals_plot <- perm_pvals %>%
  mutate(fill_group = case_when(
    timepoint == "t0" & strain == "petite" ~ "t0_petite",
    timepoint == "t0" & strain == "grande" ~ "t0_grande",
    TRUE ~ as.character(timepoint)
  ))

ggplot(perm_pvals_plot, aes(x = strain, y = slope_obs, fill = fill_group)) +
  # y=correlation or slope
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  geom_jitter(data=perm_pvals_plot, width = 0, height = 0, size = 2, aes(col=p_slope<(0.05/125))) + # col = 'grey', alpha = 0.75
  facet_wrap(~timepoint, ncol = 5, scales = "free_x") +
  theme_classic(base_size = 8) +
  guides(fill = 'none') +
  force_panelsizes(cols = c(2, 5, 5, 5, 5)) +
  scale_fill_manual(values = c(
    "t0_petite" = "#F0D77BFF",   # your custom color for petite at t0
    "t0_grande" = "#AE93BEFF",   # your custom color for grande at t0
    "t200" = "#E8C4A2FF",
    "t400" = "#D8AF39FF",
    "t600" = "#DE7862FF",
    "t1000" = "#E75B64FF"
  )) +
  scale_x_discrete(labels = function(x) {
    x[x == "petite"] <- "P"
    x[x == "grande"] <- "G"
    x
  }) +
  labs(x = "Strain", y = "Slope") +
  NULL


supp_fig_1_row1=plot_grid(p_dt_asp_r, p_dt_diam,
                          ncol=2, labels=LETTERS[1:2], label_size=10)
supp_fig_1_row1

supp_fig_1_merged=plot_grid(supp_fig_1_row1, p_t_vs_dt, p_correlation,
                            ncol=1, labels=c("", "C", "D"), label_size=10)
supp_fig_1_merged

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_1_cell_growth_properties_9sep2026.png',
       plot=supp_fig_1_merged, dpi='retina', width=6.5, height=9, bg='white')










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



