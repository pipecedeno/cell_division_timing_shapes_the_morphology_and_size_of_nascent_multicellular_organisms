
#Date: 18 September 2026
# Code used to analyze the results of cell_division_synchronicity_simulation_TCD_edit_save_csv.m

library(plyr)  # load before dplyr
library(dplyr)
library(ggplot2)
library(stringi)
library(purrr)
library(stringr)
library(ggpubr)
library(tidyverse)
library(cowplot)


setwd("~/work_dir/observed_synchrony/paper_results_edge_degree_15_2025may29/corrections_supreme_version_sep2026/tom_burst_sim_supplementary/")


BINWIDTH <- 10
BINEDGES <- seq(0, 1000, by = BINWIDTH)


## Bin raw division-event times into 10-min windows, per group,
## filling empty bins with 0 (matching MATLAB's histcounts behavior),
## then replacing 0 with 1e-2 so log-scale plotting doesn't drop points
## (this mirrors `Timetrack(Timetrack==0) = 1e-2;` in the .m file).
bin_counts <- function(df, time_col, group_vars = character(0)) {
  bin_left_vals <- head(BINEDGES, -1)
  
  df <- df %>%
    mutate(bin_left = BINEDGES[findInterval(.data[[time_col]], BINEDGES,
                                            rightmost.closed = TRUE)])
  
  df %>%
    count(across(all_of(group_vars)), bin_left, name = "n") %>%
    group_by(across(all_of(group_vars))) %>%
    complete(bin_left = bin_left_vals, fill = list(n = 0)) %>%
    ungroup() %>%
    mutate(
      bin_center = bin_left + BINWIDTH / 2,
      n_plot     = ifelse(n == 0, 1e-2, n)
    )
}


# Simulation 1: default sampling ####

sim1_df=read.csv("sim1_divisions.csv", header = TRUE)
summary(sim1_df)

sim1_binned <- bin_counts(sim1_df, "divTime")

p_sim1=ggplot(sim1_binned[sim1_binned$n_plot!=0.01,], 
              aes(x = bin_center, y = n_plot)) +
  geom_col(fill = "grey", width = BINWIDTH, color = NA) +
  labs(x="Time (min)", y = "No. Cell Divisions")+
  coord_cartesian(xlim = c(0, 600)) +
  scale_y_log10() +
  theme_classic(base_size = 8) +
  NULL
p_sim1


# Simualtion 2: Variance ####

sim2_df <- read.csv("sim2_divisions.csv", header=TRUE)
summary(sim2_df)

sim2_binned <- bin_counts(sim2_df, "divTime", group_vars = "std_value")

std_labels <- as.character(unique(sim2_df$std_value))

p_sim2=ggplot(sim2_binned,
              aes(x = bin_center, y = n_plot,
                  color = factor(std_value, labels=std_labels))) +
  geom_line() +
  scale_color_grey(start = 0.7, end = 0.1) +
  labs(x="Time (min)", y = "No. Cell Divisions")+
  guides(col=guide_legend(title ="Standard\nDeviation"))+
  coord_cartesian(xlim = c(0, 600)) +
  scale_y_log10() +
  theme_classic(base_size = 8) +
  NULL
p_sim2


# Simulation 3: Delay ####

sim3_df <- read.csv("sim3_divisions.csv", header=TRUE)
sim3_binned <- bin_counts(sim3_df, "divTime", group_vars = "delay")
sim3_binned$delay_percentage=sim3_binned$delay*100/105

p_sim3=ggplot(sim3_binned,
              aes(x = bin_center, y = n_plot,
                  color = factor(delay_percentage))) +
  geom_line() +
  scale_color_brewer(palette = "Purples")+
  labs(x="Time (min)", y = "No. Cell Divisions")+
  guides(col=guide_legend(title ="Delay %"))+
  coord_cartesian(xlim = c(0, 600)) +
  scale_y_log10() +
  theme_classic(base_size = 8) +
  NULL
p_sim3


# Simulation 4: Combined effects ancestors ####

sim4_df <- read.csv("sim4_divisions.csv", header=TRUE)
sim4_binned <- bin_counts(sim4_df, "divTime", group_vars = "condition")

p_sim4=ggplot(sim4_binned,
              aes(x = bin_center, y = n_plot, color = condition)) +
  geom_line() +
  scale_color_manual(values = c(Grande = "#AE93BEFF", Petite = "#F0D77BFF")) +
  labs(x="Time (min)", y = "No. Cell Divisions")+
  guides(col=guide_legend(title ="Ancestor"))+
  coord_cartesian(xlim = c(0, 600)) +
  scale_y_log10() +
  theme_classic(base_size = 8) +
  NULL
p_sim4


# Supplementary Figure ####

supp_fig=plot_grid(p_sim1, p_sim2, p_sim3, p_sim4,
                   labels=LETTERS[1:4], label_size=10)
supp_fig


ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_2-temp_sim_burstiness_18sep2026.png',
       plot=supp_fig, dpi='retina', width=8, height=5.6, bg='white')

