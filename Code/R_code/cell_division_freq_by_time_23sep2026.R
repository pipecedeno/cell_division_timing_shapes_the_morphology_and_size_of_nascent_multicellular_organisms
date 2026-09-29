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
BINWIDTH_stats = 5
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
  geom_col(fill = "#6699CC", width = BINWIDTH, color = NA) +
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

sigma_labels <- as.character(unique(sim2_df$std_value))

p_sim2=ggplot(sim2_binned,
              aes(x = bin_center, y = n_plot,
                  color = factor(std_value, labels=sigma_labels))) +
  geom_line(linewidth = 1) +
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
  geom_line(linewidth = 1) +
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
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(Grande = "#AE93BEFF", Petite = "#F0D77BFF")) +
  labs(x="Time (min)", y = "No. Cell Divisions")+
  guides(col=guide_legend(title = ""))+
  coord_cartesian(xlim = c(0, 600)) +
  scale_y_log10() +
  theme_classic(base_size = 8) +
  NULL
p_sim4


# STATS: ACF-based synchrony index ####
# 
# Quantifies "burstiness" per simulation run using the autocorrelation
# function (ACF) of that run's per-capita division rate, then compares
# the resulting scalar across conditions.
#
# Design notes:
#   - Unit of replication = one simulated run (NSIMS=10 per condition),
#     not a lag bin within a single curve. So instead of comparing whole
#     ACF curves directly, each run is collapsed to one number first,
#     then that number is compared across conditions with a standard
#     between-group test. This avoids treating within-curve lag bins
#     (which aren't independent observations) as if they were a sample.
#   - Raw division counts grow ~exponentially with population size, which
#     would swamp any periodic signal in a raw ACF. Counts are converted
#     to a per-capita rate (divisions in a bin / cells alive at the start
#     of that bin) before computing the ACF.
#   - MAX_EVENTS in the .m script can make a run stop before SIM_RUNTIME.
#     common_window() finds the longest window still fully covered by
#     every run in a given dataset, so all runs/conditions are compared
#     on equal footing rather than on whatever length each happened to
#     reach.

N0          <- 30    # NCELLS_mo + NCELLS_da in the .m script
MU_EXPECTED <- 105   # nominal doubling time, minutes (`mu` in the .m script)

# Bin one condition's raw events *per run* (bin_counts() already supports
# this -- just fold "run" into the grouping variables instead of pooling
# runs together).
bin_counts_per_run <- function(df, time_col, group_vars) {
  bin_counts(df, time_col, group_vars = c(group_vars, "run"))
}

# Convert binned division counts into a per-capita division rate per bin,
# per run. pop_before = N0 + (all divisions strictly before this bin).
add_percapita_rate <- function(binned_df, group_vars) {
  binned_df %>%
    group_by(across(all_of(c(group_vars, "run")))) %>%
    arrange(bin_left, .by_group = TRUE) %>%
    mutate(
      pop_before = N0 + lag(cumsum(n), default = 0),
      rate       = n / pop_before
    ) %>%
    ungroup()
}

# Longest time window (min) for which every run in `df` has real simulated
# data, i.e. wasn't cut short by hitting MAX_EVENTS before SIM_RUNTIME.
common_window <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(c(group_vars, "run")))) %>%
    summarise(max_t = max(divTime), .groups = "drop") %>%
    pull(max_t) %>%
    min()
}

# ACF-based synchrony index for one run's per-capita rate series: the
# strongest autocorrelation found in a window of lags around the expected
# doubling time (0.5x-1.5x MU_EXPECTED), excluding lag 0. A synchronized
# run should show a clear secondary peak in this window; a desynchronized
# run should not.
synchrony_index <- function(rate_vec, binwidth = BINWIDTH_stats, mu_expected = MU_EXPECTED) {
  lag_min <- floor(0.5 * mu_expected / binwidth)
  lag_hi  <- ceiling(1.5 * mu_expected / binwidth)
  lag_hi  <- min(lag_hi, length(rate_vec) - 1)
  if (lag_hi <= lag_min || length(rate_vec) < 3) return(NA_real_)
  acf_out <- acf(rate_vec, lag.max = lag_hi, plot = FALSE)
  max(acf_out$acf[(lag_min + 1):(lag_hi + 1)], na.rm = TRUE)
}

# One synchrony index per run, for one condition's full dataset. Caps the
# analysis window at tmax_ceiling but will use less if the shortest run
# in the dataset doesn't reach that far.
compute_synchrony_table <- function(df, time_col, group_vars, tmax_ceiling = 700) {
  win <- min(tmax_ceiling, common_window(df, group_vars))
  message(sprintf("  Common analysis window: 0-%.0f min (capped by shortest run)", win))
  
  result <- bin_counts_per_run(df, time_col, group_vars) %>%
    add_percapita_rate(group_vars) %>%
    filter(bin_left < win) %>%
    group_by(across(all_of(c(group_vars, "run")))) %>%
    summarise(sync_index = synchrony_index(rate), .groups = "drop")
  
  n_na <- sum(is.na(result$sync_index))
  if (n_na > 0) {
    message(sprintf("  %d run(s) had too little data for a synchrony index and were dropped.", n_na))
  }
  result %>% filter(!is.na(sync_index))
}

# Variant of compute_synchrony_table() that allows a different expected
# doubling time per condition level. Needed for Sim 4: petite's
# steady-state period is set by mu1_petite = mu+DT for EVERY division
# (muFirst and muSecond both use it in the .m script), not just the first
# one like in Sim 3's delay sweep -- so grande and petite genuinely have
# different ground-truth periods baked into the generative model, and
# searching both at the baseline 105 min window would under-detect
# petite's periodicity.
#
# `group_vars` must be a single column name here (e.g. "condition"), since
# `mu_map` is a named vector/list keyed on that column's levels, e.g.
#   c(Grande = 105, Petite = 135)
compute_synchrony_table_mu_map <- function(df, time_col, group_vars, mu_map, tmax_ceiling = 700) {
  win <- min(tmax_ceiling, common_window(df, group_vars))
  message(sprintf("  Common analysis window: 0-%.0f min (capped by shortest run)", win))
  
  group_col <- group_vars[1]
  
  result <- bin_counts_per_run(df, time_col, group_vars) %>%
    add_percapita_rate(group_vars) %>%
    filter(bin_left < win) %>%
    group_by(across(all_of(c(group_vars, "run")))) %>%
    summarise(
      sync_index = synchrony_index(
        rate,
        mu_expected = mu_map[[as.character(first(.data[[group_col]]))]]
      ),
      .groups = "drop"
    )
  
  n_na <- sum(is.na(result$sync_index))
  if (n_na > 0) {
    message(sprintf("  %d run(s) had too little data for a synchrony index and were dropped.", n_na))
  }
  result %>% filter(!is.na(sync_index))
}

# --- Diagnostic: plot one run's raw division counts vs. its normalized
# per-capita rate, stacked, to check that the population-growth de-trending
# step (add_percapita_rate()) is doing what it's supposed to: the top panel
# (raw n) should show counts growing over time as the population grows;
# the bottom panel (rate) should look roughly flat/stationary, with any
# periodic bumps from synchronized divisions now visible instead of being
# swamped by the growth trend.
#
#   df            - raw division-event dataframe (sim2_df, sim3_df, or sim4_df)
#   group_var     - name of the category column to filter on, e.g. "std_value", "delay", "condition"
#   group_value   - value of that category to look at, e.g. 20, 10, or "Petite"
#   run_number    - which replicate run to plot (matches the "run" column, 1-NSIMS)
#   mu_expected   - expected doubling time (min), used to draw dashed reference
#                   lines every mu_expected minutes. For Sim 4, pass the
#                   condition-specific value, e.g. sim4_mu_map[["Petite"]].
#
# Example:
#   plot_run_rate(sim2_df, "std_value", 20, run_number = 3)
#   plot_run_rate(sim4_df, "condition", "Petite", run_number = 1,
#                 mu_expected = sim4_mu_map[["Petite"]])
plot_run_rate <- function(df, group_var, group_value, run_number,
                          time_col = "divTime", mu_expected = MU_EXPECTED,
                          tmax_ceiling = 700) {
  
  win <- min(tmax_ceiling, common_window(df, group_var))
  
  rate_df <- bin_counts_per_run(df, time_col, group_var) %>%
    add_percapita_rate(group_var) %>%
    filter(.data[[group_var]] == group_value, run == run_number, bin_left < win) %>%
    arrange(bin_left)
  
  if (nrow(rate_df) < 2) {
    stop(sprintf("Not enough data for %s == %s, run %s (only %d bin(s) after truncating to 0-%.0f min).",
                 group_var, group_value, run_number, nrow(rate_df), win))
  }
  
  op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  on.exit(par(op))
  
  ref_lines <- seq(mu_expected, win, by = mu_expected)
  
  # Top: raw division counts per bin -- expect this to trend upward as the
  # population grows.
  plot(rate_df$bin_left, rate_df$n, type = "l",
       xlab = "Time (min)", ylab = "Divisions per bin (n)",
       main = sprintf("%s = %s, run %s -- raw counts", group_var, group_value, run_number))
  abline(v = ref_lines, col = "red", lty = 2)
  
  # Bottom: per-capita rate (n / population at start of bin) -- should look
  # roughly stationary if de-trending worked, with periodic bumps visible.
  plot(rate_df$bin_left, rate_df$rate, type = "l",
       xlab = "Time (min)", ylab = "Normalized rate (n / pop_before)",
       main = "normalized bin data")
  abline(v = ref_lines, col = "red", lty = 2)
  
  invisible(rate_df)
}

# --- Diagnostic: plot the ACF of one specific run, to eyeball what the
# synchrony_index()/generation_peaks() functions above are actually seeing.
# Works on any of the raw per-event dataframes (sim2_df, sim3_df, sim4_df)
# by specifying which category column/value and run number to look at.
#
#   df            - raw division-event dataframe (sim2_df, sim3_df, or sim4_df)
#   group_var     - name of the category column to filter on, e.g. "std_value", "delay", "condition"
#   group_value   - value of that category to look at, e.g. 20, 10, or "Petite"
#   run_number    - which replicate run to plot (matches the "run" column, 1-NSIMS)
#   mu_expected   - expected doubling time (min) to mark with dashed lines. For
#                   Sim 4, pass the condition-specific value, e.g.
#                   sim4_mu_map[["Petite"]], since grande/petite differ.
#   n_generations - how many generation markers (mu, 2*mu, 3*mu, ...) to draw
#
# Example:
#   plot_run_acf(sim2_df, "std_value", 20, run_number = 3)
#   plot_run_acf(sim4_df, "condition", "Petite", run_number = 1,
#                mu_expected = sim4_mu_map[["Petite"]])
plot_run_acf <- function(df, group_var, group_value, run_number,
                         time_col = "divTime", mu_expected = MU_EXPECTED,
                         tmax_ceiling = 700, n_generations = 3) {
  
  win <- min(tmax_ceiling, common_window(df, group_var))
  
  rate_df <- bin_counts_per_run(df, time_col, group_var) %>%
    add_percapita_rate(group_var) %>%
    filter(.data[[group_var]] == group_value, run == run_number, bin_left < win) %>%
    arrange(bin_left)
  
  if (nrow(rate_df) < 3) {
    stop(sprintf("Not enough data for %s == %s, run %s (only %d bin(s) after truncating to 0-%.0f min).",
                 group_var, group_value, run_number, nrow(rate_df), win))
  }
  
  lag_max <- min(ceiling((n_generations + 0.5) * mu_expected / BINWIDTH), nrow(rate_df) - 1)
  
  acf_out <- acf(rate_df$rate, lag.max = lag_max, plot = TRUE,
                 main = sprintf("%s = %s, run %s", group_var, group_value, run_number),
                 xlab = "Lag (bins)")
  
  # Dashed reference lines at 1x, 2x, 3x... the expected doubling time, in
  # bin units, so you can see whether the ACF peaks line up with them.
  abline(v = seq_len(n_generations) * mu_expected / BINWIDTH, col = "red", lty = 2)
  
  invisible(acf_out)
}

#### Example plots of binned normalized data ####

plot_run_rate(sim2_df, "std_value", 20, run_number = 3)

plot_run_rate(sim3_df, "delay", 13.125, run_number = 1)

plot_run_rate(sim4_df, "condition", "Petite", run_number = 1,
             mu_expected = sim4_mu_map[["Petite"]])

plot_run_rate(sim4_df, "condition", "Grande", run_number = 1,
             mu_expected = sim4_mu_map[["Grande"]])



#### Example plots ACF values ####

plot_run_acf(sim2_df, "std_value", 20, run_number = 3)

unique(sim3_df$delay)
plot_run_acf(sim3_df, "delay", 13.125, run_number = 1)


sim4_mu_map <- c(Grande = MU_EXPECTED, Petite = 120)

plot_run_acf(sim4_df, "condition", "Petite", run_number = 1,
             mu_expected = sim4_mu_map[["Petite"]])

plot_run_acf(sim4_df, "condition", "Grande", run_number = 1,
             mu_expected = sim4_mu_map[["Grande"]])



#### Simulation 2: variance sweep ####
sim2_sync <- compute_synchrony_table(sim2_df, "divTime", "std_value")

print(cor.test(sim2_sync$std_value, sim2_sync$sync_index, method = "spearman"))
# p-value < 2.2e-16
# rho -0.7778512 

p_sim2_sync <- ggplot(sim2_sync, aes(x = factor(std_value), y = sync_index)) +
  geom_boxplot(fill = "grey85", outlier.shape = NA) +
  geom_jitter(width = 0.1, color = "#6699CC", size = 1.5) +
  labs(x = "Standard deviation (min)", y = "ACF synchrony index") +
  theme_classic(base_size = 8)
p_sim2_sync

#### Simulation 3: delays ####
sim3_sync <- compute_synchrony_table(sim3_df, "divTime", "delay")

print(cor.test(sim3_sync$delay, sim3_sync$sync_index, method = "spearman"))
# p-value < 2.2e-16
# rho -0.8770457

p_sim3_sync <- ggplot(sim3_sync, aes(x = factor(delay*100/105), y = sync_index)) +
  geom_boxplot(fill = "grey85", outlier.shape = NA) +
  geom_jitter(width = 0.1, color = "#9966CC", size = 1.5) +
  labs(x = "Delay %", y = "ACF synchrony index") +
  theme_classic(base_size = 8)
p_sim3_sync

#### Simulation 4: grande vs. petite ####

sim4_mu_map <- c(Grande = MU_EXPECTED, Petite = 120)

sim4_sync <- compute_synchrony_table_mu_map(sim4_df, "divTime", "condition", mu_map = sim4_mu_map)

print(wilcox.test(sync_index ~ condition, data = sim4_sync))
# p-value < 2.2e-16

sim4_sync %>%
  group_by(condition) %>%
  summarize(median_acf=median(sync_index))

p_sim4_sync <- ggplot(sim4_sync, aes(x = condition, y = sync_index, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5) +
  scale_color_manual(values = c(Grande = "#AE93BEFF", Petite = "#F0D77BFF")) +
  labs(x = "Genotype", y = "ACF synchrony index") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")
p_sim4_sync

## --- Combined synchrony-index summary plot ------------------------------
supp_fig_sync <- plot_grid(p_sim2_sync, p_sim3_sync, p_sim4_sync,
                           labels = LETTERS[1:3], label_size = 10, nrow = 1)
supp_fig_sync

ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_3-temp_sim_burstiness_18sep2026_ACF_dist.png',
       plot=supp_fig_sync, dpi='retina', width=8, height=5.6, bg='white')



# Supplementary Figure ####

supp_fig=plot_grid(p_sim1, p_sim2, p_sim3, p_sim4,
                   labels=LETTERS[1:4], label_size=10)
supp_fig


ggsave(filename='~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Paper_figures/supp_fig_3-temp_sim_burstiness_18sep2026.png',
       plot=supp_fig, dpi='retina', width=8, height=5.6, bg='white')


# ACF for empirical data ####

df_all=read.csv("~/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/timelapse_doubling_time_inf_2025may5.csv", header = TRUE)

unique(df_all[df_all$strain=='grande',]$id_file)
# "2024nov5_gob8_4"   "2022nov22_gob8_6"  "2024oct31_gob8_2"  "2024oct30_gob8_11" "2022nov22_gob8_7"  "2024nov5_gob8_12" "2022nov23_gob8_7"  "2024oct30_gob8_6"

df_all=df_all[! df_all$id_file %in% c("2022nov22_gob8_6", "2022nov22_gob8_7", "2022nov23_gob8_7"),]

df_all$division_number=as.factor(df_all$division_number)
df_all$timepoint=factor(df_all$timepoint, levels=c('t0','t200','t400','t600','t800','t1000'))

summary(df_all)


doubled=df_all[df_all$N_PREDECESSORS==1 & df_all$N_SUCCESSORS==2,]

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


ancestor_tl_dt=normal_frames[normal_frames$strain=='petite' | normal_frames$strain=='grande',]
table(ancestor_tl_dt$strain)
ancestor_tl_dt$strain=factor(ancestor_tl_dt$strain, levels=c("grande", "petite"))

# Count divisions per existing 5-min bin, per group, filling in bins with
# zero divisions (which simply won't appear as rows in the raw event data).
bin_counts_empirical <- function(df, time_col, group_vars, binwidth = BINWIDTH_stats) {
  df <- df %>% mutate(bin_left = .data[[time_col]])
  
  df %>%
    count(across(all_of(group_vars)), bin_left, name = "n") %>%
    group_by(across(all_of(group_vars))) %>%
    complete(bin_left = seq(min(bin_left), max(bin_left), by = binwidth),
             fill = list(n = 0)) %>%
    ungroup()
}

# Longest time window (min) covered by every group's data. Analogous to
# common_window(), but grouped on strain only (no run column).
common_window_empirical <- function(df, group_vars, time_col) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(max_t = max(.data[[time_col]]), .groups = "drop") %>%
    pull(max_t) %>%
    min()
}

# Look up a per-group value that may be supplied either as one shared
# number or as a named vector/list keyed on the group level (e.g.
# c(WT = 105, mutant = 130)) -- lets mu_expected differ by strain, the
# same way compute_synchrony_table_mu_map() does for Sim 4.
lookup_per_group <- function(value, group_level) {
  if (length(value) == 1 && is.null(names(value))) return(value)
  value[[as.character(group_level)]]
}

# Per-capita division rate per bin, per group, using each group's own
# first-bin division count as its population-size proxy (n0) -- see note
# above. Analogous to add_percapita_rate(), but keyed on group_vars only
# (no run) and with n0 derived from the data instead of a fixed constant.
add_percapita_rate_empirical <- function(binned_df, group_vars) {
  binned_df %>%
    group_by(across(all_of(group_vars))) %>%
    arrange(bin_left, .by_group = TRUE) %>%
    mutate(
      n0         = first(n),
      pop_before = n0 + lag(cumsum(n), default = 0),
      rate       = n / pop_before
    ) %>%
    ungroup() %>%
    dplyr::select(-n0)
}

# One synchrony index per group (strain) for the empirical dataset.
# Returns one row per group_vars level -- NOT one row per run, since
# empirical data has no run/replicate structure to summarise across.
compute_synchrony_table_empirical <- function(df, time_col = "FRAME_minutes", group_vars = "strain",
                                              mu_expected = MU_EXPECTED,
                                              binwidth = BINWIDTH_stats, tmax_ceiling = 700) {
  win <- min(tmax_ceiling, common_window_empirical(df, group_vars, time_col))
  message(sprintf("  Common analysis window: 0-%.0f min (capped by shorter group)", win))
  
  group_col <- group_vars[1]
  
  result <- bin_counts_empirical(df, time_col, group_vars, binwidth = binwidth) %>%
    add_percapita_rate_empirical(group_vars) %>%
    filter(bin_left < win) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      sync_index = synchrony_index(
        rate,
        binwidth    = binwidth,
        mu_expected = lookup_per_group(mu_expected, first(.data[[group_col]]))
      ),
      .groups = "drop"
    )
  
  n_na <- sum(is.na(result$sync_index))
  if (n_na > 0) {
    message(sprintf("  %d group(s) had too little data for a synchrony index and were dropped.", n_na))
  }
  result %>% filter(!is.na(sync_index))
}

ancestor_tl_dt$FRAME_minutes=ancestor_tl_dt$FRAME*5
compute_synchrony_table_empirical(ancestor_tl_dt)
# strain sync_index
# 1 grande     0.0898
# 2 petite     0.0514


