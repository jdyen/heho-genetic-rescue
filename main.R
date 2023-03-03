# Analysis of helmeted honeyeater fitness under genetic rescue

# load packages
library(qs)
library(rstanarm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(patchwork)

# load some helper functions
source("R/utils.R")

# flag to refit models (load saved models otherwise)
refit_models <- FALSE

# load data
breeding_data <- read.csv("data/TableS4_nests.csv")
nesting_data <- read.csv("data/TableS3_first_time_pairs.csv")

# replacing NAs with "UNKNOWN" for FOSTERED_EGG and FOSTERED_EGG_CHICK so model will work
breeding_data <- breeding_data %>%
  mutate(
    FOSTERED_EGG = ifelse(is.na(FOSTERED_EGG), "UNKNOWN", FOSTERED_EGG),
    FOSTERED_EGG_CHICK = ifelse(is.na(FOSTERED_EGG_CHICK), "UNKNOWN", FOSTERED_EGG_CHICK)
  )

# MCMC settings
seed <- 35151121
num_chains <- 2
num_iter <- 50
num_cores <- num_chains

# refit models if required
models_exist <- grepl("mod-", dir("outputs/fitted/"))
if (refit_models | !any(models_exist)) {
  
  # fit poisson regression models for clutch size data
  mod_clutchsize <- stan_glmer(Clutch_size  ~ SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
                                 (1 | SIRE_ORIGIN_AGE) + 
                                 (1 | Breeding_Season)
                               + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire), 
                               data = breeding_data,
                               family = stats::poisson,
                               prior = NULL,
                               seed = seed,
                               adapt_delta = 0.9, 
                               chains = num_chains, 
                               iter = num_iter, 
                               cores = num_cores
  )
  
  # save fitted model
  qsave(mod_clutchsize, file = "outputs/fitted/mod-clutchsize.qs")
  
  # fit binomial model for count of fertile eggs relative to all eggs
  #   with known fertility
  mod_eggsfertile <- stan_glmer(
    cbind(N_CONFIRMED_FERTILE_EGGS,
          N_EGGS_KNOWN_FERTILITY - N_CONFIRMED_FERTILE_EGGS)  ~
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire),
    data = breeding_data,
    family = stats::binomial,
    prior = NULL,
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(mod_eggsfertile, file = "outputs/fitted/mod-eggsfertile.qs")
  
  # fit binomial model for count of hatched eggs relative to clutch size
  mod_hatchedeggs <- stan_glmer(
    cbind(Hatched_eggs, Clutch_size - Hatched_eggs) ~ 
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire)
    + (1 | FOSTERED_EGG), 
    data = breeding_data,
    family = stats::binomial,
    prior = NULL,
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(mod_hatchedeggs, file = "outputs/fitted/mod-hatchedeggs.qs")
  
  # fit binomial model for count of fledged chicks relative to hatched eggs
  mod_chicksfledged <- stan_glmer(
    cbind(Fledged_chicks, Hatched_eggs - Fledged_chicks) ~
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire)
    + (1 | FOSTERED_EGG_CHICK), 
    data = breeding_data, 
    family = stats::binomial,
    prior = NULL,
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(mod_chicksfledged, file = "outputs/fitted/mod-chicksfledged.qs")
  
  # fit binomial model for count of independent chicks relative
  #   to fledged chicks
  mod_independentchicks <- stan_glmer(
    cbind(Indep_chicks, Fledged_chicks - Indep_chicks) ~
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire)
    + (1 | FOSTERED_EGG_CHICK),
    data = breeding_data, 
    family = stats::binomial,
    prior = NULL,
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(
    mod_independentchicks,
    file = "outputs/fitted/mod-independentchicks.qs"
  )
  
  # fit binomial model for count of male independent chicks 
  #   relative to independent chicks
  mod_male_independentchicks <- stan_glmer(
    cbind(MALE_INDEPENDENT, FEMALE_INDEPENDENT) ~
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID) + (1 | Dam_Sire)
    + (1 | FOSTERED_EGG_CHICK),
    data = breeding_data, 
    family = stats::binomial,
    prior = NULL,
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(
    mod_male_independentchicks, 
    file = "outputs/fitted/mod-male-independentchicks.qs"
  )
  
  # fit logistic regression (Bernoulli model) to data on nest building success
  nesting_data <- nesting_data %>%
    mutate(NEST_BUILT_YES = ifelse(NEST_BUILT == "Built", 1, 0),
           NEST_BUILT_NO = ifelse(NEST_BUILT == "Not_built", 1, 0))
  mod_nestbuilt <- stan_glmer(
    cbind(NEST_BUILT_YES, NEST_BUILT_NO) ~
      SEX_SPECIFIC_PAIR_TYPE + (1 | DAM_ORIGIN_AGE) + 
      (1 | SIRE_ORIGIN_AGE) +
      (1 | Breeding_Season)
    + (1 | Sire_ID) + (1 | Dam_ID)
    + Sire_pairing_for_season +
      Dam_pairing_for_season,
    data = nesting_data, 
    family = stats::binomial,
    prior = NULL, 
    seed = seed,
    adapt_delta = 0.9, 
    chains = num_chains, 
    iter = num_iter, 
    cores = num_cores
  )
  
  # save fitted model
  qsave(mod_nestbuilt, file = "outputs/fitted/mod-nestbuilt.qs")
  
} else {
  
  # load saved models from files otherwise
  mod_clutchsize <- qread("outputs/fitted/mod-clutchsize.qs")
  mod_eggsfertile <- qread("outputs/fitted/mod-eggsfertile.qs")
  mod_hatchedeggs <- qread("outputs/fitted/mod-hatchedeggs.qs")
  mod_chicksfledged <- qread("outputs/fitted/mod-chicksfledged.qs")
  mod_independentchicks <- qread("outputs/fitted/mod-independentchicks.qs")
  mod_male_independentchicks <- qread("outputs/fitted/mod-male-independentchicks.qs")
  mod_nestbuilt <- qread("outputs/fitted/mod-nestbuilt.qs")
  
}

# extract and plot outcomes on the link scale (relative to CC cross type)
p1 <- plot_parameters(
  mod_nestbuilt,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "A. nest built"
)

p2 <- plot_parameters(
  mod_nestbuilt,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "A. nest built"
)

p3 <- plot_parameters(
  mod_clutchsize,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "B. clutch size"
)

p4 <- plot_parameters(
  mod_eggsfertile,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "C. fertile eggs"
)

p5 <- plot_parameters(
  mod_eggsfertile,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1"),
  title = "C. fertile eggs"
)

p6 <- plot_parameters(
  mod_hatchedeggs,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "D. hatched eggs"
)

p7 <- plot_parameters(
  mod_chicksfledged,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "E. chicks fledged"
)

p8 <- plot_parameters(
  mod_chicksfledged,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "E. chicks fledged"
)

p9 <- plot_parameters(
  mod_independentchicks,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "F. independent chicks"
)

p10 <- plot_parameters(
  mod_independentchicks,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfGm"),
  title = "F. independent chicks"
)

p11 <- plot_parameters(
  mod_male_independentchicks,
  pars = c("SEX_SPECIFIC_PAIR_TYPECfF1m", 
           "SEX_SPECIFIC_PAIR_TYPECfGm",
           "SEX_SPECIFIC_PAIR_TYPEF1F1",
           "SEX_SPECIFIC_PAIR_TYPEF1fCm",
           "SEX_SPECIFIC_PAIR_TYPEGfCm"),
  title = "G. male independent chicks"
)

# combine these into a single plot
fig1_link_scale <- (p1 | p2 ) /
  (p3 | plot_spacer()) /
  (p4 | p5) / 
  (p6 | plot_spacer()) /
  (p7 | p8 ) /
  (p9 | p10) /
  (p11 | plot_spacer()) 

ggsave(
  file = "outputs/figures/Fig1-link-scale.png",
  plot = fig1_link_scale,
  device = png,
  height = 12,
  width = 6,
  units = "in",
  res = 600
)

# recreate this plot on the observation scale (used in main text)
p1_alt <- plot_parameters(
  mod_nestbuilt,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "A. nest built",
  transform = plogis
)

p2_alt <- plot_parameters(
  mod_clutchsize,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "B. clutch size",
  transform = exp
)

p3_alt <- plot_parameters(
  mod_eggsfertile,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "C. fertile eggs",
  transform = plogis
)

p4_alt <- plot_parameters(
  mod_hatchedeggs,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "D. hatched eggs",
  transform = plogis
)

p5_alt <- plot_parameters(
  mod_chicksfledged,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "E. chicks fledged",
  transform = plogis
)

p6_alt <- plot_parameters(
  mod_independentchicks,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "F. independent chicks",
  transform = plogis
)

p7_alt <- plot_parameters(
  mod_male_independentchicks,
  pars = c(
    "(Intercept)",
    "SEX_SPECIFIC_PAIR_TYPECfF1m", 
    "SEX_SPECIFIC_PAIR_TYPECfGm",
    "SEX_SPECIFIC_PAIR_TYPEF1F1",
    "SEX_SPECIFIC_PAIR_TYPEF1fCm",
    "SEX_SPECIFIC_PAIR_TYPEGfCm"
  ),
  title = "G. male independent chicks",
  transform = plogis
)

# combine these into a single plot
fig1_obs_scale <- (p1_alt | p2_alt ) /
  (p3_alt | p4_alt) /
  (p5_alt | p6_alt) / 
  (p7_alt | plot_spacer()) 

ggsave(
  file = "outputs/figures/Fig1-observation-scale.png",
  plot = fig1_obs_scale,
  device = png,
  height = 12,
  width = 6,
  units = "in",
  res = 600
)

# extract parameters for each model
all_pars <- c(
  "(Intercept)",
  "SEX_SPECIFIC_PAIR_TYPECfF1m", 
  "SEX_SPECIFIC_PAIR_TYPECfGm",
  "SEX_SPECIFIC_PAIR_TYPEF1F1",
  "SEX_SPECIFIC_PAIR_TYPEF1fCm",
  "SEX_SPECIFIC_PAIR_TYPEGfCm"
)
pars_nestbuilt <- extract_parameters(
  mod_nestbuilt,
  pars = all_pars,
  model = "nestbuilt",
  transform = plogis
)
pars_clutchsize <- extract_parameters(
  mod_clutchsize,
  pars = all_pars,
  model = "clutchsize",
  transform = exp
)
pars_eggsfertile <- extract_parameters(
  mod_eggsfertile,
  pars = all_pars,
  model = "eggsfertile",
  transform = plogis
)
pars_hatchedeggs <- extract_parameters(
  mod_hatchedeggs,
  pars = all_pars,
  model = "hatchedeggs",
  transform = plogis
)
pars_chicksfledged <- extract_parameters(
  mod_chicksfledged,
  pars = all_pars,
  model = "chicksfledged",
  transform = plogis
)
pars_independentchicks <- extract_parameters(
  mod_independentchicks,
  pars = all_pars,
  model = "independentchicks",
  transform = plogis
)

# combine these all together so we can multiply the different
#    stages together
pars_all <- bind_rows(
  pars_nestbuilt,
  pars_clutchsize,
  pars_eggsfertile,
  pars_hatchedeggs,
  pars_chicksfledged, 
  pars_independentchicks
)

# claculate expected number of independent chicks by cross type
pars_all <- pars_all %>%
  pivot_wider(
    id_cols = c(iter, group),
    names_from = model,
    values_from = value
  ) %>%
  mutate(
    expected_independent_chicks = 
      nestbuilt * clutchsize * eggsfertile * 
      hatchedeggs * chicksfledged * independentchicks,
    expected_independent_chicks_reduced_stages =
      clutchsize * hatchedeggs * chicksfledged * independentchicks,
  ) %>%
  pivot_longer(
    cols = c(
      nestbuilt, 
      clutchsize, 
      eggsfertile, 
      hatchedeggs, 
      chicksfledged, 
      independentchicks, 
      expected_independent_chicks,
      expected_independent_chicks_reduced_stages
    )
  )

# calculate a table of plotted values for each parameter in the model
raw_values <- pars_all %>%
  group_by(group, name) %>%
  summarise(
    q2.5 = quantile(value, probs = 0.025),
    q10 = quantile(value, probs = 0.1),
    q50 = quantile(value, probs = 0.5),
    q90 = quantile(value, probs = 0.9),
    q97.5 = quantile(value, probs = 0.975)
  ) %>%
  arrange(name, group)

# save to a file
write.csv(raw_values, file = "outputs/tables/Table_raw_values.csv")

# calculate probability that overall outcomes are greater than CC cross type
expected_chicks_by_crosstype <- compare_cross_type(
  pars_all,
  var = "expected_independent_chicks"
)

# repeat, for reduced stages (conditional on nest formation and fertile eggs)
expected_chicks_reduced_by_crosstype <- compare_cross_type(
  pars_all,
  var = "expected_independent_chicks_reduced_stages"
)

# plot the expected number of independent chicks by cross type
expected_independent_chicks <- plot_expected_independent_chicks(
  pars_all,
  var = "expected_independent_chicks",
  log_x = TRUE
)
expected_independent_chicks_reduced <- plot_expected_independent_chicks(
  pars_all,
  var = "expected_independent_chicks_reduced_stages",
  log_x = TRUE
)

# save these plots to files
ggsave(
  filename = "outputs/figures/Fig2-all-stages.png",
  plot = expected_independent_chicks,
  device = png,
  width = 5,
  height = 5,
  units = "in",
  res = 600
)
ggsave(
  filename = "outputs/figures/Fig2.png",
  plot = expected_independent_chicks_reduced,
  device = png,
  width = 5,
  height = 5,
  units = "in",
  res = 600
)

# list random effects included in each model
random_effect_list <- c(
  "DAM_ORIGIN_AGE",  # nest built
  "SIRE_ORIGIN_AGE", # nest built
  "Breeding_Season",  # nest built
  "Sire_ID",         # nest built
  "Dam_ID",          # nest built
  "Dam_Sire",
  "FOSTERED_EGG",   # hatched eggs only
  "FOSTERED_EGG_CHICK"   # fledged/independent chicks only
)

# extract random effects for each model
re_nestbuilt <- lapply(
  random_effect_list[1:5],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_nestbuilt,
  .y = "nestbuilt"
)
re_clutchsize <- lapply(
  random_effect_list[1:6],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_clutchsize,
  .y = "clutchsize"
)
re_eggsfertile <- lapply(
  random_effect_list[1:6],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_eggsfertile,
  .y = "eggsfertile"
)
re_hatchedeggs <- lapply(
  random_effect_list[1:7],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_hatchedeggs,
  .y = "hatchedeggs"
)
re_chicksfledged <- lapply(
  random_effect_list[c(1:6, 8)],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_chicksfledged,
  .y = "chicksfledged"
)
re_independent_chicks <- lapply(
  random_effect_list[c(1:6, 8)],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_independentchicks,
  .y = "independent_chicks"
)
re_male_independent_chicks <- lapply(
  random_effect_list[c(1:6, 8)],
  function(x, .x, .y) extract_random_effects(x = .x, term = x, model = .y),
  .x = mod_male_independentchicks,
  .y = "male_independent_chicks"
)

# and plot these for each model
random_plots_nestbuilt <- lapply(re_nestbuilt, plot_random_effects)
random_plots_clutchsize <- lapply(re_clutchsize, plot_random_effects)
random_plots_eggsfertile <- lapply(re_eggsfertile, plot_random_effects)
random_plots_hatchedeggs <- lapply(re_hatchedeggs, plot_random_effects)
random_plots_chicksfledged <- lapply(
  re_chicksfledged, 
  plot_random_effects
)
random_plots_independent_chicks <- lapply(
  re_independent_chicks,
  plot_random_effects
)
random_plots_male_independent_chicks <- lapply(
  re_male_independent_chicks,
  plot_random_effects
)

# plot nestbuilt terms in a single plot
re_plot_nestbuilt <- 
  (random_plots_nestbuilt[[1]] / random_plots_nestbuilt[[2]] / random_plots_nestbuilt[[3]]) | 
  random_plots_nestbuilt[[4]] |
  random_plots_nestbuilt[[5]]
ggsave(
  filename = "outputs/figures/re_nestbuilt.png",
  plot = re_plot_nestbuilt,
  device = png,
  width = 8,
  height = 8,
  units = "in",
  res = 600
)

# plot clutchsize terms in a single plot
re_plot_clutchsize <- 
  (random_plots_clutchsize[[1]] / random_plots_clutchsize[[2]] / random_plots_clutchsize[[3]]) |
  (random_plots_clutchsize[[4]] / random_plots_clutchsize[[5]]) |
  random_plots_clutchsize[[6]]
ggsave(
  filename = "outputs/figures/re_clutchsize.png",
  plot = re_plot_clutchsize,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# plot eggsfertile terms in a single plot
re_plot_eggsfertile <- 
  (random_plots_eggsfertile[[1]] / random_plots_eggsfertile[[2]] / random_plots_eggsfertile[[3]]) |
  (random_plots_eggsfertile[[4]] / random_plots_eggsfertile[[5]]) |
  random_plots_eggsfertile[[6]]
ggsave(
  filename = "outputs/figures/re_eggsfertile.png",
  plot = re_plot_eggsfertile,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# plot hatchedeggs terms in a single plot
re_plot_hatchedeggs <- 
  (random_plots_hatchedeggs[[1]] / random_plots_hatchedeggs[[2]] / 
     random_plots_hatchedeggs[[3]] / random_plots_hatchedeggs[[7]]) |
  (random_plots_hatchedeggs[[4]] / random_plots_hatchedeggs[[5]]) |
  random_plots_hatchedeggs[[6]]
ggsave(
  filename = "outputs/figures/re_hatchedeggs.png",
  plot = re_plot_hatchedeggs,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# plot chicksfledged terms in a single plot
re_plot_chicksfledged <- 
  (random_plots_chicksfledged[[1]] / random_plots_chicksfledged[[2]] / 
     random_plots_chicksfledged[[3]] / random_plots_chicksfledged[[7]]) |
  (random_plots_chicksfledged[[4]] / random_plots_chicksfledged[[5]]) |
  random_plots_chicksfledged[[6]]
ggsave(
  filename = "outputs/figures/re_chicksfledged.png",
  plot = re_plot_chicksfledged,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# plot independent_chicks terms in a single plot
re_plot_independent_chicks <- 
  (random_plots_independent_chicks[[1]] / random_plots_independent_chicks[[2]] / 
     random_plots_independent_chicks[[3]] / random_plots_independent_chicks[[7]]) |
  (random_plots_independent_chicks[[4]] / random_plots_independent_chicks[[5]]) |
  random_plots_independent_chicks[[6]]
ggsave(
  filename = "outputs/figures/re_independent_chicks.png",
  plot = re_plot_independent_chicks,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# plot male_independent_chicks terms in a single plot
re_plot_male_independent_chicks <- 
  (random_plots_male_independent_chicks[[1]] / random_plots_male_independent_chicks[[2]] / 
     random_plots_male_independent_chicks[[3]] / random_plots_male_independent_chicks[[7]]) |
  (random_plots_male_independent_chicks[[4]] / random_plots_male_independent_chicks[[5]]) |
  random_plots_male_independent_chicks[[6]]
ggsave(
  filename = "outputs/figures/re_male_independent_chicks.png",
  plot = re_plot_male_independent_chicks,
  device = png,
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

# calculate probability that each cross type is greater than the others
pr_nestbuilt <- probability_greater_than_others(mod_nestbuilt)
pr_clutchsize <- probability_greater_than_others(mod_clutchsize)
pr_eggsfertile <- probability_greater_than_others(mod_eggsfertile)
pr_hatchedeggs <- probability_greater_than_others(mod_hatchedeggs)
pr_chicksfledged <- probability_greater_than_others(mod_chicksfledged)
pr_independentchicks <- probability_greater_than_others(mod_independentchicks)
pr_male_independentchicks <- probability_greater_than_others(mod_male_independentchicks)

# extract comparisons against CC, combine all models, and save
pr_all <- bind_rows(
  pr_nestbuilt[1:6, ],
  pr_clutchsize[1:6, ],
  pr_eggsfertile[1:6, ],
  pr_hatchedeggs[1:6, ],
  pr_chicksfledged[1:6, ],
  pr_independentchicks[1:6, ],
  pr_male_independentchicks[1:6, ]
) %>%
  mutate(
    model = rep(
      c("nestbuilt", "clutchsize", "eggsfertile", 
        "hatchedeggs", "chicksfledged",
        "independentchicks", "male_independentchicks"),
      each = 6
    ),
    value = 1 - value,
    group = gsub(">", "<", group)
  )
write.csv(pr_all, file = "outputs/tables/Table3.csv")
