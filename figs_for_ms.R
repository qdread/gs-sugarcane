# GS Sugarcane: figures formatted for manuscript
# QDR / 28 Mar 2022

library(tidyverse)

# Lookup to convert some of the trait names to standardized abbreviations. Otherwise convert to uppercase.
trait_lookup <- c('diam' = 'SD', 'stalk_ha' = 'SP', 'Fiber' = 'FC', 'Sucrose' = 'SC', 'stkwt_kg' = 'SW')
model_lookup <- c('SVMlinear' = 'SVM', 'BayesA' = 'Bayes A', 'BayesB' = 'Bayes B', 'RF' = 'Random forest')
traits_exclude <- c('TRS', 'Purity') # Do not report these traits.
models_exclude <- c('SVMradial', 'SVMsigmoid') # Do not report these models.



metrics_gs <- read_csv('project/output/GS_all_metrics.csv') %>%
  filter(!model %in% models_exclude, !trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)),
         model = str_replace_all(model, pattern = model_lookup))
metrics_tags <- read_csv('project/output/TAGS_all_metrics.csv') %>%
  filter(!trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)))

metrics_ts <- read_csv('project/output/TS_all_metrics.csv') %>%
  filter(!model %in% models_exclude, !trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)),
         model = str_replace_all(model, pattern = model_lookup))
metrics_md <- read_csv('project/output/MD_all_metrics.csv') %>%
  filter(!model %in% models_exclude, !trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)),
         model = str_replace_all(model, pattern = model_lookup))

metrics_ts <- bind_rows(
  metrics_ts,
  metrics_gs %>% mutate(training_size = 1)
)
metrics_md <- bind_rows(
  metrics_md,
  metrics_gs %>% mutate(marker_density = 1)
)

physical_traits <- c("SW", "SD", "BRIX", "FC", "POL", "SC", "SP")
economic_traits <- c("TCH", "CRS", "TSH", "EI")

theme_set(theme_bw() +
            theme(strip.background = element_blank(),
                  legend.position = 'bottom',
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text = element_text(color = 'black')))

colorblind_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fill_scale <- scale_fill_manual(values = colorblind_colors)
ci_null_line <- geom_hline(yintercept = 0.167, size = 1.2, linetype = 'dashed')
r_null_line <- geom_hline(yintercept = 0, size = 1.2, linetype = 'dashed')
color_scale <- scale_color_manual(values = colorblind_colors)
x_scale_ts <- scale_x_continuous(name = 'training set size', breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1), labels = scales::percent)
x_scale_md <- scale_x_continuous(name = 'marker density', breaks = c(0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1), labels = scales::percent)


# Scale intercepts --------------------------------------------------

# To plot them all on the same plot scale the intercepts within each trait (divide by SD)
metrics_gs <- metrics_gs %>%
  group_by(trait) %>%
  mutate(intercept = intercept/sd(intercept)) %>%
  ungroup
metrics_tags <- metrics_tags %>%
  group_by(trait) %>%
  mutate(intercept = intercept/sd(intercept)) %>%
  ungroup
metrics_ts <- metrics_ts %>%
  group_by(trait) %>%
  mutate(intercept = intercept/sd(intercept)) %>%
  ungroup
metrics_md <- metrics_md %>%
  group_by(trait) %>%
  mutate(intercept = intercept/sd(intercept)) %>%
  ungroup

# Global y scales ---------------------------------------------------------

ci_range <- range(c(metrics_gs$CI, metrics_tags$CI))
r_range <- range(c(metrics_gs$r, metrics_tags$r))
intercept_range <- range(c(metrics_gs$intercept, metrics_tags$intercept))
slope_range <- range(c(metrics_gs$slope, metrics_tags$slope))

ci_y_scale <- scale_y_continuous(name = 'coincidence index', limits = ci_range)
r_y_scale <- scale_y_continuous(name = 'prediction accuracy', limits = r_range)
intercept_y_scale <- scale_y_continuous(name = 'intercept', limits = intercept_range)
slope_y_scale <- scale_y_continuous(name = 'slope', limits = slope_range)

# Physical traits CI ------------------------------------------------------

p_ci_phys_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% physical_traits), 
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: physical traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + ci_null_line + ci_y_scale

p_ci_phys_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% physical_traits),
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: physical traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + ci_null_line + ci_y_scale

p_ci_phys_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% physical_traits), 
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: physical traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + ci_null_line + ci_y_scale


# Economic traits CI ------------------------------------------------------

p_ci_econ_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% economic_traits), 
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: economic traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + ci_null_line + ci_y_scale

p_ci_econ_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% economic_traits),
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: economic traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + ci_null_line + ci_y_scale

p_ci_econ_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% economic_traits), 
       aes(x = trait, y = CI, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Coincidence index: economic traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + ci_null_line + ci_y_scale


# Physical traits r -------------------------------------------------------

p_r_phys_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% physical_traits), 
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: physical traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + r_null_line + r_y_scale

p_r_phys_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% physical_traits),
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: physical traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + r_null_line + r_y_scale

p_r_phys_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% physical_traits), 
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: physical traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + r_null_line + r_y_scale


# Economic traits r -------------------------------------------------------

p_r_econ_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% economic_traits), 
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: economic traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + r_null_line + r_y_scale

p_r_econ_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% economic_traits),
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: economic traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + r_null_line + r_y_scale

p_r_econ_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% economic_traits), 
       aes(x = trait, y = r, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted correlation: economic traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + r_null_line + r_y_scale


# Physical traits slope ---------------------------------------------------

p_slope_phys_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% physical_traits), 
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: physical traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + slope_y_scale

p_slope_phys_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% physical_traits),
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: physical traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + slope_y_scale

p_slope_phys_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% physical_traits), 
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: physical traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + slope_y_scale


# Economic traits slope ---------------------------------------------------

p_slope_econ_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% economic_traits), 
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: economic traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + slope_y_scale

p_slope_econ_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% economic_traits),
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: economic traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + slope_y_scale

p_slope_econ_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% economic_traits), 
       aes(x = trait, y = slope, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Observed-predicted slope: economic traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + slope_y_scale

# Physical traits intercept ---------------------------------------------------

p_int_phys_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% physical_traits), 
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: physical traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + intercept_y_scale

p_int_phys_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% physical_traits),
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: physical traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + intercept_y_scale

p_int_phys_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% physical_traits), 
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: physical traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + intercept_y_scale


# Economic traits intercept ---------------------------------------------------

p_int_econ_pc <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'PlantCane', trait %in% economic_traits), 
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: economic traits', 'Crop cycle: plant cane (2017)') +
  fill_scale + intercept_y_scale

p_int_econ_1r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon1', trait %in% economic_traits),
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: economic traits', 'Crop cycle: first ratoon (2018)') +
  fill_scale + intercept_y_scale

p_int_econ_2r <- ggplot(metrics_gs %>% filter(crop_cycle %in% 'Ratoon2', trait %in% economic_traits), 
       aes(x = trait, y = intercept, fill = model, group = interaction(trait, model))) +
  geom_boxplot(position = 'dodge') +
  ggtitle('Scaled observed-predicted intercept: economic traits', 'Crop cycle: second ratoon (2019)') +
  fill_scale + intercept_y_scale


# Trait-assisted gs -------------------------------------------------------

fill_scale_tags <- scale_fill_brewer(palette = 'Dark2', name = 'predicted crop cycle', labels = c('first ratoon', 'second ratoon'))

p_ci_phys_tags <- ggplot(metrics_tags %>% filter(trait %in% physical_traits), aes(x = trait, y = CI, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ci_null_line +
  ci_y_scale +
  ggtitle('Coincidence index: physical traits', 'Trait-assisted GS')

p_ci_econ_tags <- ggplot(metrics_tags %>% filter(trait %in% economic_traits), aes(x = trait, y = CI, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ci_null_line +
  ci_y_scale +
  ggtitle('Coincidence index: economic traits', 'Trait-assisted GS')

p_r_phys_tags <- ggplot(metrics_tags %>% filter(trait %in% physical_traits), aes(x = trait, y = r, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  r_null_line +
  r_y_scale +
  ggtitle('Observed-predicted correlation: physical traits', 'Trait-assisted GS')

p_r_econ_tags <- ggplot(metrics_tags %>% filter(trait %in% economic_traits), aes(x = trait, y = r, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  r_null_line +
  r_y_scale +
  ggtitle('Observed-predicted correlation: economic traits', 'Trait-assisted GS')

p_slope_phys_tags <- ggplot(metrics_tags %>% filter(trait %in% physical_traits), aes(x = trait, y = slope, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ggtitle('Observed-predicted slope: physical traits', 'Trait-assisted GS')

p_slope_econ_tags <- ggplot(metrics_tags %>% filter(trait %in% economic_traits), aes(x = trait, y = slope, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ggtitle('Observed-predicted slope: economic traits', 'Trait-assisted GS')

p_int_phys_tags <- ggplot(metrics_tags %>% filter(trait %in% physical_traits), aes(x = trait, y = intercept, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ggtitle('Scaled observed-predicted intercept: physical traits', 'Trait-assisted GS')

p_int_econ_tags <- ggplot(metrics_tags %>% filter(trait %in% economic_traits), aes(x = trait, y = intercept, fill = crop_cycle, group = interaction(trait, crop_cycle))) +
  geom_boxplot() +
  fill_scale_tags +
  ggtitle('Scaled observed-predicted intercept: economic traits', 'Trait-assisted GS')


# Training size trends ----------------------------------------------------

pd <- position_dodge(width = 0.07)
bplw <- 0.2 # Boxplot line width should be thin.
r_y_extended <- scale_y_continuous(name = 'prediction accuracy', limits = c(-0.182, 0.505))
ci_y_extended <- scale_y_continuous(name = 'coincidence index', limits = c(0.06, 0.406))

### pc
p_r_phys_ts_pc <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'PlantCane'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing training set proportion', 'physical traits: plant cane (2017)')

p_r_econ_ts_pc <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'PlantCane'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing training set proportion', 'economic traits: plant cane (2017)')

p_ci_phys_ts_pc <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'PlantCane'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing training set proportion', 'physical traits: plant cane (2017)')

p_ci_econ_ts_pc <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'PlantCane'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing training set proportion', 'economic traits: plant cane (2017)')

### 1r

p_r_phys_ts_1r <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon1'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing training set proportion', 'physical traits: first ratoon (2018)')

p_r_econ_ts_1r <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon1'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing training set proportion', 'economic traits: first ratoon (2018)')

p_ci_phys_ts_1r <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon1'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing training set proportion', 'physical traits: first ratoon (2018)')

p_ci_econ_ts_1r <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon1'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing training set proportion', 'economic traits: first ratoon (2018)')

### 2r

p_r_phys_ts_2r <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon2'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing training set proportion', 'physical traits: second ratoon (2019)')

p_r_econ_ts_2r <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon2'), aes(x = training_size, group = interaction(training_size, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing training set proportion', 'economic traits: second ratoon (2019)')

p_ci_phys_ts_2r <- ggplot(metrics_ts %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon2'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing training set proportion', 'physical traits: second ratoon (2019)')

p_ci_econ_ts_2r <- ggplot(metrics_ts %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon2'), aes(x = training_size, group = interaction(training_size, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_ts + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing training set proportion', 'economic traits: second ratoon (2019)')


# Marker density trends ---------------------------------------------------

p_r_phys_md_pc <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'PlantCane'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing marker density', 'physical traits: plant cane (2017)')

p_r_econ_md_pc <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'PlantCane'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing marker density', 'economic traits: plant cane (2017)')

p_ci_phys_md_pc <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'PlantCane'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing marker density', 'physical traits: plant cane (2017)')

p_ci_econ_md_pc <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'PlantCane'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing marker density', 'economic traits: plant cane (2017)')

### 1r

p_r_phys_md_1r <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon1'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing marker density', 'physical traits: first ratoon (2018)')

p_r_econ_md_1r <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon1'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing marker density', 'economic traits: first ratoon (2018)')

p_ci_phys_md_1r <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon1'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing marker density', 'physical traits: first ratoon (2018)')

p_ci_econ_md_1r <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon1'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing marker density', 'economic traits: first ratoon (2018)')

### 2r

p_r_phys_md_2r <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon2'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Prediction accuracy with increasing marker density', 'physical traits: second ratoon (2019)')

p_r_econ_md_2r <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon2'), aes(x = marker_density, group = interaction(marker_density, model), y = r, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  r_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + r_y_extended +
  ggtitle('Prediction accuracy with increasing marker density', 'economic traits: second ratoon (2019)')

p_ci_phys_md_2r <- ggplot(metrics_md %>% filter(trait %in% physical_traits, crop_cycle %in% 'Ratoon2'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  theme(legend.position = c(0.7, 0.15)) +
  guides(fill = guide_legend(nrow = 3)) +
  ggtitle('Coincidence index with increasing marker density', 'physical traits: second ratoon (2019)')

p_ci_econ_md_2r <- ggplot(metrics_md %>% filter(trait %in% economic_traits, crop_cycle %in% 'Ratoon2'), aes(x = marker_density, group = interaction(marker_density, model), y = CI, fill = model)) +
  stat_summary(fun = median, geom = 'line', aes(group = model, color = model), alpha = 0.6, position = pd) +
  geom_boxplot(position = pd, lwd = bplw) +
  ci_null_line +
  facet_wrap(~ trait) +
  x_scale_md + fill_scale + color_scale + ci_y_extended +
  ggtitle('Coincidence index with increasing marker density', 'economic traits: second ratoon (2019)')


# Save everything ---------------------------------------------------------

h1 <- 4
w1 <- 5
h2 <- 4
w2 <- 4.5

ggsave('project/figs/mainanalysis_coincidenceindex_physical_plantcane.png', p_ci_phys_pc, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_coincidenceindex_physical_1stratoon.png', p_ci_phys_1r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_coincidenceindex_physical_2ndratoon.png', p_ci_phys_2r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_coincidenceindex_economic_plantcane.png', p_ci_econ_pc, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_coincidenceindex_economic_1stratoon.png', p_ci_econ_1r, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_coincidenceindex_economic_2ndratoon.png', p_ci_econ_2r, height = h2, width = w2, dpi = 400)

ggsave('project/figs/mainanalysis_predaccuracy_physical_plantcane.png', p_r_phys_pc, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_predaccuracy_physical_1stratoon.png', p_r_phys_1r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_predaccuracy_physical_2ndratoon.png', p_r_phys_2r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_predaccuracy_economic_plantcane.png', p_r_econ_pc, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_predaccuracy_economic_1stratoon.png', p_r_econ_1r, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_predaccuracy_economic_2ndratoon.png', p_r_econ_2r, height = h2, width = w2, dpi = 400)

ggsave('project/figs/mainanalysis_intercept_physical_plantcane.png', p_int_phys_pc, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_intercept_physical_1stratoon.png', p_int_phys_1r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_intercept_physical_2ndratoon.png', p_int_phys_2r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_intercept_economic_plantcane.png', p_int_econ_pc, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_intercept_economic_1stratoon.png', p_int_econ_1r, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_intercept_economic_2ndratoon.png', p_int_econ_2r, height = h2, width = w2, dpi = 400)

ggsave('project/figs/mainanalysis_slope_physical_plantcane.png', p_slope_phys_pc, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_slope_physical_1stratoon.png', p_slope_phys_1r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_slope_physical_2ndratoon.png', p_slope_phys_2r, height = h1, width = w1, dpi = 400)
ggsave('project/figs/mainanalysis_slope_economic_plantcane.png', p_slope_econ_pc, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_slope_economic_1stratoon.png', p_slope_econ_1r, height = h2, width = w2, dpi = 400)
ggsave('project/figs/mainanalysis_slope_economic_2ndratoon.png', p_slope_econ_2r, height = h2, width = w2, dpi = 400)

ggsave('project/figs/traitassisted_coincidenceindex_physical.png', p_ci_phys_tags, height = h1, width = w1, dpi = 400)
ggsave('project/figs/traitassisted_coincidenceindex_economic.png', p_ci_econ_tags, height = h2, width = w2, dpi = 400)
ggsave('project/figs/traitassisted_predaccuracy_physical.png', p_r_phys_tags, height = h1, width = w1, dpi = 400)
ggsave('project/figs/traitassisted_predaccuracy_economic.png', p_r_econ_tags, height = h2, width = w2, dpi = 400)
ggsave('project/figs/traitassisted_intercept_physical.png', p_int_phys_tags, height = h1, width = w1, dpi = 400)
ggsave('project/figs/traitassisted_intercept_economic.png', p_int_econ_tags, height = h2, width = w2, dpi = 400)
ggsave('project/figs/traitassisted_slope_physical.png', p_slope_phys_tags, height = h1, width = w1, dpi = 400)
ggsave('project/figs/traitassisted_slope_economic.png', p_slope_econ_tags, height = h2, width = w2, dpi = 400)

h3 <- 8
w3 <- 10
h4 <- 6
w4 <- 8

ggsave('project/figs/markerdensity_coincidenceindex_physical_plantcane.png', p_ci_phys_md_pc, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_coincidenceindex_physical_1stratoon.png', p_ci_phys_md_1r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_coincidenceindex_physical_2ndratoon.png', p_ci_phys_md_2r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_coincidenceindex_economic_plantcane.png', p_ci_econ_md_pc, height = h4, width = w4, dpi = 400)
ggsave('project/figs/markerdensity_coincidenceindex_economic_1stratoon.png', p_ci_econ_md_1r, height = h4, width = w4, dpi = 400)
ggsave('project/figs/markerdensity_coincidenceindex_economic_2ndratoon.png', p_ci_econ_md_2r, height = h4, width = w4, dpi = 400)

ggsave('project/figs/markerdensity_predaccuracy_physical_plantcane.png', p_r_phys_md_pc, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_predaccuracy_physical_1stratoon.png', p_r_phys_md_1r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_predaccuracy_physical_2ndratoon.png', p_r_phys_md_2r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/markerdensity_predaccuracy_economic_plantcane.png', p_r_econ_md_pc, height = h4, width = w4, dpi = 400)
ggsave('project/figs/markerdensity_predaccuracy_economic_1stratoon.png', p_r_econ_md_1r, height = h4, width = w4, dpi = 400)
ggsave('project/figs/markerdensity_predaccuracy_economic_2ndratoon.png', p_r_econ_md_2r, height = h4, width = w4, dpi = 400)

ggsave('project/figs/trainingsize_coincidenceindex_physical_plantcane.png', p_ci_phys_ts_pc, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_coincidenceindex_physical_1stratoon.png', p_ci_phys_ts_1r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_coincidenceindex_physical_2ndratoon.png', p_ci_phys_ts_2r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_coincidenceindex_economic_plantcane.png', p_ci_econ_ts_pc, height = h4, width = w4, dpi = 400)
ggsave('project/figs/trainingsize_coincidenceindex_economic_1stratoon.png', p_ci_econ_ts_1r, height = h4, width = w4, dpi = 400)
ggsave('project/figs/trainingsize_coincidenceindex_economic_2ndratoon.png', p_ci_econ_ts_2r, height = h4, width = w4, dpi = 400)

ggsave('project/figs/trainingsize_predaccuracy_physical_plantcane.png', p_r_phys_ts_pc, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_predaccuracy_physical_1stratoon.png', p_r_phys_ts_1r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_predaccuracy_physical_2ndratoon.png', p_r_phys_ts_2r, height = h3, width = w3, dpi = 400)
ggsave('project/figs/trainingsize_predaccuracy_economic_plantcane.png', p_r_econ_ts_pc, height = h4, width = w4, dpi = 400)
ggsave('project/figs/trainingsize_predaccuracy_economic_1stratoon.png', p_r_econ_ts_1r, height = h4, width = w4, dpi = 400)
ggsave('project/figs/trainingsize_predaccuracy_economic_2ndratoon.png', p_r_econ_ts_2r, height = h4, width = w4, dpi = 400)