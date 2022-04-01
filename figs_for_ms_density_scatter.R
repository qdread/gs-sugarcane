# Density plots with no y-axis labels and all on the same image
# Also coefficients of variation for all traits

library(tidyverse)
library(readxl)
library(ggrepel)
library(grid)
library(gridExtra)

# Lookup to convert some of the trait names to standardized abbreviations. Otherwise convert to uppercase.
trait_lookup <- c('diam' = 'SD', 'stalk_ha' = 'SP', 'Fiber' = 'FC', 'Sucrose' = 'SC', 'stkwt_kg' = 'SW')
model_lookup <- c('SVMlinear' = 'SVM', 'BayesA' = 'Bayes A', 'BayesB' = 'Bayes B', 'RF' = 'Random forest')
traits_exclude <- c('TRS', 'Purity') # Do not report these traits.
models_exclude <- c('SVMradial', 'SVMsigmoid') # Do not report these models.

metrics_gs <- read_csv('project/output/GS_all_metrics.csv') %>%
  filter(!model %in% models_exclude, !trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)),
         model = str_replace_all(model, pattern = model_lookup))

h2_all <- read_csv('project/output/h2_all_traits.csv') %>%
  filter(!trait %in% traits_exclude) %>%
  mutate(trait = toupper(str_replace_all(trait, pattern = trait_lookup)))

physical_traits <- c("SW", "SD", "BRIX", "FC", "POL", "SC", "SP")
economic_traits <- c("TCH", "CRS", "TSH", "EI")

theme_set(theme_bw() +
            theme(strip.background = element_blank(),
                  legend.position = 'bottom'))


pheno_file <- 'project/data/Phenotype_updated_2017-18_IL_analysis_120921.xlsx'
sheet_names <- excel_sheets(pheno_file)
phenotype_sheets <- map(sheet_names[1:3], ~ cbind(crop_cycle = ., read_xlsx(pheno_file, sheet = ., na = c('.', '-', '-!'))))

# Clean up so they can be concatenated
names(phenotype_sheets[[2]])[10] <- 'diam'
names(phenotype_sheets[[3]])[10] <- 'diam'
phenotype_sheets[[1]]$Crop <- as.character(phenotype_sheets[[1]]$Crop)
phenotype_sheets[[2]]$Crop <- as.character(phenotype_sheets[[2]]$Crop)

phenotypes <- bind_rows(phenotype_sheets)

# Convert trait names to new names
names(phenotypes)[match(names(trait_lookup), names(phenotypes))] <- trait_lookup
names(phenotypes)[9:21] <- toupper(names(phenotypes)[9:21])


# Density plots -----------------------------------------------------------

dens_plots <- map2(c(sort(physical_traits), sort(economic_traits)), LETTERS[1:11], ~
                    ggplot(phenotypes, aes_string(x = .x, fill = 'crop_cycle')) +
                    geom_density(alpha = 0.6) +
                    geom_text(aes(x = -Inf, y = Inf, label = .y, hjust = -1, vjust = 1.5), size = 6) +
                    scale_fill_viridis_d(name = 'Crop cycle', labels = c('plant cane', 'ratoon 1', 'ratoon 2')) +
                    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
                    theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                          legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

dens_plots[[9]] <- dens_plots[[9]] + theme(legend.position = c(0.7, 0.7), legend.background = element_blank())



png('project/figs/densityplots.png', height = 7, width = 10, res = 400, units = 'in')
do.call(grid.arrange, c(dens_plots, nrow = 3))
dev.off()


# SD and CV ---------------------------------------------------------------

phenotypes_long <- phenotypes %>%
  pivot_longer(cols = SW:EI, names_to = 'trait', values_to = 'value') %>%
  select(crop_cycle, trait, value) %>%
  filter(trait %in% c(physical_traits, economic_traits)) %>%
  mutate(crop_cycle = case_when(
    crop_cycle == 'Plant Cane_2017' ~ 'PlantCane',
    crop_cycle == 'Ratoon 1_2018' ~ 'Ratoon1',
    crop_cycle == 'Ratoon 2_2019' ~ 'Ratoon2'
  ))

CV_bycycle <- phenotypes_long %>%
  group_by(crop_cycle, trait) %>%
  summarize(trait_sd = sd(value, na.rm = TRUE),
            trait_mean = mean(value, na.rm = TRUE),
            CV = trait_sd / trait_mean)

CV_combined <- phenotypes_long %>%
  group_by(trait) %>%
  summarize(trait_sd = sd(value, na.rm = TRUE),
            trait_mean = mean(value, na.rm = TRUE),
            CV = trait_sd / trait_mean) %>%
  mutate(crop_cycle = 'combined')

CV_data <- bind_rows(CV_bycycle, CV_combined)

avg_r_bycycle <- metrics_gs %>%
  group_by(crop_cycle, trait) %>%
  summarize(r = mean(r))

avg_r_combined <- metrics_gs %>%
  group_by(trait) %>%
  summarize(r = mean(r)) %>%
  mutate(crop_cycle = 'combined')

avg_r <- bind_rows(avg_r_bycycle, avg_r_combined)

# Combine summary stats
h2_all <- h2_all %>%
  mutate(crop_cycle = case_when(
    crop_cycle == 'Plant Cane_2017' ~ 'PlantCane',
    crop_cycle == 'Ratoon 1_2018' ~ 'Ratoon1',
    crop_cycle == 'Ratoon 2_2019' ~ 'Ratoon2',
    TRUE ~ 'combined'
  ))

summary_stats <- CV_data %>% left_join(avg_r) %>% left_join(h2_all)

write_csv(summary_stats, 'project/output/summary_stats.csv')


# Scatterplots of summary stats -------------------------------------------

ph2r <- ggplot(summary_stats, aes(x = h2_narrow, y = r, label = trait)) +
  geom_point() +
  geom_text_repel(alpha = 0.6) +
  facet_wrap(~ crop_cycle) +
  labs(x = 'Narrow-sense heritability', y = 'Prediction accuracy')

pcvr <- ggplot(summary_stats, aes(x = CV, y = r, label = trait)) +
  geom_point() +
  geom_text_repel(alpha = 0.6) +
  facet_wrap(~ crop_cycle) +
  labs(x = 'Trait coefficient of variation', y = 'Prediction accuracy')

ggsave('project/figs/scatter_predaccuracy_vs_h2.png', ph2r, height = 6, width = 6, dpi = 400)
ggsave('project/figs/scatter_predaccuracy_vs_cv.png', pcvr, height = 6, width = 6, dpi = 400)
