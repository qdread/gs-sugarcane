# Plots to show the relationships between traits in different crop cycles

source('GS_sugarcane_2022_loaddata.R')

# Reshape so we have crop cycles in separate columns
phenotypes_cc <- dcast(phenotypes_long, Clone + Rep + Row + Column + trait ~ crop_cycle)

# Paired plots for the three possible comparisons
library(ggplot2)

yx <- geom_abline(slope=1, intercept=0, color='red')
fit <- stat_smooth(method='lm', se=FALSE)

ggplot(phenotypes_cc, aes(x = `Plant Cane_2017`, y = `Ratoon 1_2018`)) +
  geom_point() + facet_wrap(~ trait, scales = 'free') + yx + fit

ggplot(phenotypes_cc, aes(x = `Plant Cane_2017`, y = `Ratoon 2_2019`)) +
  geom_point() + facet_wrap(~ trait, scales = 'free') + yx + fit

ggplot(phenotypes_cc, aes(x = `Ratoon 1_2018`, y = `Ratoon 2_2019`)) +
  geom_point() + facet_wrap(~ trait, scales = 'free') + yx + fit

# What are the correlations?

phenotypes_cc[, .(
    cor12 = cor(`Plant Cane_2017`, `Ratoon 1_2018`, use = 'pairwise.complete.obs'),
    cor13 = cor(`Plant Cane_2017`, `Ratoon 2_2019`, use = 'pairwise.complete.obs'),
    cor23 = cor(`Ratoon 1_2018`, `Ratoon 2_2019`, use = 'pairwise.complete.obs')
    ), by = .(trait)]
