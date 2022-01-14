# Fit mixed model to get BLUPs
# Separately for each crop cycle, and also one with crop cycle as a fixed effect.

library(lme4)



lmm_stkwt_2017 <- lmer(stkwt_kg ~ 1 + (1|Clone) + (1|Row) + (1|Column), data = phenotypes, subset = (crop_cycle == 'Plant Cane_2017'))

blup_stkwt_2017 <- ranef(lmm_stkwt_2017)[['Clone']] + fixef(lmm_stkwt_2017)
all.equal(blup_stkwt_2017, coef(lmm_stkwt_2017)[['Clone']])
blup_stkwt_2017 <- coef(lmm_stkwt_2017)[['Clone']]

lmm_stkwt_2018 <- lmer(stkwt_kg ~ 1 + (1|Clone) + (1|Row) + (1|Column), data = phenotypes, subset = (crop_cycle == 'Ratoon 1_2018'))
blup_stkwt_2018 <- coef(lmm_stkwt_2018)[['Clone']]

# Get the BLUPs for all put together, with crop cycle as fixed effect.
phenotypes_3crops <- phenotypes[!crop_cycle %in% 'Ratoonability']
lmm_stkwt <- lmer(stkwt_kg ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = phenotypes_3crops)

blup_stkwt <- outer(ranef(lmm_stkwt)[['Clone']][['(Intercept)']], fixef(lmm_stkwt), `+`)
row.names(blup_stkwt) <- row.names(ranef(lmm_stkwt)[['Clone']])


blup_trait <- function(trait, crop_cycle_to_use) {
  lmm <- lmer(as.formula(glue('{trait} ~ 1 + (1|Clone) + (1|Row) + (1|Column)')), data = phenotypes, subset = (crop_cycle == crop_cycle_to_use))
  blup <- coef(lmm)[['Clone']]
  data.frame(Clone = row.names(blup), BLUP = blup[,1])
}

blup_trait('stkwt_kg', 'Ratoon 1_2018')
