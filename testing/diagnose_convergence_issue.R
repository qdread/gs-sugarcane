# diagnose convergence

mods <- list()
trts <- unique(phenotypes_long$trait)
trts <- trts[!trts %in% 'diam']
warn <- numeric(length(trts))

for (i in 1:length(trts)) {
  trt <- trts[i]
  message(trt)
  tryCatch(
    mods[[i]] <- lmer(value ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = phenotypes_long[trait == trt]),
    warning = function(w) {
      print(w)
      warn[i] <<- 1
    })
}

data.frame(trts, warn)
# Doesn't converge for stalk_ha

library(ggplot2)
ggplot(phenotypes[!crop_cycle %in% 'Ratoonability'], aes(x = stalk_ha, group = crop_cycle, fill = crop_cycle)) + geom_density(alpha = 0.5)

# If we select a different optimizer it works better
lmer(value ~ 0 + crop_cycle + (1|Clone) + (1|Row) + (1|Column), data = phenotypes_long[trait == 'stalk_ha'],
     control = lmerControl(optimizer = 'bobyqa'))

