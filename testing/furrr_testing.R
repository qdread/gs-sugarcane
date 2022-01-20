# Test parallel

library(furrr)


# Set up parallel processing if this code is running on ceres
if (grepl('ceres', system2('hostname', stdout = TRUE))) {
  options(mc.cores = 64)
  plan(
    list(tweak(multicore, workers = 16),
         tweak(multicore, workers = 16),
         tweak(multicore, workers = 16),
         tweak(multicore, workers = 16))
  )
}

fake_fn <- function(x) {
  t1 <- Sys.time()
  Sys.sleep(x)
  t2 <- Sys.time()
  hn <- system2('hostname', stdout = TRUE)
  nn <- Sys.info()[['nodename']]
  data.frame(start=t1, end=t2, hostname=hn, nodename=nn)
}

stuff <- data.frame(x = rep(5, 64))

output <- future_pmap_dfr(stuff, function(x) fake_fn(x))

# Other tries

# plan(
#   list(tweak(multisession, workers = 16),
#        tweak(multisession, workers = 16),
#        tweak(multisession, workers = 16),
#        tweak(multisession, workers = 16))
# )
# 
# plan(cluster(workers = availableWorkers()))