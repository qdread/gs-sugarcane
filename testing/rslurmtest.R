# Some code to test running code in parallel across multiple nodes and multiple cores

library(rslurm)

fake_fun <- function(iter,trait,crop_cycle,train_size) {
  Sys.sleep(5)
  data.frame(stuff = paste(iter,trait,crop_cycle,train_size))
}

sjob <- slurm_apply(fake_fun, combos, jobname = 'test', nodes = 4, cpus_per_node = 32, 
                    slurm_options = list(time = "1:00:00", mem = "10mb"))