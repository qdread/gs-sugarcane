mdtest <- pmap(.f = gs_md_fun, .l = combos[1:2])

mdjobshort <- slurm_apply(f = gs_md_fun, params = combos[3:6], jobname = 'testgsmd', nodes = 1, cpus_per_node = 4,
                     global_objects = c('geno_mat', 'pheno_blups', 'n_folds'),
                     slurm_options = list(partition = 'short', time = '8:00:00', mem = '80gb'))