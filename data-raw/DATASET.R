Z_missing = read.csv("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/celia.escribe/git_repositories/BDcocolasso/data-raw/simulated_data.csv", header=TRUE)
Z_missing = Z_missing[,2:dim(Z_missing)[2]]
simulated_data_missing = cbind(data.frame(y = Z_missing[,1]),Z_missing[,2:dim(Z_missing)[2]])

usethis::use_data(simulated_data_missing, overwrite = TRUE)

Z_add = read.csv("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/celia.escribe/git_repositories/BDcocolasso/data-raw/simulated_data_additive.csv", header = TRUE)

Z_add = Z_add[,2:dim(Z_add)[2]]
simulated_data_additive = cbind(data.frame(y = Z_add[,1]),Z_add[,2:dim(Z_add)[2]])

usethis::use_data(simulated_data_additive, overwrite = TRUE)

Z_missing_block = read.csv("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/celia.escribe/git_repositories/BDcocolasso/data-raw/simulated_data_block_missing.csv", header = TRUE)

Z_missing_block = Z_missing_block[,2:dim(Z_missing_block)[2]]
simulated_data_missing_block = cbind(data.frame(y = Z_missing_block[,1]),Z_missing_block[,2:dim(Z_missing_block)[2]])

usethis::use_data(simulated_data_missing_block, overwrite = TRUE)

Z_additive_block = read.csv("/mnt/GREENWOOD_JBOD1/GREENWOOD_SCRATCH/celia.escribe/git_repositories/BDcocolasso/data-raw/simulated_data_block_additive.csv", header = TRUE)

Z_additive_block = Z_additive_block[,2:dim(Z_additive_block)[2]]
simulated_data_additive_block = cbind(data.frame(y = Z_additive_block[,1]),Z_additive_block[,2:dim(Z_additive_block)[2]])

usethis::use_data(simulated_data_additive_block, overwrite = TRUE)