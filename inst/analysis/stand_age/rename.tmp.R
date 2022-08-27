load("midnorth.stand_age.simulations.RData.bkp")
names(boot_at_mle$loglik_surfaces)[2:4] <- paste0("true_", names(boot_at_mle$loglik_surfaces)[2:4])
names(boot_at_mle$loglik_surfaces)[5:7] <- paste0("sim_", names(boot_at_mle$loglik_surfaces)[5:7])
names(boot_at_null$loglik_surfaces)[5:7] <- paste0("sim_", names(boot_at_null$loglik_surfaces)[5:7])
names(boot_at_null$loglik_surfaces)[2:4] <- paste0("true_", names(boot_at_null$loglik_surfaces)[2:4])
save(boot_at_mle, boot_at_null, file="midnorth.stand_age.simulations.RData")

