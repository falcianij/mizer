params <- newCommunityParams()
forcing <- list(
    times = c(0, 5),
    T_vec = c(10, 12),
    pO2_vec = c(20, 12)
)
params <- enable_tpo2_community(params, forcing = forcing)
out <- project_tpo2(params, t_max = 1, dt = 0.1, effort = 0)
diag0 <- getTPO2Diagnostics(params, t = 0)
print(diag0$f)
print(diag0$U)
print(diag0$g)
plot(getTimes(out$sim), getBiomass(out$sim), type = "l",
     xlab = "Time", ylab = "Biomass")
