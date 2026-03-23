# Minimal day/night DVM example using one site and one scenario.
#
# Run this after installing/loading mizer:
# source(system.file("examples", "tpo2_dvm_minimal_example.R", package = "mizer"))

library(mizer)

csv <- system.file(
    "extdata",
    "profiles_emblematic_sites_hist_vs_ssp585_long.csv",
    package = "mizer"
)
profiles <- utils::read.csv(csv)
site <- unique(profiles$site)[1]
scenario <- "hist"

params <- newCommunityParams()
params_dvm <- enable_tpo2_dvm(
    params,
    profiles = profiles,
    site = site,
    scenario = scenario
)

out <- project_tpo2_dvm(params_dvm, t_max = 2, dt = 0.1, effort = 0)

cfg <- other_params(out$params_final)$tpo2_dvm
step1 <- out$profiles_saved[[1]]

# Visualise the day/night depth probabilities for the first species and size.
par(mfrow = c(2, 2))
matplot(
    cfg$geometry$z_mid,
    t(step1$p_day[1, 1:3, ]),
    type = "l",
    lty = 1,
    xlab = "Depth midpoint (m)",
    ylab = "Probability",
    main = paste("Day depth probabilities:", site, scenario)
)
matplot(
    cfg$geometry$z_mid,
    t(step1$p_night[1, 1:3, ]),
    type = "l",
    lty = 1,
    xlab = "Depth midpoint (m)",
    ylab = "Probability",
    main = paste("Night depth probabilities:", site, scenario)
)
matplot(
    cfg$geometry$z_mid,
    t(step1$realized_day$nu[1, 1:3, ]),
    type = "l",
    lty = 1,
    xlab = "Depth midpoint (m)",
    ylab = "nu",
    main = "Realised day nu profiles"
)
matplot(
    cfg$geometry$z_mid,
    t(step1$realized_night$g[1, 1:3, ]),
    type = "l",
    lty = 1,
    xlab = "Depth midpoint (m)",
    ylab = "g",
    main = "Realised night g profiles"
)

print(dim(step1$p_day))
print(dim(step1$realized_day$nu))
print(step1$effective$n_pp_eff[1:10])
