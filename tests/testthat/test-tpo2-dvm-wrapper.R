make_dvm_profiles <- function() {
    utils::read.csv(test_path("..", "..", "data",
                              "profiles_emblematic_sites_hist_vs_ssp585_long.csv"))
}

test_that("initial day and night probabilities sum to one over depth", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- enable_tpo2_dvm(newCommunityParams(), profiles = profiles,
                              site = site, scenario = "hist")
    cfg <- other_params(params)$tpo2_dvm
    expect_equal(apply(cfg$p_day, c(1, 2), sum),
                 array(1, dim = dim(getMaxIntakeRate(params))))
    expect_equal(apply(cfg$p_night, c(1, 2), sum),
                 array(1, dim = dim(getMaxIntakeRate(params))))
})

test_that("integrated to local reconstruction integrates back with dz", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- enable_tpo2_dvm(newCommunityParams(), profiles = profiles,
                              site = site, scenario = "hist")
    cfg <- other_params(params)$tpo2_dvm
    N_local <- .dvm_integrated_to_local(initialN(params), cfg$p_day,
                                        cfg$geometry$H_total, cfg$geometry$dz)
    reconstructed <- apply(sweep(N_local, 3, cfg$geometry$dz, "*"), c(1, 2), sum)
    expect_equal(reconstructed, initialN(params) * cfg$geometry$H_total)
})

test_that("resource spectra are non-negative and integrate correctly", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- newCommunityParams()
    dvm_cfg <- .tpo2_dvm_defaults()
    prof <- profiles[profiles$site == site & profiles$scenario == "hist", ]
    resource <- .dvm_build_resource_spectra_from_cobalt(params, prof, dvm_cfg)
    geom <- .dvm_depth_geometry(prof)
    expect_true(all(resource$n_pp_local >= 0))
    expect_equal(rowSums(resource$b_pp_local * rep(dw_full(params), each = nrow(resource$b_pp_local))),
                 resource$B_micro + resource$B_meso, tolerance = 1e-8)
    eff <- .dvm_effective_npp(resource$n_pp_local, geom$w_depth)
    manual <- as.numeric(crossprod(geom$w_depth, resource$n_pp_local))
    expect_equal(eff$n_pp_eff, manual)
})

test_that("resources are not redistributed with fish occupancy", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- enable_tpo2_dvm(newCommunityParams(), profiles = profiles,
                              site = site, scenario = "hist")
    cfg <- other_params(params)$tpo2_dvm
    shifted <- cfg$p_day
    shifted[] <- 0
    shifted[, , 1] <- 1
    expect_equal(cfg$resource$n_pp_local, cfg$resource$n_pp_local)
    expect_equal(other_params(params)$tpo2_dvm$resource$n_pp_local,
                 cfg$resource$n_pp_local)
})

test_that("uniform abiotic and biotic fields give approximately uniform movement", {
    params <- newCommunityParams()
    geom <- list(z_mid = c(10, 20, 30))
    nu <- array(1, dim = c(nrow(initialN(params)), ncol(initialN(params)), 3))
    mu <- array(1, dim = c(nrow(initialN(params)), ncol(initialN(params)), 3))
    dimnames(nu) <- c(dimnames(initialN(params)), list(depth = c("1", "2", "3")))
    dimnames(mu) <- dimnames(nu)
    dvm_cfg <- .tpo2_dvm_defaults()
    dvm_cfg$dmax_const <- 1e6
    dvm_cfg$eta_move <- 1
    P <- .dvm_movement_kernel(params, nu, mu, geom, dvm_cfg)
    expect_equal(P[1, 1, 1, ], rep(1 / 3, 3), tolerance = 1e-6)
})

test_that("movement kernel is row-stochastic over destinations", {
    params <- newCommunityParams()
    geom <- list(z_mid = c(5, 15, 25))
    nu <- array(runif(nrow(initialN(params)) * ncol(initialN(params)) * 3),
                dim = c(nrow(initialN(params)), ncol(initialN(params)), 3))
    mu <- array(runif(nrow(initialN(params)) * ncol(initialN(params)) * 3),
                dim = c(nrow(initialN(params)), ncol(initialN(params)), 3))
    dimnames(nu) <- c(dimnames(initialN(params)), list(depth = c("1", "2", "3")))
    dimnames(mu) <- dimnames(nu)
    P <- .dvm_movement_kernel(params, nu, mu, geom, .tpo2_dvm_defaults())
    expect_equal(apply(P, c(1, 2, 3), sum), array(1, dim = dim(P)[1:3]), tolerance = 1e-8)
})

test_that("realised local fields can differ from forecast fields after reshuffling", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- enable_tpo2_dvm(newCommunityParams(), profiles = profiles,
                              site = site, scenario = "hist")
    cfg <- other_params(params)$tpo2_dvm
    N_int <- initialN(params)
    abiotic_day <- list(
        temp = cfg$profiles$temp_C[match(cfg$geometry$depth$depth_idx, cfg$profiles$depth_idx)],
        pO2 = cfg$profiles$pO2_kPa[match(cfg$geometry$depth$depth_idx, cfg$profiles$depth_idx)],
        light = cfg$profiles$I_day_rel[match(cfg$geometry$depth$depth_idx, cfg$profiles$depth_idx)]
    )
    fish_prev <- .dvm_integrated_to_local(N_int, cfg$prev_day$p,
                                          cfg$geometry$H_total, cfg$geometry$dz)
    forecast <- .dvm_local_fields(params, fish_prev, cfg$resource$n_pp_local,
                                  abiotic_day, cfg$prev_day$g, effort = 0, t = 0)
    P <- .dvm_movement_kernel(params, forecast$nu, forecast$mu, cfg$geometry, cfg)
    p_new <- .dvm_apply_kernel(cfg$p_night, P)
    fish_real <- .dvm_integrated_to_local(N_int, p_new, cfg$geometry$H_total, cfg$geometry$dz)
    realised <- .dvm_local_fields(params, fish_real, cfg$resource$n_pp_local,
                                  abiotic_day, .dvm_transport_g(cfg$p_night, P, cfg$g_state_night),
                                  effort = 0, t = 0)
    expect_false(isTRUE(all.equal(forecast$nu, realised$nu, tolerance = 1e-12)))
})


test_that("enable_tpo2_dvm uses robust scalar forcing and expected setup dimensions", {
    profiles <- make_dvm_profiles()
    site <- unique(profiles$site)[1]
    params <- enable_tpo2_dvm(newCommunityParams(), profiles = profiles,
                              site = site, scenario = "hist")
    cfg <- other_params(params)$tpo2_dvm
    tpo2_cfg <- other_params(params)$tpo2
    expect_true(is.function(tpo2_cfg$forcing$T$fun))
    expect_true(is.function(tpo2_cfg$forcing$pO2$fun))
    expect_identical(dim(cfg$resource$n_pp_local),
                     c(length(cfg$geometry$dz), length(w_full(params))))
    expect_identical(dim(cfg$p_day),
                     c(dim(getMaxIntakeRate(params)), length(cfg$geometry$dz)))
})
