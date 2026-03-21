
make_dvm_forcing <- function(params, z_mid,
                             T_day = c(10, 10, 10),
                             T_night = c(10, 10, 10),
                             pO2_day = c(20, 20, 20),
                             pO2_night = c(20, 20, 20),
                             L_day = c(1, 1, 1),
                             L_night = c(0, 0, 0),
                             n_pp_day = NULL,
                             n_pp_night = NULL) {
    if (is.null(n_pp_day)) {
        n_pp_day <- matrix(rep(initialNResource(params), each = length(z_mid)),
                           nrow = length(z_mid))
    }
    if (is.null(n_pp_night)) {
        n_pp_night <- n_pp_day
    }
    list(
        z_mid = z_mid,
        T_day = function(t, depth) T_day,
        T_night = function(t, depth) T_night,
        pO2_day = function(t, depth) pO2_day,
        pO2_night = function(t, depth) pO2_night,
        L_day = function(t, depth) L_day,
        L_night = function(t, depth) L_night,
        n_pp_day = function(t, depth) n_pp_day,
        n_pp_night = function(t, depth) n_pp_night
    )
}

test_that("day and night probabilities remain normalised over depth", {
    params <- newCommunityParams(no_w = 6)
    z_mid <- c(5, 25, 90)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid))
    out <- project_tpo2_dvm(params, t_max = 2 / 365, dt = 1 / 365, effort = 0)
    prof <- out$history
    expect_equal(apply(prof$p_day, c(1, 2, 3), sum),
                 array(1, dim = dim(prof$p_day)[1:3]), tolerance = 1e-8)
    expect_equal(apply(prof$p_night, c(1, 2, 3), sum),
                 array(1, dim = dim(prof$p_night)[1:3]), tolerance = 1e-8)
})

test_that("hard reachability cutoff is enforced", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(0, 100, 400)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid),
                              dvm_pars = list(dmax_const = 50,
                                              hard_reachability_cutoff = TRUE))
    cfg <- other_params(params)$tpo2_dvm
    nu <- array(rep(c(0, 0, 1), each = prod(dim(getMaxIntakeRate(params)))),
                dim = c(dim(getMaxIntakeRate(params)), length(z_mid)))
    mu <- nu * 0
    P <- .dvm_movement_kernel(params, nu, mu, "day")
    expect_equal(P[, , 1, 3], array(0, dim = dim(P[, , 1, 3])))
})

test_that("uniform depth profiles keep movement approximately uniform", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(10, 30, 60)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid),
                              dvm_pars = list(beta_move = 0.1,
                                              hard_reachability_cutoff = FALSE,
                                              dmax_const = 1e6))
    cfg <- other_params(params)$tpo2_dvm
    nu <- array(1, dim = dim(cfg$p_day), dimnames = dimnames(cfg$p_day))
    mu <- array(1, dim = dim(cfg$p_day), dimnames = dimnames(cfg$p_day))
    P <- .dvm_movement_kernel(params, nu, mu, "day")
    p_new <- .dvm_apply_kernel(cfg$p_day, P)
    expect_equal(p_new, cfg$p_day, tolerance = 1e-2)
})

test_that("light affects encounter but not mortality directly", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(5, 25, 80)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(
        params, z_mid, L_day = c(1, 0.1, 0.01), L_night = c(0, 0, 0)
    ))
    cfg <- other_params(params)$tpo2_dvm
    local_day <- .dvm_local_fish_array(initialN(params), cfg$p_day)
    day <- .dvm_local_realized_rates(params, 0, "day", local_day,
                                     .dvm_get_phase_resource(params, 0, "day"),
                                     cfg$g_prev_day, effort = 0)
    expect_true(day$encounter[1, 1, 1] > day$encounter[1, 1, 3])
    expect_equal(day$mu - day$pred_mort,
                 array(rep(params@mu_b, length(z_mid)), dim = dim(day$mu)) +
                     array(rep(0, length(day$mu)), dim = dim(day$mu)),
                 tolerance = 1e-8)
})

test_that("same-phase memory is kept separate between day and night", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(5, 25, 80)
    p_day0 <- array(0, dim = c(dim(getMaxIntakeRate(params)), length(z_mid)))
    p_night0 <- p_day0
    p_day0[, , 1] <- 1
    p_night0[, , 3] <- 1
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid),
                              init = list(p_day0 = p_day0, p_night0 = p_night0))
    cfg <- other_params(params)$tpo2_dvm
    expect_false(isTRUE(all.equal(cfg$p_day, cfg$p_night)))
})

test_that("realised local fields respond to reshuffling", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(5, 25, 80)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid,
        L_day = c(1, 0.5, 0.01), pO2_day = c(3, 10, 10)
    ))
    cfg <- other_params(params)$tpo2_dvm
    local_prev <- .dvm_local_fish_array(initialN(params), cfg$p_day)
    fore <- .dvm_local_forecast_rates(params, 0, "day", local_prev,
                                      .dvm_get_phase_resource(params, 0, "day"),
                                      cfg$g_prev_day, effort = 0)
    P <- .dvm_movement_kernel(params, fore$nu, fore$mu, "day")
    p_new <- .dvm_apply_kernel(cfg$p_day, P)
    real <- .dvm_local_realized_rates(params, 0, "day",
                                      .dvm_local_fish_array(initialN(params), p_new),
                                      .dvm_get_phase_resource(params, 0, "day"),
                                      .dvm_transport_g(cfg$g_prev_day, cfg$p_day, P),
                                      effort = 0)
    expect_false(isTRUE(all.equal(fore$nu, real$nu)))
})

test_that("effective rates have the expected dimensions", {
    params <- newCommunityParams(no_w = 4)
    z_mid <- c(5, 25, 80)
    params <- enable_tpo2_dvm(params, make_dvm_forcing(params, z_mid))
    out <- project_tpo2_dvm(params, t_max = 1 / 365, dt = 1 / 365, effort = 0)
    cfg <- other_params(out$params_final)$tpo2_dvm
    expect_identical(dim(cfg$effective$e_eff), dim(getMaxIntakeRate(params)))
    expect_identical(dim(cfg$effective$feeding_level_eff), dim(getMaxIntakeRate(params)))
    expect_identical(dim(cfg$effective$mort_eff), dim(getMaxIntakeRate(params)))
    expect_identical(dim(cfg$effective$pred_rate_eff), dim(getPredRate(params)))
    expect_identical(length(cfg$effective$n_pp_eff), length(initialNResource(params)))
})

test_that("depth weights integrate uneven grids correctly", {
    w <- .dvm_depth_weights(c(5, 20, 80))
    expect_equal(sum(w), 1)
    expect_equal(w, c(15, 37.5, 60) / 112.5)
})


test_that("depth weights use supplied top and bottom bounds", {
    w <- .dvm_depth_weights(c(5, 20, 80), z_top = c(0, 10, 40), z_bot = c(10, 40, 120))
    expect_equal(w, c(10, 30, 80) / 120)
})

test_that("profile_data can generate a depth-resolved resource spectrum", {
    params <- newCommunityParams(no_w = 4)
    profile_data <- data.frame(
        depth_mid_m = c(5, 25, 80),
        depth_top_m = c(0, 10, 40),
        depth_bot_m = c(10, 40, 120),
        T_mean = c(12, 10, 8),
        pO2_mean = c(8, 6, 4),
        zmeso_mean = c(0.1, 0.5, 1),
        zmicro_mean = c(1, 0.5, 0.1),
        I_day_rel = c(1, 0.4, 0.05),
        I_night_rel = c(0, 0, 0)
    )
    params <- enable_tpo2_dvm(params, list(profile_data = profile_data))
    day_resource <- .dvm_get_phase_resource(params, 0, "day")
    cfg <- other_params(params)$tpo2_dvm
    expect_equal(cfg$z_top, profile_data$depth_top_m)
    expect_equal(cfg$z_bot, profile_data$depth_bot_m)
    expect_identical(dim(day_resource), c(nrow(profile_data), length(initialNResource(params))))
    expect_true(sum(day_resource[1, ]) > 0)
    expect_false(isTRUE(all.equal(day_resource[1, ], day_resource[3, ])))
})


test_that("depth slicing preserves species-by-size matrices for single-species models", {
    x <- array(seq_len(12), dim = c(1, 4, 3),
               dimnames = list(species = "Community", size = 1:4, depth = 1:3))
    expect_identical(dim(.dvm_slice_species_size(x, 2)), c(1L, 4L))
    y <- array(seq_len(15), dim = c(1, 5, 3),
               dimnames = list(species = "Community", prey_size = 1:5, depth = 1:3))
    expect_identical(dim(.dvm_slice_pred_rate(y, 2)), c(1L, 5L))
})
