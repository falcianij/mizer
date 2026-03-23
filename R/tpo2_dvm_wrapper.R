# Day/night DVM wrapper built on top of the temperature-oxygen extension.

.tpo2_dvm_defaults <- function() {
    list(
        enabled = TRUE,
        phase_weight_day = 0.5,
        phase_weight_night = 0.5,
        beta_move = 2,
        eta_move = 0.05,
        eps_norm = 1e-8,
        light_half_sat = 0.5,
        light_hill = 1,
        dmax_const = 150,
        dmax_beta = 0,
        hard_reachability_cutoff = TRUE,
        resource_input_units = "molN_m3",
        input_to_gN = 14,
        N_to_C = 5.625,
        C_to_wet = 10,
        w_group = c(micro = 1e-6, meso = 1e-3),
        sigma_group = c(micro = 0.7, meso = 1.0)
    )
}

.tpo2_dvm_get_config <- function(params) {
    cfg <- other_params(params)$tpo2_dvm
    if (is.null(cfg) || !isTRUE(cfg$enabled)) {
        stop("The tpo2 DVM extension is not enabled for this params object.")
    }
    cfg
}


.dvm_validate_setup <- function(params, geometry, resource, p_day0, p_night0,
                                g_state_day0, g_state_night0) {
    target_sp_size <- dim(getMaxIntakeRate(params))
    target_sp_size_depth <- c(target_sp_size, length(geometry$dz))
    if (!identical(dim(resource$n_pp_local), c(length(geometry$dz), length(w_full(params))))) {
        stop("n_pp_local must have dimensions depth x resource_size = ",
             paste(c(length(geometry$dz), length(w_full(params))), collapse = " x "),
             "; got ", paste(dim(resource$n_pp_local), collapse = " x "), ".")
    }
    arrays_to_check <- list(
        p_day0 = p_day0,
        p_night0 = p_night0,
        g_state_day0 = g_state_day0,
        g_state_night0 = g_state_night0
    )
    for (nm in names(arrays_to_check)) {
        x <- arrays_to_check[[nm]]
        if (!identical(dim(x), target_sp_size_depth)) {
            got_dim <- dim(x)
            stop(nm, " must have dimensions species x size x depth = ",
                 paste(target_sp_size_depth, collapse = " x "),
                 "; got ",
                 if (is.null(got_dim)) "<no dim>" else paste(got_dim, collapse = " x "),
                 ".")
        }
    }
    invisible(TRUE)
}

.dvm_normalise_probabilities <- function(p, name = "probability array") {
    dims <- dim(p)
    flat <- matrix(p, nrow = prod(dims[1:2]), ncol = dims[3])
    rs <- rowSums(flat)
    if (any(abs(rs - 1) > 1e-8)) {
        stop(name, " must sum to 1 over depth for each species and size.")
    }
    p
}


.dvm_prepare_profiles <- function(profile_df) {
    required_cols <- c(
        "depth_idx", "depth_mid_m", "depth_top_m", "depth_bot_m", "dz_m",
        "temp_C", "pO2_kPa", "zmicro", "zmeso", "I_day_rel", "I_night_rel"
    )
    missing_cols <- setdiff(required_cols, names(profile_df))
    if (length(missing_cols) > 0) {
        stop("Profile data frame is missing columns: ",
             paste(missing_cols, collapse = ", "))
    }
    keep <- stats::complete.cases(profile_df[, required_cols])
    trimmed <- profile_df[keep, , drop = FALSE]
    if (nrow(trimmed) == 0) {
        stop("No valid depth cells remain after removing rows with missing geometry, abiotic, light, and resource values.")
    }
    trimmed[order(trimmed$depth_idx), , drop = FALSE]
}

.dvm_depth_geometry <- function(profile_df) {
    cols <- c("depth_idx", "depth_mid_m", "depth_top_m", "depth_bot_m", "dz_m")
    missing_cols <- setdiff(cols, names(profile_df))
    if (length(missing_cols) > 0) {
        stop("Profile data frame is missing columns: ", paste(missing_cols, collapse = ", "))
    }
    depth_df <- unique(profile_df[, cols])
    depth_df <- depth_df[order(depth_df$depth_idx), ]
    z_mid <- depth_df$depth_mid_m
    z_top <- depth_df$depth_top_m
    z_bot <- depth_df$depth_bot_m
    dz <- depth_df$dz_m
    H_total <- sum(dz)
    list(
        depth = depth_df,
        z_mid = z_mid,
        z_top = z_top,
        z_bot = z_bot,
        dz = dz,
        H_total = H_total,
        w_depth = dz / H_total
    )
}

.dvm_build_resource_spectra_from_cobalt <- function(params, profile_df, dvm_cfg) {
    geom <- .dvm_depth_geometry(profile_df)
    w_pp <- w_full(params)
    dw_pp <- dw_full(params)
    phi <- lapply(names(dvm_cfg$w_group), function(group) {
        log_raw <- - (log(w_pp) - log(dvm_cfg$w_group[[group]]))^2 /
            (2 * dvm_cfg$sigma_group[[group]]^2)
        log_raw <- log_raw - max(log_raw)
        raw <- exp(log_raw)
        denom <- sum(raw * dw_pp)
        if (!is.finite(denom) || denom <= 0) {
            raw[] <- 0
            raw[which.min(abs(log(w_pp) - log(dvm_cfg$w_group[[group]])))] <- 1
            denom <- sum(raw * dw_pp)
        }
        raw / denom
    })
    names(phi) <- names(dvm_cfg$w_group)

    conv <- dvm_cfg$input_to_gN * dvm_cfg$N_to_C * dvm_cfg$C_to_wet
    B_micro <- profile_df$zmicro[match(geom$depth$depth_idx, profile_df$depth_idx)] * conv
    B_meso <- profile_df$zmeso[match(geom$depth$depth_idx, profile_df$depth_idx)] * conv
    b_pp_micro <- outer(B_micro, phi$micro)
    b_pp_meso <- outer(B_meso, phi$meso)
    b_pp_local <- b_pp_micro + b_pp_meso
    n_pp_local <- sweep(b_pp_local, 2, w_pp, "/")
    dimnames(n_pp_local) <- list(depth = geom$depth$depth_idx, w = names(w_pp))
    list(
        B_micro = B_micro,
        B_meso = B_meso,
        phi = phi,
        b_pp_micro = b_pp_micro,
        b_pp_meso = b_pp_meso,
        b_pp_local = b_pp_local,
        n_pp_local = n_pp_local
    )
}

.dvm_effective_npp <- function(n_pp_local, depth_weights, phase_weight_day = 0.5,
                               phase_weight_night = 0.5) {
    n_pp_day_eff <- as.numeric(crossprod(depth_weights, n_pp_local))
    n_pp_night_eff <- as.numeric(crossprod(depth_weights, n_pp_local))
    n_pp_eff <- phase_weight_day * n_pp_day_eff + phase_weight_night * n_pp_night_eff
    if (any(!is.finite(n_pp_eff))) {
        stop("n_pp_eff contains non-finite values; check the local resource spectrum construction.")
    }
    list(
        n_pp_day_eff = n_pp_day_eff,
        n_pp_night_eff = n_pp_night_eff,
        n_pp_eff = n_pp_eff
    )
}

.dvm_integrated_to_local <- function(N_int, p_phase, H_total, dz) {
    out <- array(0, dim = dim(p_phase), dimnames = dimnames(p_phase))
    for (k in seq_along(dz)) {
        out[, , k] <- N_int * H_total * p_phase[, , k] / dz[k]
    }
    out
}

.dvm_light_scalar <- function(light, dvm_cfg) {
    light^dvm_cfg$light_hill /
        (light^dvm_cfg$light_hill + dvm_cfg$light_half_sat^dvm_cfg$light_hill)
}

.dvm_pred_rate_local <- function(params, n_local, n_pp_local, n_other, feeding_level, t) {
    mizerPredRate(params, n = n_local, n_pp = n_pp_local, n_other = n_other,
                  t = t, feeding_level = feeding_level)
}

.dvm_local_fields <- function(params, fish_local, n_pp_local, abiotic, g_prev_local,
                              effort, t) {
    n_species <- dim(fish_local)[1]
    n_w <- dim(fish_local)[2]
    n_z <- dim(fish_local)[3]
    no_w_full <- length(w_full(params))
    pred_rate_local <- array(0, dim = c(n_species, no_w_full, n_z),
                             dimnames = list(species = dimnames(fish_local)[[1]],
                                             prey_size = names(w_full(params)),
                                             depth = dimnames(fish_local)[[3]]))
    pred_mort_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    f_mort_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    mort_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    encounter_raw_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    encounter_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    feeding_level_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    g_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    nu_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))
    f_real_local <- array(0, dim = dim(fish_local), dimnames = dimnames(fish_local))

    cfg <- .tpo2_get_config(params)
    dvm_cfg <- .tpo2_dvm_get_config(params)

    for (k in seq_len(n_z)) {
        T_k <- abiotic$temp[k]
        h_t <- getMaxIntakeRate(params) *
            .tpo2_theta_q10(T_k, cfg$T_ref, cfg$Q10_C)
        M_M <- sweep(outer(params@species_params$p, params@w, function(x, y) y^x),
                     1, params@species_params$ks, "*") *
            .tpo2_theta_q10(T_k, cfg$T_ref, cfg$Q10_M)
        M_A_max <- outer(params@species_params$k, params@w) *
            .tpo2_theta_q10(T_k, cfg$T_ref, cfg$Q10_A)
        G_O <- cfg$a_O * outer(rep(1, nrow(params@species_params)),
                               params@w^cfg$b_O) *
            .tpo2_theta_q10(T_k, cfg$T_ref, cfg$Q10_O)
        encounter_raw_k <- mizerEncounter(params, n = fish_local[, , k],
                                          n_pp = n_pp_local[k, ],
                                          n_other = n_other, t = t)
        phi_L <- .dvm_light_scalar(abiotic$light[k], dvm_cfg)
        encounter_k <- encounter_raw_k * phi_L
        f_k <- encounter_k / (encounter_k + h_t)
        U_k <- .tpo2_u_from_f(f_k)
        M_A_star <- M_A_max * U_k
        D_approx <- M_M + g_prev_local[, , k] * (M_A_star + cfg$epsilon_SDA * h_t * f_k)
        pO2_int_approx <- abiotic$pO2[k] - D_approx / G_O
        pO2_pos <- pmax(pO2_int_approx, 0)
        g_k <- pO2_pos^cfg$h_g / (pO2_pos^cfg$h_g + cfg$K_g^cfg$h_g)
        M_A <- M_A_star * g_k
        C <- h_t * f_k * g_k
        SDA <- cfg$epsilon_SDA * C
        nu_k <- C - SDA - M_M - M_A
        f_real_k <- f_k * g_k
        pred_rate_k <- .dvm_pred_rate_local(params, fish_local[, , k], n_pp_local[k, ],
                                            n_other, feeding_level = f_real_k, t = t)
        pred_mort_k <- mizerPredMort(params, n = fish_local[, , k],
                                     n_pp = n_pp_local[k, ], n_other = n_other,
                                     t = t, pred_rate = pred_rate_k)
        f_mort_k <- mizerFMort(params, n = fish_local[, , k],
                               n_pp = n_pp_local[k, ], n_other = n_other,
                               t = t, effort = effort,
                               e_growth = pmax(nu_k, 0), pred_mort = pred_mort_k)
        mort_k <- mizerMort(params, n = fish_local[, , k],
                            n_pp = n_pp_local[k, ], n_other = n_other,
                            t = t, f_mort = f_mort_k, pred_mort = pred_mort_k)
        encounter_raw_local[, , k] <- encounter_raw_k
        encounter_local[, , k] <- encounter_k
        feeding_level_local[, , k] <- f_k
        g_local[, , k] <- g_k
        nu_local[, , k] <- nu_k
        f_real_local[, , k] <- f_real_k
        pred_rate_local[, , k] <- pred_rate_k
        pred_mort_local[, , k] <- pred_mort_k
        f_mort_local[, , k] <- f_mort_k
        mort_local[, , k] <- mort_k
    }
    list(
        encounter_raw = encounter_raw_local,
        encounter = encounter_local,
        f = feeding_level_local,
        g = g_local,
        nu = nu_local,
        mu = mort_local,
        pred_rate = pred_rate_local,
        pred_mort = pred_mort_local,
        f_mort = f_mort_local,
        f_real = f_real_local
    )
}

.dvm_movement_kernel <- function(params, nu_fore, mu_fore, geometry, dvm_cfg) {
    dims <- dim(nu_fore)
    n_species <- dims[1]
    n_w <- dims[2]
    n_z <- dims[3]
    P <- array(0, dim = c(n_species, n_w, n_z, n_z),
               dimnames = c(dimnames(nu_fore), list(dest_depth = dimnames(nu_fore)[[3]])))
    dist_mat <- abs(outer(geometry$z_mid, geometry$z_mid, "-"))
    Dmax <- outer(rep(1, n_species), params@w^dvm_cfg$dmax_beta) * dvm_cfg$dmax_const
    for (i in seq_len(n_species)) {
        for (j in seq_len(n_w)) {
            nu_vec <- nu_fore[i, j, ]
            mu_vec <- mu_fore[i, j, ]
            nu_rng <- range(nu_vec)
            mu_rng <- range(mu_vec)
            nu_tilde <- if (diff(nu_rng) < dvm_cfg$eps_norm) rep(0, n_z) else {
                (nu_vec - nu_rng[1]) / (diff(nu_rng) + dvm_cfg$eps_norm)
            }
            mu_tilde <- if (diff(mu_rng) < dvm_cfg$eps_norm) rep(0, n_z) else {
                (mu_vec - mu_rng[1]) / (diff(mu_rng) + dvm_cfg$eps_norm)
            }
            for (k in seq_len(n_z)) {
                d_tilde <- dist_mat[k, ] / (Dmax[i, j] + dvm_cfg$eps_norm)
                score <- nu_tilde - mu_tilde - d_tilde
                feasible <- is.finite(score)
                if (isTRUE(dvm_cfg$hard_reachability_cutoff)) {
                    feasible <- feasible & d_tilde <= 1
                    score[!feasible] <- -Inf
                }
                if (!any(feasible)) {
                    P[i, j, k, k] <- 1
                    next
                }
                stable <- score[feasible] - max(score[feasible])
                probs <- exp(dvm_cfg$beta_move * stable)
                probs <- probs / sum(probs)
                uniform <- rep(1 / sum(feasible), sum(feasible))
                probs <- (1 - dvm_cfg$eta_move) * probs + dvm_cfg$eta_move * uniform
                P[i, j, k, feasible] <- probs
            }
        }
    }
    P
}

.dvm_apply_kernel <- function(p_old, P) {
    dims <- dim(p_old)
    out <- array(0, dim = dims, dimnames = dimnames(p_old))
    for (i in seq_len(dims[1])) {
        for (j in seq_len(dims[2])) {
            out[i, j, ] <- as.vector(p_old[i, j, ] %*% P[i, j, , ])
        }
    }
    out
}

.dvm_transport_g <- function(p_old, P, g_old) {
    dims <- dim(g_old)
    out <- array(1, dim = dims, dimnames = dimnames(g_old))
    for (i in seq_len(dims[1])) {
        for (j in seq_len(dims[2])) {
            weights <- p_old[i, j, ]
            for (l in seq_len(dims[3])) {
                kernel_l <- P[i, j, , l]
                denom <- sum(weights * kernel_l)
                if (denom > 0) {
                    out[i, j, l] <- sum(weights * kernel_l * g_old[i, j, ]) / denom
                }
            }
        }
    }
    out
}

.dvm_local_to_phase_effective <- function(p_phase, local_fields) {
    e_phase <- apply(p_phase * local_fields$nu, c(1, 2), sum)
    feeding_phase <- apply(p_phase * local_fields$f_real, c(1, 2), sum)
    mort_phase <- apply(p_phase * local_fields$mu, c(1, 2), sum)
    pred_rate_phase <- apply(local_fields$pred_rate, c(1, 2), sum)
    list(
        e_phase = e_phase,
        feeding_level_phase = feeding_phase,
        mort_phase = mort_phase,
        pred_rate_phase = pred_rate_phase
    )
}

.dvm_phase_to_daily_effective <- function(day_phase, night_phase, dvm_cfg) {
    list(
        e_eff = dvm_cfg$phase_weight_day * day_phase$e_phase +
            dvm_cfg$phase_weight_night * night_phase$e_phase,
        feeding_level_eff = dvm_cfg$phase_weight_day * day_phase$feeding_level_phase +
            dvm_cfg$phase_weight_night * night_phase$feeding_level_phase,
        mort_eff = dvm_cfg$phase_weight_day * day_phase$mort_phase +
            dvm_cfg$phase_weight_night * night_phase$mort_phase,
        pred_rate_eff = dvm_cfg$phase_weight_day * day_phase$pred_rate_phase +
            dvm_cfg$phase_weight_night * night_phase$pred_rate_phase
    )
}

.tpo2_dvm_make_effective_rates <- function(params, t) {
    cfg <- .tpo2_dvm_get_config(params)
    eff <- cfg$effective
    feeding <- eff$feeding_level_eff
    encounter <- feeding / pmax(1 - feeding, cfg$eps_norm) * getMaxIntakeRate(params)
    pred_mort <- pmax(eff$mort_eff - getExtMort(params), 0)
    list(
        encounter = encounter,
        feeding_level = feeding,
        e = eff$e_eff,
        e_repro = mizerERepro(params, n = NULL, n_pp = eff$n_pp_eff,
                              n_other = initialNOther(params), t = t,
                              e = eff$e_eff),
        e_growth = mizerEGrowth(params, n = NULL, n_pp = eff$n_pp_eff,
                                n_other = initialNOther(params), t = t,
                                e_repro = mizerERepro(params, n = NULL,
                                                      n_pp = eff$n_pp_eff,
                                                      n_other = initialNOther(params),
                                                      t = t, e = eff$e_eff),
                                e = eff$e_eff),
        pred_rate = eff$pred_rate_eff,
        pred_mort = pred_mort,
        f_mort = array(0, dim = dim(feeding), dimnames = dimnames(feeding)),
        mort = eff$mort_eff,
        rdi = NULL,
        rdd = NULL,
        resource_mort = mizerResourceMort(params, n = initialN(params), n_pp = eff$n_pp_eff,
                                          n_other = initialNOther(params), t = t,
                                          pred_rate = eff$pred_rate_eff)
    )
}

tpo2DVMMizerRates <- function(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...) {
    cfg <- .tpo2_dvm_get_config(params)
    eff <- cfg$effective
    encounter <- eff$feeding_level_eff / pmax(1 - eff$feeding_level_eff, cfg$eps_norm) *
        getMaxIntakeRate(params)
    e_repro <- rates_fns$ERepro(params, n = n, n_pp = eff$n_pp_eff, n_other = n_other,
                                e = eff$e_eff, t = t, ...)
    e_growth <- rates_fns$EGrowth(params, n = n, n_pp = eff$n_pp_eff, n_other = n_other,
                                  e_repro = e_repro, e = eff$e_eff, t = t, ...)
    pred_mort <- rates_fns$PredMort(params, n = n, n_pp = eff$n_pp_eff, n_other = n_other,
                                    pred_rate = eff$pred_rate_eff, t = t, ...)
    f_mort <- array(0, dim = dim(n), dimnames = dimnames(n))
    mort <- eff$mort_eff
    rdi <- rates_fns$RDI(params, n = n, n_pp = eff$n_pp_eff, n_other = n_other,
                         e_growth = e_growth, mort = mort, e_repro = e_repro,
                         t = t, ...)
    rdd <- rates_fns$RDD(rdi = rdi, species_params = params@species_params,
                         params = params, t = t, ...)
    resource_mort <- rates_fns$ResourceMort(params, n = n, n_pp = eff$n_pp_eff,
                                            n_other = n_other,
                                            pred_rate = eff$pred_rate_eff,
                                            t = t, ...)
    list(
        encounter = encounter,
        feeding_level = eff$feeding_level_eff,
        e = eff$e_eff,
        e_repro = e_repro,
        e_growth = e_growth,
        pred_rate = eff$pred_rate_eff,
        pred_mort = pred_mort,
        f_mort = f_mort,
        mort = mort,
        rdi = rdi,
        rdd = rdd,
        resource_mort = resource_mort
    )
}

#' Enable day/night DVM on top of the temperature-oxygen extension
#'
#' @param params A `MizerParams` object.
#' @param profiles A data frame with depth-resolved forcing columns, or a path
#'   to the CSV profile file.
#' @param site Site selector matching the `site` column.
#' @param scenario Scenario selector matching the `scenario` column.
#' @param tpo2_pars Optional overrides for the temperature-oxygen parameters.
#' @param dvm_pars Optional overrides for DVM and resource mapping parameters.
#' @param p_day0,p_night0 Optional species x size x depth probability arrays.
#' @param g_state_day0,g_state_night0 Optional species x size x depth arrays.
#'
#' @return A modified `MizerParams` object with the DVM wrapper enabled.
#' @export
#' @examples
#' \donttest{
#' params <- newCommunityParams()
#' csv <- system.file("extdata",
#'     "profiles_emblematic_sites_hist_vs_ssp585_long.csv",
#'     package = "mizer")
#' if (nzchar(csv)) {
#'     params <- enable_tpo2_dvm(params, profiles = csv,
#'         site = unique(utils::read.csv(csv)$site)[1], scenario = "hist")
#' }
#' }
enable_tpo2_dvm <- function(params, profiles, site, scenario = "hist",
                            tpo2_pars = list(), dvm_pars = list(),
                            p_day0 = NULL, p_night0 = NULL,
                            g_state_day0 = NULL, g_state_night0 = NULL) {
    params <- validParams(params)
    profile_df <- if (is.character(profiles)) utils::read.csv(profiles) else profiles
    profile_df <- profile_df[profile_df$site == site & profile_df$scenario == scenario, ]
    if (nrow(profile_df) == 0) {
        stop("No profile rows matched the requested site and scenario.")
    }
    profile_df <- .dvm_prepare_profiles(profile_df)
    forcing <- list(
        T_fun = function(t) mean(profile_df$temp_C, na.rm = TRUE),
        pO2_fun = function(t) mean(profile_df$pO2_kPa, na.rm = TRUE)
    )
    resource_dynamics(params) <- "resource_constant"
    params <- enable_tpo2_community(params, forcing = forcing, tpo2_pars = tpo2_pars)
    dvm_cfg <- utils::modifyList(.tpo2_dvm_defaults(), dvm_pars)
    geometry <- .dvm_depth_geometry(profile_df)
    resource <- .dvm_build_resource_spectra_from_cobalt(params, profile_df, dvm_cfg)
    dims <- c(dim(getMaxIntakeRate(params)), length(geometry$dz))
    dn <- c(dimnames(getMaxIntakeRate(params)), list(depth = as.character(geometry$depth$depth_idx)))
    if (is.null(p_day0)) p_day0 <- array(1 / dims[3], dim = dims, dimnames = dn)
    if (is.null(p_night0)) p_night0 <- array(1 / dims[3], dim = dims, dimnames = dn)
    if (is.null(g_state_day0)) g_state_day0 <- array(1, dim = dims, dimnames = dn)
    if (is.null(g_state_night0)) g_state_night0 <- array(1, dim = dims, dimnames = dn)
    .dvm_validate_setup(params, geometry, resource, p_day0, p_night0,
                        g_state_day0, g_state_night0)
    .dvm_normalise_probabilities(p_day0, "p_day0")
    .dvm_normalise_probabilities(p_night0, "p_night0")
    effective_npp <- .dvm_effective_npp(resource$n_pp_local, geometry$w_depth,
                                        dvm_cfg$phase_weight_day,
                                        dvm_cfg$phase_weight_night)
    dvm_cfg$site <- site
    dvm_cfg$scenario <- scenario
    dvm_cfg$profiles <- profile_df
    dvm_cfg$geometry <- geometry
    dvm_cfg$resource <- resource
    dvm_cfg$p_day <- p_day0
    dvm_cfg$p_night <- p_night0
    dvm_cfg$g_state_day <- g_state_day0
    dvm_cfg$g_state_night <- g_state_night0
    dvm_cfg$prev_day <- list(p = p_day0, g = g_state_day0)
    dvm_cfg$prev_night <- list(p = p_night0, g = g_state_night0)
    dvm_cfg$effective <- c(effective_npp, list(
        e_eff = array(0, dim = dim(getMaxIntakeRate(params)),
                      dimnames = dimnames(getMaxIntakeRate(params))),
        feeding_level_eff = array(0, dim = dim(getMaxIntakeRate(params)),
                                  dimnames = dimnames(getMaxIntakeRate(params))),
        mort_eff = array(0, dim = dim(getMaxIntakeRate(params)),
                         dimnames = dimnames(getMaxIntakeRate(params))),
        pred_rate_eff = array(0, dim = c(dim(getMaxIntakeRate(params))[1], length(w_full(params))),
                              dimnames = list(species = rownames(getMaxIntakeRate(params)),
                                              prey_size = names(w_full(params))))
    ))
    other <- other_params(params)
    other$tpo2_dvm <- dvm_cfg
    other_params(params) <- other
    params <- setRateFunction(params, "Rates", "tpo2DVMMizerRates")
    initialNResource(params) <- effective_npp$n_pp_eff
    params
}

#' Project a DVM-enabled temperature-oxygen model
#'
#' @param params A DVM-enabled `MizerParams` object.
#' @param t_max Number of biological days to project.
#' @param dt Time step passed to [project()].
#' @param effort Fishing effort passed to [project()].
#' @param ... Additional arguments passed to [project()].
#'
#' @return A named list containing `sim`, `params_final`, `history` and
#'   `profiles_saved`.
#' @export
project_tpo2_dvm <- function(params, t_max, dt = 0.1, effort = 0, ...) {
    params <- validParams(params)
    dvm_cfg <- .tpo2_dvm_get_config(params)
    day_times <- seq_len(t_max)
    sim_combined <- NULL
    profile_history <- vector("list", length(day_times))
    current_params <- params
    for (step in day_times) {
        cfg <- .tpo2_dvm_get_config(current_params)
        geometry <- cfg$geometry
        resource <- cfg$resource
        N_int <- initialN(current_params)
        n_other <- initialNOther(current_params)
        abiotic_day <- list(temp = cfg$profiles$temp_C[match(geometry$depth$depth_idx, cfg$profiles$depth_idx)],
                            pO2 = cfg$profiles$pO2_kPa[match(geometry$depth$depth_idx, cfg$profiles$depth_idx)],
                            light = cfg$profiles$I_day_rel[match(geometry$depth$depth_idx, cfg$profiles$depth_idx)])
        abiotic_night <- list(temp = abiotic_day$temp,
                              pO2 = abiotic_day$pO2,
                              light = cfg$profiles$I_night_rel[match(geometry$depth$depth_idx, cfg$profiles$depth_idx)])
        fish_prev_day <- .dvm_integrated_to_local(N_int, cfg$prev_day$p, geometry$H_total, geometry$dz)
        forecast_day <- .dvm_local_fields(current_params, fish_prev_day, resource$n_pp_local,
                                          abiotic_day, cfg$prev_day$g, effort, t = step - 1)
        P_day <- .dvm_movement_kernel(current_params, forecast_day$nu, forecast_day$mu,
                                      geometry, cfg)
        p_day_new <- .dvm_apply_kernel(cfg$p_night, P_day)
        g_day_new <- .dvm_transport_g(cfg$p_night, P_day, cfg$g_state_night)
        fish_day_real <- .dvm_integrated_to_local(N_int, p_day_new, geometry$H_total, geometry$dz)
        realized_day <- .dvm_local_fields(current_params, fish_day_real, resource$n_pp_local,
                                          abiotic_day, g_day_new, effort, t = step - 1)

        fish_prev_night <- .dvm_integrated_to_local(N_int, cfg$prev_night$p, geometry$H_total, geometry$dz)
        forecast_night <- .dvm_local_fields(current_params, fish_prev_night, resource$n_pp_local,
                                            abiotic_night, cfg$prev_night$g, effort, t = step - 1)
        P_night <- .dvm_movement_kernel(current_params, forecast_night$nu, forecast_night$mu,
                                        geometry, cfg)
        p_night_new <- .dvm_apply_kernel(p_day_new, P_night)
        g_night_new <- .dvm_transport_g(p_day_new, P_night, g_day_new)
        fish_night_real <- .dvm_integrated_to_local(N_int, p_night_new, geometry$H_total, geometry$dz)
        realized_night <- .dvm_local_fields(current_params, fish_night_real, resource$n_pp_local,
                                            abiotic_night, g_night_new, effort, t = step - 1)

        day_phase <- .dvm_local_to_phase_effective(p_day_new, realized_day)
        night_phase <- .dvm_local_to_phase_effective(p_night_new, realized_night)
        eff_daily <- .dvm_phase_to_daily_effective(day_phase, night_phase, cfg)
        eff_npp <- .dvm_effective_npp(resource$n_pp_local, geometry$w_depth,
                                      cfg$phase_weight_day, cfg$phase_weight_night)
        other <- other_params(current_params)
        other$tpo2_dvm$effective <- c(eff_npp, eff_daily)
        other$tpo2_dvm$p_day <- p_day_new
        other$tpo2_dvm$p_night <- p_night_new
        other$tpo2_dvm$g_state_day <- g_day_new
        other$tpo2_dvm$g_state_night <- g_night_new
        other$tpo2_dvm$prev_day <- list(p = p_day_new, g = realized_day$g)
        other$tpo2_dvm$prev_night <- list(p = p_night_new, g = realized_night$g)
        other_params(current_params) <- other
        initialNResource(current_params) <- eff_npp$n_pp_eff
        sim_step <- project(current_params, t_start = step - 1, t_max = 1, dt = dt,
                            t_save = 1, effort = effort, progress_bar = FALSE, ...)
        current_params <- setInitialValues(current_params, sim_step)
        other <- other_params(current_params)
        other$tpo2_dvm$effective <- c(eff_npp, eff_daily)
        other$tpo2_dvm$p_day <- p_day_new
        other$tpo2_dvm$p_night <- p_night_new
        other$tpo2_dvm$g_state_day <- g_day_new
        other$tpo2_dvm$g_state_night <- g_night_new
        other$tpo2_dvm$prev_day <- list(p = p_day_new, g = realized_day$g)
        other$tpo2_dvm$prev_night <- list(p = p_night_new, g = realized_night$g)
        other_params(current_params) <- other
        initialNResource(current_params) <- eff_npp$n_pp_eff
        profile_history[[step]] <- list(
            p_day = p_day_new,
            p_night = p_night_new,
            forecast_day = forecast_day,
            realized_day = realized_day,
            forecast_night = forecast_night,
            realized_night = realized_night,
            P_day = P_day,
            P_night = P_night,
            effective = c(eff_npp, eff_daily)
        )
        if (is.null(sim_combined)) {
            sim_combined <- sim_step
        } else {
            old_n <- sim_combined@n
            old_pp <- sim_combined@n_pp
            old_effort <- sim_combined@effort
            old_other <- sim_combined@n_other
            merged_times <- c(as.numeric(dimnames(old_n)[[1]]), step)
            sim_new <- MizerSim(current_params, t_dimnames = merged_times)
            n_old <- dim(old_n)[1]
            sim_new@n[1:n_old, , ] <- old_n
            sim_new@n_pp[1:n_old, ] <- old_pp
            sim_new@effort[1:n_old, ] <- old_effort
            sim_new@n_other[1:n_old, ] <- old_other
            sim_new@n[n_old + 1, , ] <- sim_step@n[2, , ]
            sim_new@n_pp[n_old + 1, ] <- sim_step@n_pp[2, ]
            sim_new@effort[n_old + 1, ] <- sim_step@effort[2, ]
            sim_new@n_other[n_old + 1, ] <- sim_step@n_other[2, ]
            sim_combined <- sim_new
        }
    }
    list(sim = sim_combined, params_final = current_params,
         history = profile_history, profiles_saved = profile_history)
}
