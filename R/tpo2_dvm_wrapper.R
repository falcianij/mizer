
# Temperature-oxygen-light diel vertical migration wrapper for mizer.

.tpo2_dvm_defaults <- function() {
    utils::modifyList(.tpo2_defaults(), list(
        a_O = 2,
        b_O = 0.67,
        h_g = 2,
        phase_weights = c(day = 0.5, night = 0.5),
        beta_move = 5,
        eps_norm = 1e-8,
        light_half_sat = 0.5,
        light_hill = 1,
        dmax_const = 150,
        dmax_beta = 0,
        hard_reachability_cutoff = TRUE,
        save_forecast_profiles = TRUE,
        save_realized_profiles = TRUE
    ))
}

.tpo2_dvm_get_config <- function(params) {
    cfg <- other_params(params)$tpo2_dvm
    if (is.null(cfg) || !isTRUE(cfg$enabled)) {
        stop("The tpo2 DVM wrapper is not enabled for this params object.")
    }
    cfg
}

.dvm_species_size_dimnames <- function(params) {
    dimnames(getMaxIntakeRate(params))
}

.dvm_species_size_depth_dimnames <- function(params, z_mid) {
    c(.dvm_species_size_dimnames(params), list(depth = as.character(z_mid)))
}

.dvm_array3 <- function(value, params, z_mid, name) {
    target_dim <- c(dim(getMaxIntakeRate(params)), length(z_mid))
    target_dimnames <- .dvm_species_size_depth_dimnames(params, z_mid)
    if (length(value) == 1) {
        return(array(value, dim = target_dim, dimnames = target_dimnames))
    }
    if (is.array(value) && identical(dim(value), target_dim)) {
        if (is.null(dimnames(value))) {
            dimnames(value) <- target_dimnames
        }
        return(value)
    }
    stop(name, " must be a scalar or a species x size x depth array.")
}

.dvm_prob_array <- function(value, params, z_mid, name) {
    out <- .dvm_array3(value, params, z_mid, name)
    totals <- apply(out, c(1, 2), sum)
    if (any(abs(totals - 1) > 1e-8)) {
        stop(name, " must sum to 1 over depth for every species and size.")
    }
    out
}

.dvm_depth_weights <- function(z_mid, z_top = NULL, z_bot = NULL) {
    stopifnot(is.numeric(z_mid), length(z_mid) > 0)
    if (!is.null(z_top) || !is.null(z_bot)) {
        if (is.null(z_top) || is.null(z_bot)) {
            stop("z_top and z_bot must either both be NULL or both be supplied.")
        }
        widths <- z_bot - z_top
        if (any(widths <= 0)) {
            stop("Depth cell widths must be strictly positive.")
        }
        return(widths / sum(widths))
    }
    if (length(z_mid) == 1) {
        return(1)
    }
    bounds <- numeric(length(z_mid) + 1)
    bounds[2:length(z_mid)] <- (z_mid[-1] + z_mid[-length(z_mid)]) / 2
    bounds[1] <- z_mid[1] - (bounds[2] - z_mid[1])
    bounds[length(bounds)] <- z_mid[length(z_mid)] +
        (z_mid[length(z_mid)] - bounds[length(bounds) - 1])
    widths <- diff(bounds)
    widths / sum(widths)
}


.dvm_esd_um_to_mass_g <- function(esd_um, density_g_m3 = 1e6) {
    radius_m <- (esd_um * 1e-6) / 2
    (4 / 3) * pi * radius_m^3 * density_g_m3
}

.dvm_bin_overlap_weights <- function(w_full, dw_full, lower, upper) {
    bin_lower <- pmax(w_full - dw_full / 2, .Machine$double.eps)
    bin_upper <- w_full + dw_full / 2
    overlap <- pmax(0, pmin(bin_upper, upper) - pmax(bin_lower, lower))
    if (sum(overlap) <= 0) {
        closest <- which.min(abs(log(w_full) - mean(log(c(lower, upper)))))
        overlap[closest] <- 1
    }
    overlap / sum(overlap)
}

.dvm_profile_resource_spectrum <- function(params, micro_biomass_g_m3,
                                           meso_biomass_g_m3,
                                           resource_pars = list()) {
    rp <- utils::modifyList(list(
        micro_esd_um = c(20, 200),
        meso_esd_um = c(200, 20000),
        biomass_density_g_m3 = 1e6,
        mol_to_g = 12
    ), resource_pars)
    w_full <- params@w_full
    dw_full <- params@dw_full
    micro_range <- .dvm_esd_um_to_mass_g(rp$micro_esd_um, rp$biomass_density_g_m3)
    meso_range <- .dvm_esd_um_to_mass_g(rp$meso_esd_um, rp$biomass_density_g_m3)
    micro_weights <- .dvm_bin_overlap_weights(w_full, dw_full, min(micro_range), max(micro_range))
    meso_weights <- .dvm_bin_overlap_weights(w_full, dw_full, min(meso_range), max(meso_range))
    micro_biomass_by_bin <- micro_biomass_g_m3 * micro_weights
    meso_biomass_by_bin <- meso_biomass_g_m3 * meso_weights
    total_biomass_by_bin <- micro_biomass_by_bin + meso_biomass_by_bin
    total_biomass_by_bin / (w_full * dw_full)
}

.dvm_make_profile_forcing <- function(params, profile_data, profile_pars = list()) {
    needed <- c("depth_mid_m", "depth_top_m", "depth_bot_m", "T_mean", "pO2_mean",
                "zmeso_mean", "zmicro_mean", "I_day_rel", "I_night_rel")
    missing <- needed[!needed %in% names(profile_data)]
    if (length(missing) > 0) {
        stop("profile_data is missing required columns: ", paste(missing, collapse = ", "))
    }
    prof <- profile_data[order(profile_data$depth_mid_m), , drop = FALSE]
    pars <- utils::modifyList(list(mol_to_g = 12), profile_pars)
    z_mid <- prof$depth_mid_m
    z_top <- prof$depth_top_m
    z_bot <- prof$depth_bot_m
    micro_g_m3 <- prof$zmicro_mean * pars$mol_to_g
    meso_g_m3 <- prof$zmeso_mean * pars$mol_to_g
    n_pp_profile <- vapply(seq_len(nrow(prof)), function(i) {
        .dvm_profile_resource_spectrum(params,
                                       micro_biomass_g_m3 = micro_g_m3[i],
                                       meso_biomass_g_m3 = meso_g_m3[i],
                                       resource_pars = pars)
    }, numeric(length(params@w_full)))
    n_pp_profile <- t(n_pp_profile)
    dimnames(n_pp_profile) <- list(depth = as.character(z_mid), resource = names(params@initial_n_pp))
    list(
        z_mid = z_mid,
        z_top = z_top,
        z_bot = z_bot,
        T_day = function(t, depth) prof$T_mean,
        T_night = function(t, depth) prof$T_mean,
        pO2_day = function(t, depth) prof$pO2_mean,
        pO2_night = function(t, depth) prof$pO2_mean,
        L_day = function(t, depth) prof$I_day_rel,
        L_night = function(t, depth) prof$I_night_rel,
        n_pp_day = function(t, depth) n_pp_profile,
        n_pp_night = function(t, depth) n_pp_profile,
        profile_data = prof,
        resource_profile = n_pp_profile
    )
}

.dvm_eval_input <- function(x, t, z_mid, expected_cols = NULL, name) {
    if (is.function(x)) {
        value <- x(t, z_mid)
    } else if (is.array(x)) {
        dims <- dim(x)
        idx <- max(1, min(dims[1], floor(t) + 1))
        if (length(dims) == 2) {
            value <- x[idx, , drop = TRUE]
        } else if (length(dims) == 3) {
            value <- x[idx, , , drop = FALSE]
            value <- matrix(value[1, , ], nrow = dims[2], ncol = dims[3])
        } else {
            stop(name, " must have 2 or 3 dimensions when supplied as an array.")
        }
    } else if (is.numeric(x)) {
        value <- x
    } else {
        stop("Unsupported input for ", name, ".")
    }

    if (is.null(expected_cols)) {
        if (length(value) == 1) {
            value <- rep(value, length(z_mid))
        }
        if (length(value) != length(z_mid)) {
            stop(name, " must evaluate to one value per depth cell.")
        }
        return(as.numeric(value))
    }

    if (is.null(dim(value))) {
        stop(name, " must evaluate to a depth x resource array.")
    }
    if (nrow(value) != length(z_mid) || ncol(value) != expected_cols) {
        stop(name, " must evaluate to a depth x resource array.")
    }
    value
}

.dvm_get_phase_env <- function(params, t, phase) {
    cfg <- .tpo2_dvm_get_config(params)
    list(
        T = .dvm_eval_input(cfg$forcing[[paste0("T_", phase)]], t, cfg$z_mid,
                            name = paste0("T_", phase)),
        pO2 = .dvm_eval_input(cfg$forcing[[paste0("pO2_", phase)]], t, cfg$z_mid,
                              name = paste0("pO2_", phase)),
        L = .dvm_eval_input(cfg$forcing[[paste0("L_", phase)]], t, cfg$z_mid,
                            name = paste0("L_", phase))
    )
}

.dvm_get_phase_resource <- function(params, t, phase) {
    cfg <- .tpo2_dvm_get_config(params)
    out <- .dvm_eval_input(cfg$forcing[[paste0("n_pp_", phase)]], t, cfg$z_mid,
                           expected_cols = length(initialNResource(params)),
                           name = paste0("n_pp_", phase))
    dimnames(out) <- list(depth = as.character(cfg$z_mid),
                          resource = dimnames(initialNResource(params))[[1]])
    out
}

.dvm_local_fish_array <- function(N_integrated, p_phase) {
    sweep(p_phase, c(1, 2), N_integrated, "*", check.margin = FALSE)
}

.dvm_slice_species_size <- function(x, z) {
    out <- x[, , z, drop = FALSE]
    array(out, dim = dim(out)[1:2], dimnames = dimnames(out)[1:2])
}

.dvm_slice_pred_rate <- function(x, z) {
    out <- x[, , z, drop = FALSE]
    array(out, dim = dim(out)[1:2], dimnames = dimnames(out)[1:2])
}

.dvm_light_scalar <- function(L, light_half_sat, light_hill) {
    L_pos <- pmax(L, 0)
    L_pos^light_hill / (L_pos^light_hill + light_half_sat^light_hill)
}

.dvm_temperature_multiplier <- function(T_t, T_ref, Q10) {
    Q10^((T_t - T_ref) / 10)
}

.dvm_phase_metabolism <- function(params, T, pO2) {
    cfg <- .tpo2_dvm_get_config(params)
    ns <- nrow(getMaxIntakeRate(params))
    nw <- ncol(getMaxIntakeRate(params))
    nz <- length(T)
    dims <- c(ns, nw, nz)
    dn <- .dvm_species_size_depth_dimnames(params, cfg$z_mid)
    h_t <- array(0, dim = dims, dimnames = dn)
    M_M <- array(0, dim = dims, dimnames = dn)
    M_A_max <- array(0, dim = dims, dimnames = dn)
    G_O <- array(0, dim = dims, dimnames = dn)
    pO2_arr <- array(0, dim = dims, dimnames = dn)
    base_mm <- sweep(outer(params@species_params$p, params@w, function(x, y) y^x),
                     1, params@species_params$ks, "*")
    base_ma <- outer(params@species_params$k, params@w)
    base_go <- cfg$a_O * outer(rep(1, ns), params@w^cfg$b_O)
    for (z in seq_len(nz)) {
        h_t[, , z] <- getMaxIntakeRate(params) *
            .dvm_temperature_multiplier(T[z], cfg$T_ref, cfg$Q10_C)
        M_M[, , z] <- base_mm * .dvm_temperature_multiplier(T[z], cfg$T_ref, cfg$Q10_M)
        M_A_max[, , z] <- base_ma * .dvm_temperature_multiplier(T[z], cfg$T_ref, cfg$Q10_A)
        G_O[, , z] <- base_go * .dvm_temperature_multiplier(T[z], cfg$T_ref, cfg$Q10_O)
        pO2_arr[, , z] <- pO2[z]
    }
    list(h_t = h_t, M_M = M_M, M_A_max = M_A_max, G_O = G_O, pO2 = pO2_arr)
}

.dvm_local_encounter <- function(params, local_fish, local_resource, t, phase) {
    cfg <- .tpo2_dvm_get_config(params)
    env <- .dvm_get_phase_env(params, t, phase)
    nz <- length(cfg$z_mid)
    out <- array(0, dim = dim(local_fish), dimnames = dimnames(local_fish))
    light <- .dvm_light_scalar(env$L, cfg$light_half_sat, cfg$light_hill)
    n_other <- initialNOther(params)
    for (z in seq_len(nz)) {
        enc <- mizerEncounter(params,
                              n = .dvm_slice_species_size(local_fish, z),
                              n_pp = local_resource[z, ],
                              n_other = n_other,
                              t = t)
        out[, , z] <- enc * light[z]
    }
    out
}

.dvm_local_pred_rate <- function(params, local_fish, f_real, t) {
    nz <- dim(local_fish)[3]
    res <- array(0, dim = c(dim(local_fish)[1], length(params@w_full), nz),
                 dimnames = list(species = dimnames(local_fish)[[1]],
                                 size = dimnames(params@initial_n_pp)[[1]],
                                 depth = dimnames(local_fish)[[3]]))
    n_other <- initialNOther(params)
    for (z in seq_len(nz)) {
        res[, , z] <- mizerPredRate(params,
                                    n = .dvm_slice_species_size(local_fish, z),
                                    n_pp = initialNResource(params),
                                    n_other = n_other,
                                    t = t,
                                    feeding_level = .dvm_slice_species_size(f_real, z))
    }
    res
}

.dvm_local_mortality <- function(params, local_fish, pred_rate, effort, t) {
    nz <- dim(local_fish)[3]
    n_other <- initialNOther(params)
    f_mort <- mizerFMort(params,
                         n = .dvm_slice_species_size(local_fish, 1),
                         n_pp = initialNResource(params),
                         n_other = n_other,
                         t = t,
                         effort = validEffortVector(effort, params),
                         e_growth = array(0, dim = dim(getMaxIntakeRate(params))),
                         pred_mort = array(0, dim = dim(getMaxIntakeRate(params))))
    pred_mort <- array(0, dim = dim(local_fish), dimnames = dimnames(local_fish))
    mort <- array(0, dim = dim(local_fish), dimnames = dimnames(local_fish))
    for (z in seq_len(nz)) {
        pred_mort[, , z] <- mizerPredMort(params,
                                          n = .dvm_slice_species_size(local_fish, z),
                                          n_pp = initialNResource(params),
                                          n_other = n_other,
                                          t = t,
                                          pred_rate = .dvm_slice_pred_rate(pred_rate, z))
        mort[, , z] <- pred_mort[, , z] + params@mu_b + f_mort
    }
    list(pred_mort = pred_mort, mort = mort)
}

.dvm_build_rate_state <- function(params, encounter, env_state, g_prev_phase, effort, local_fish, t) {
    cfg <- .tpo2_dvm_get_config(params)
    f <- encounter / (encounter + env_state$h_t)
    U <- pmin(pmax(4 * f * (1 - f), 0), 1)
    M_A_star <- env_state$M_A_max * U
    D_approx <- env_state$M_M + g_prev_phase * (M_A_star + cfg$epsilon_SDA * env_state$h_t * f)
    pO2_int_approx <- env_state$pO2 - D_approx / env_state$G_O
    pO2_pos <- pmax(pO2_int_approx, 0)
    g <- pO2_pos^cfg$h_g / (pO2_pos^cfg$h_g + cfg$K_g^cfg$h_g)
    M_A <- M_A_star * g
    C <- env_state$h_t * f * g
    SDA <- cfg$epsilon_SDA * C
    nu <- C - SDA - env_state$M_M - M_A
    f_real <- f * g
    pred_rate <- .dvm_local_pred_rate(params, local_fish, f_real, t)
    mort_info <- .dvm_local_mortality(params, local_fish, pred_rate, effort, t)
    list(encounter = encounter, f = f, U = U, M_A_star = M_A_star,
         D_approx = D_approx, pO2_int_approx = pO2_int_approx, g = g,
         M_A = M_A, C = C, SDA = SDA, nu = nu, f_real = f_real,
         pred_rate = pred_rate, pred_mort = mort_info$pred_mort, mu = mort_info$mort)
}

.dvm_local_forecast_rates <- function(params, t, phase, local_fish_prev_same_phase,
                                      local_resource_prev_same_phase, g_prev_phase,
                                      effort = 0) {
    env <- .dvm_get_phase_env(params, t, phase)
    env_state <- .dvm_phase_metabolism(params, env$T, env$pO2)
    encounter <- .dvm_local_encounter(params, local_fish_prev_same_phase,
                                      local_resource_prev_same_phase, t, phase)
    .dvm_build_rate_state(params, encounter, env_state, g_prev_phase,
                          effort, local_fish_prev_same_phase, t)
}

.dvm_movement_kernel <- function(params, nu_fore, mu_fore, phase) {
    cfg <- .tpo2_dvm_get_config(params)
    ns <- dim(nu_fore)[1]
    nw <- dim(nu_fore)[2]
    nz <- dim(nu_fore)[3]
    z_mid <- cfg$z_mid
    P <- array(0, dim = c(ns, nw, nz, nz),
               dimnames = c(dimnames(nu_fore),
                            list(depth_to = dimnames(nu_fore)[[3]])))
    for (i in seq_len(ns)) {
        for (j in seq_len(nw)) {
            nu_vec <- nu_fore[i, j, ]
            mu_vec <- mu_fore[i, j, ]
            nu_rng <- range(nu_vec)
            mu_rng <- range(mu_vec)
            if (abs(diff(nu_rng)) < cfg$eps_norm) {
                nu_tilde <- rep(0, nz)
            } else {
                nu_tilde <- (nu_vec - nu_rng[1]) / (diff(nu_rng) + cfg$eps_norm)
            }
            if (abs(diff(mu_rng)) < cfg$eps_norm) {
                mu_tilde <- rep(0, nz)
            } else {
                mu_tilde <- (mu_vec - mu_rng[1]) / (diff(mu_rng) + cfg$eps_norm)
            }
            Dmax <- cfg$dmax_const * params@w[j]^cfg$dmax_beta
            for (z1 in seq_len(nz)) {
                d_tilde <- abs(z_mid - z_mid[z1]) / (Dmax + cfg$eps_norm)
                score <- nu_tilde - mu_tilde - d_tilde
                if (isTRUE(cfg$hard_reachability_cutoff)) {
                    score[d_tilde > 1] <- -Inf
                }
                shifted <- score - max(score[is.finite(score)])
                probs <- exp(cfg$beta_move * shifted)
                probs[!is.finite(score)] <- 0
                if (sum(probs) <= 0) {
                    probs[z1] <- 1
                }
                P[i, j, z1, ] <- probs / sum(probs)
            }
        }
    }
    P
}

.dvm_apply_kernel <- function(p_old, P) {
    ns <- dim(p_old)[1]
    nw <- dim(p_old)[2]
    nz <- dim(p_old)[3]
    p_new <- array(0, dim = dim(p_old), dimnames = dimnames(p_old))
    for (i in seq_len(ns)) {
        for (j in seq_len(nw)) {
            for (z2 in seq_len(nz)) {
                p_new[i, j, z2] <- sum(p_old[i, j, ] * P[i, j, , z2])
            }
            total <- sum(p_new[i, j, ])
            if (total > 0) {
                p_new[i, j, ] <- p_new[i, j, ] / total
            }
        }
    }
    p_new
}

.dvm_transport_g <- function(g_prev_old, p_old, P) {
    ns <- dim(g_prev_old)[1]
    nw <- dim(g_prev_old)[2]
    nz <- dim(g_prev_old)[3]
    out <- array(1, dim = dim(g_prev_old), dimnames = dimnames(g_prev_old))
    for (i in seq_len(ns)) {
        for (j in seq_len(nw)) {
            for (z2 in seq_len(nz)) {
                weights <- p_old[i, j, ] * P[i, j, , z2]
                denom <- sum(weights)
                out[i, j, z2] <- if (denom > 0) {
                    sum(weights * g_prev_old[i, j, ]) / denom
                } else {
                    1
                }
            }
        }
    }
    out
}

.dvm_local_realized_rates <- function(params, t, phase, local_fish_realized,
                                      local_resource_realized,
                                      g_prev_phase_arrive, effort = 0) {
    env <- .dvm_get_phase_env(params, t, phase)
    env_state <- .dvm_phase_metabolism(params, env$T, env$pO2)
    encounter <- .dvm_local_encounter(params, local_fish_realized,
                                      local_resource_realized, t, phase)
    .dvm_build_rate_state(params, encounter, env_state, g_prev_phase_arrive,
                          effort, local_fish_realized, t)
}

.dvm_collapse_effective_rates <- function(params, day_realized, night_realized,
                                          p_day, p_night, n_pp_day, n_pp_night) {
    cfg <- .tpo2_dvm_get_config(params)
    w_phase <- cfg$phase_weights / sum(cfg$phase_weights)
    e_eff_day <- apply(p_day * day_realized$nu, c(1, 2), sum)
    e_eff_night <- apply(p_night * night_realized$nu, c(1, 2), sum)
    feeding_level_eff_day <- apply(p_day * day_realized$f_real, c(1, 2), sum)
    feeding_level_eff_night <- apply(p_night * night_realized$f_real, c(1, 2), sum)
    mort_eff_day <- apply(p_day * day_realized$mu, c(1, 2), sum)
    mort_eff_night <- apply(p_night * night_realized$mu, c(1, 2), sum)
    pred_rate_eff <- w_phase[["day"]] * apply(day_realized$pred_rate, c(1, 2), sum) +
        w_phase[["night"]] * apply(night_realized$pred_rate, c(1, 2), sum)
    depth_w <- .dvm_depth_weights(cfg$z_mid, cfg$z_top, cfg$z_bot)
    n_pp_day_eff <- colSums(sweep(n_pp_day, 1, depth_w, "*"))
    n_pp_night_eff <- colSums(sweep(n_pp_night, 1, depth_w, "*"))
    list(
        e_eff_day = e_eff_day,
        e_eff_night = e_eff_night,
        e_eff = w_phase[["day"]] * e_eff_day + w_phase[["night"]] * e_eff_night,
        feeding_level_eff_day = feeding_level_eff_day,
        feeding_level_eff_night = feeding_level_eff_night,
        feeding_level_eff = w_phase[["day"]] * feeding_level_eff_day +
            w_phase[["night"]] * feeding_level_eff_night,
        mort_eff_day = mort_eff_day,
        mort_eff_night = mort_eff_night,
        mort_eff = w_phase[["day"]] * mort_eff_day + w_phase[["night"]] * mort_eff_night,
        pred_rate_eff = pred_rate_eff,
        n_pp_day_eff = n_pp_day_eff,
        n_pp_night_eff = n_pp_night_eff,
        n_pp_eff = w_phase[["day"]] * n_pp_day_eff + w_phase[["night"]] * n_pp_night_eff
    )
}

.tpo2_dvm_effective_array <- function(params, name) {
    cfg <- .tpo2_dvm_get_config(params)
    eff <- cfg$effective[[name]]
    if (is.null(eff)) {
        stop("Effective DVM rate '", name, "' has not been populated yet.")
    }
    eff
}

#' Effective DVM feeding level.
#' @inheritParams mizerFeedingLevel
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2DVMFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    .tpo2_dvm_effective_array(params, "feeding_level_eff")
}

#' Effective DVM predation rate.
#' @inheritParams mizerPredRate
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2DVMPredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
    .tpo2_dvm_effective_array(params, "pred_rate_eff")
}

#' Effective DVM total mortality.
#' @inheritParams mizerMort
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2DVMMort <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
    .tpo2_dvm_effective_array(params, "mort_eff")
}

#' Effective DVM energy for reproduction and growth.
#' @inheritParams mizerEReproAndGrowth
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2DVMEReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter,
                                   feeding_level, ...) {
    .tpo2_dvm_effective_array(params, "e_eff")
}

#' Enable the temperature-oxygen DVM wrapper.
#'
#' Enables a diel vertical migration wrapper around the existing temperature-
#' oxygen extension logic. The wrapper owns day and night vertical occupancy and
#' lagged oxygen-limitation state, computes local forecast and realised profiles
#' by depth and phase, collapses them to effective vertically integrated rates,
#' and then lets [project()] advance the integrated fish community state for one
#' short step.
#'
#' @param params A `MizerParams` object.
#' @param dvm_forcing A named list with either explicit `z_mid`, optional
#'   `z_top`/`z_bot`, and phase-specific `T`, `pO2`, `L`, and `n_pp` forcings,
#'   or a `profile_data` data frame with `depth_mid_m`, `depth_top_m`,
#'   `depth_bot_m`, `T_mean`, `pO2_mean`, `zmeso_mean`, `zmicro_mean`,
#'   `I_day_rel`, and `I_night_rel` so the wrapper can derive a depth-resolved
#'   resource spectrum.
#' @param dvm_pars Named list overriding DVM and T-pO2 default parameters.
#' @param init Named list with optional `p_day0`, `p_night0`, `g_prev_day0`, and
#'   `g_prev_night0` arrays.
#'
#' @return A modified `MizerParams` object with the DVM wrapper enabled.
#' @export
#' @examples
#' params <- newCommunityParams(no_w = 8)
#' z_mid <- c(5, 25, 80)
#' nres <- length(initialNResource(params))
#' dvm_forcing <- list(
#'     z_mid = z_mid,
#'     T_day = function(t, depth) c(12, 10, 8),
#'     T_night = function(t, depth) c(11, 9, 8),
#'     pO2_day = function(t, depth) c(8, 6, 4),
#'     pO2_night = function(t, depth) c(8, 6, 4),
#'     L_day = function(t, depth) c(1, 0.4, 0.1),
#'     L_night = function(t, depth) c(0, 0, 0),
#'     n_pp_day = function(t, depth) matrix(rep(initialNResource(params), each = 3),
#'                                         nrow = 3),
#'     n_pp_night = function(t, depth) matrix(rep(initialNResource(params), each = 3),
#'                                           nrow = 3)
#' )
#' params <- enable_tpo2_dvm(params, dvm_forcing)
enable_tpo2_dvm <- function(params, dvm_forcing, dvm_pars = list(), init = list()) {
    params <- validParams(params)
    cfg <- utils::modifyList(.tpo2_dvm_defaults(), dvm_pars)
    cfg$enabled <- TRUE
    if (!is.null(dvm_forcing$profile_data)) {
        dvm_forcing <- utils::modifyList(
            .dvm_make_profile_forcing(params, dvm_forcing$profile_data,
                                      profile_pars = dvm_forcing$profile_pars %||% list()),
            dvm_forcing
        )
    }
    if (is.null(dvm_forcing$z_mid)) {
        stop("dvm_forcing must provide z_mid or profile_data.")
    }
    cfg$z_mid <- dvm_forcing$z_mid
    cfg$z_top <- dvm_forcing$depth_top_m %||% dvm_forcing$z_top
    cfg$z_bot <- dvm_forcing$depth_bot_m %||% dvm_forcing$z_bot
    needed <- c("T_day", "T_night", "pO2_day", "pO2_night", "L_day", "L_night",
                "n_pp_day", "n_pp_night")
    missing <- needed[!needed %in% names(dvm_forcing)]
    if (length(missing) > 0) {
        stop("Missing DVM forcing entries: ", paste(missing, collapse = ", "))
    }
    cfg$forcing <- dvm_forcing[c(needed, "profile_data", "resource_profile")]
    uniform <- array(1 / length(cfg$z_mid),
                     dim = c(dim(getMaxIntakeRate(params)), length(cfg$z_mid)),
                     dimnames = .dvm_species_size_depth_dimnames(params, cfg$z_mid))
    cfg$p_day <- .dvm_prob_array(init$p_day0 %||% uniform, params, cfg$z_mid, "p_day0")
    cfg$p_night <- .dvm_prob_array(init$p_night0 %||% uniform, params, cfg$z_mid, "p_night0")
    cfg$g_prev_day <- .dvm_array3(init$g_prev_day0 %||% 1, params, cfg$z_mid, "g_prev_day0")
    cfg$g_prev_night <- .dvm_array3(init$g_prev_night0 %||% 1, params, cfg$z_mid, "g_prev_night0")
    cfg$effective <- list(
        e_eff = getEReproAndGrowth(params),
        feeding_level_eff = getFeedingLevel(params),
        mort_eff = getMort(params),
        pred_rate_eff = getPredRate(params),
        n_pp_eff = initialNResource(params)
    )
    cfg$history <- NULL
    other <- other_params(params)
    other$tpo2_dvm <- cfg
    other_params(params) <- other
    params <- setRateFunction(params, "FeedingLevel", "tpo2DVMFeedingLevel")
    params <- setRateFunction(params, "PredRate", "tpo2DVMPredRate")
    params <- setRateFunction(params, "Mort", "tpo2DVMMort")
    params <- setRateFunction(params, "EReproAndGrowth", "tpo2DVMEReproAndGrowth")
    params
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.dvm_store_history_array <- function(store, name, index, value, times) {
    if (is.null(store[[name]])) {
        value_dim <- dim(value)
        if (is.null(value_dim)) {
            store[[name]] <- array(NA, dim = c(length(times), length(value)),
                                   dimnames = list(time = as.character(times),
                                                   value = names(value)))
        } else {
            store[[name]] <- array(NA, dim = c(length(times), value_dim),
                                   dimnames = c(list(time = as.character(times)), dimnames(value)))
        }
    }
    idx <- c(list(index), rep(list(quote(expr = )), length(dim(store[[name]])) - 1))
    store[[name]] <- do.call("[<-", c(list(store[[name]]), idx, list(value = value)))
    store
}

#' Project a DVM-enabled model one integrated step at a time.
#'
#' @inheritParams project
#' @param params A DVM-enabled `MizerParams` object.
#'
#' @return A named list with `sim`, `params_final`, and `history`.
#' @export
#' @examples
#' \donttest{
#' params <- newCommunityParams(no_w = 8)
#' z_mid <- c(5, 25, 80)
#' dvm_forcing <- list(
#'     z_mid = z_mid,
#'     T_day = function(t, depth) c(12, 10, 8),
#'     T_night = function(t, depth) c(11, 9, 8),
#'     pO2_day = function(t, depth) c(8, 6, 4),
#'     pO2_night = function(t, depth) c(8, 6, 4),
#'     L_day = function(t, depth) c(1, 0.5, 0.1),
#'     L_night = function(t, depth) c(0, 0, 0),
#'     n_pp_day = function(t, depth) matrix(rep(initialNResource(params), each = 3), nrow = 3),
#'     n_pp_night = function(t, depth) matrix(rep(initialNResource(params), each = 3), nrow = 3)
#' )
#' params <- enable_tpo2_dvm(params, dvm_forcing)
#' out <- project_tpo2_dvm(params, t_max = 2 / 365, dt = 1 / 365, effort = 0)
#' prof <- getTPO2DVMProfiles(out$params_final)
#' matplot(z_mid, t(prof$p_day[1, 1, 1, ]), type = "l")
#' }
project_tpo2_dvm <- function(params, t_max, dt = 1 / 365, effort = 0, ...) {
    params <- validParams(params)
    cfg <- .tpo2_dvm_get_config(params)
    times <- seq(0, t_max, by = dt)
    if (tail(times, 1) < t_max) {
        times <- c(times, t_max)
    }
    history <- list()
    sim_combined <- NULL
    current_params <- params

    for (step in seq_along(times)) {
        cfg <- .tpo2_dvm_get_config(current_params)
        N_integrated <- initialN(current_params)
        t_now <- times[step]
        n_pp_day <- .dvm_get_phase_resource(current_params, t_now, "day")
        n_pp_night <- .dvm_get_phase_resource(current_params, t_now, "night")

        local_day_prev <- .dvm_local_fish_array(N_integrated, cfg$p_day)
        day_fore <- .dvm_local_forecast_rates(current_params, t_now, "day",
                                              local_day_prev, n_pp_day,
                                              cfg$g_prev_day, effort = effort)
        P_day <- .dvm_movement_kernel(current_params, day_fore$nu, day_fore$mu, "day")
        p_day_new <- .dvm_apply_kernel(cfg$p_day, P_day)
        g_prev_day_arrive <- .dvm_transport_g(cfg$g_prev_day, cfg$p_day, P_day)
        local_day_real <- .dvm_local_fish_array(N_integrated, p_day_new)
        day_real <- .dvm_local_realized_rates(current_params, t_now, "day",
                                              local_day_real, n_pp_day,
                                              g_prev_day_arrive, effort = effort)

        local_night_prev <- .dvm_local_fish_array(N_integrated, cfg$p_night)
        night_fore <- .dvm_local_forecast_rates(current_params, t_now, "night",
                                                local_night_prev, n_pp_night,
                                                cfg$g_prev_night, effort = effort)
        P_night <- .dvm_movement_kernel(current_params, night_fore$nu, night_fore$mu, "night")
        p_night_new <- .dvm_apply_kernel(cfg$p_night, P_night)
        g_prev_night_arrive <- .dvm_transport_g(cfg$g_prev_night, cfg$p_night, P_night)
        local_night_real <- .dvm_local_fish_array(N_integrated, p_night_new)
        night_real <- .dvm_local_realized_rates(current_params, t_now, "night",
                                                local_night_real, n_pp_night,
                                                g_prev_night_arrive, effort = effort)

        effective <- .dvm_collapse_effective_rates(current_params, day_real, night_real,
                                                   p_day_new, p_night_new,
                                                   n_pp_day, n_pp_night)
        other <- other_params(current_params)
        other$tpo2_dvm$p_day <- p_day_new
        other$tpo2_dvm$p_night <- p_night_new
        other$tpo2_dvm$g_prev_day <- day_real$g
        other$tpo2_dvm$g_prev_night <- night_real$g
        other$tpo2_dvm$effective <- effective
        other_params(current_params) <- other
        initialNResource(current_params) <- effective$n_pp_eff

        history <- .dvm_store_history_array(history, "p_day", step, p_day_new, times)
        history <- .dvm_store_history_array(history, "p_night", step, p_night_new, times)
        history <- .dvm_store_history_array(history, "g_prev_day", step, day_real$g, times)
        history <- .dvm_store_history_array(history, "g_prev_night", step, night_real$g, times)
        if (isTRUE(cfg$save_forecast_profiles)) {
            for (nm in c("nu", "mu", "f", "g")) {
                history <- .dvm_store_history_array(history, paste0(nm, "_fore_day"),
                                                    step, day_fore[[nm]], times)
                history <- .dvm_store_history_array(history, paste0(nm, "_fore_night"),
                                                    step, night_fore[[nm]], times)
            }
        }
        if (isTRUE(cfg$save_realized_profiles)) {
            for (nm in c("nu", "mu", "f_real", "g")) {
                out_nm <- sub("f_real", "f", nm)
                history <- .dvm_store_history_array(history, paste0(out_nm, "_real_day"),
                                                    step, day_real[[nm]], times)
                history <- .dvm_store_history_array(history, paste0(out_nm, "_real_night"),
                                                    step, night_real[[nm]], times)
            }
        }
        for (nm in c("e_eff", "feeding_level_eff", "mort_eff", "pred_rate_eff", "n_pp_eff")) {
            history <- .dvm_store_history_array(history, nm, step, effective[[nm]], times)
        }

        if (step < length(times)) {
            interval <- times[step + 1] - times[step]
            sim_step <- project(current_params, t_start = times[step], t_max = interval,
                                dt = dt, t_save = interval, effort = effort,
                                progress_bar = FALSE, ...)
            if (is.null(sim_combined)) {
                sim_combined <- sim_step
            } else {
                all_times <- c(as.numeric(dimnames(sim_combined@n)[[1]]),
                               as.numeric(dimnames(sim_step@n)[[1]])[2])
                sim_new <- MizerSim(current_params, t_dimnames = all_times)
                n_old <- dim(sim_combined@n)[1]
                sim_new@n[1:n_old, , ] <- sim_combined@n
                sim_new@n_pp[1:n_old, ] <- sim_combined@n_pp
                sim_new@n_other[1:n_old, ] <- sim_combined@n_other
                sim_new@effort[1:n_old, ] <- sim_combined@effort
                sim_new@n[n_old + 1, , ] <- sim_step@n[2, , ]
                sim_new@n_pp[n_old + 1, ] <- sim_step@n_pp[2, ]
                sim_new@n_other[n_old + 1, ] <- sim_step@n_other[2, ]
                sim_new@effort[n_old + 1, ] <- sim_step@effort[2, ]
                sim_combined <- sim_new
            }
            current_params <- setInitialValues(current_params, sim_step)
        }
    }

    other <- other_params(current_params)
    other$tpo2_dvm$history <- history
    other_params(current_params) <- other
    list(sim = sim_combined, params_final = current_params, history = history)
}

#' Retrieve saved T-pO2-DVM profiles.
#'
#' @param params A DVM-enabled `MizerParams` object.
#' @param time_index Optional integer vector selecting saved times.
#'
#' @return A named list of saved profile and effective-rate arrays.
#' @export
getTPO2DVMProfiles <- function(params, time_index = NULL) {
    cfg <- .tpo2_dvm_get_config(validParams(params))
    history <- cfg$history
    if (is.null(history)) {
        stop("No DVM history has been saved on this params object yet.")
    }
    if (is.null(time_index)) {
        return(history)
    }
    lapply(history, function(x) {
        idx <- c(list(time_index), rep(list(quote(expr = )), length(dim(x)) - 1),
                 list(drop = FALSE))
        do.call("[", c(list(x), idx))
    })
}
