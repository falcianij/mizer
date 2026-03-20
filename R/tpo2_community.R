# Temperature-oxygen community extension for mizer.

.tpo2_defaults <- function() {
    list(
        enabled = TRUE,
        T_ref = 10,
        Q10_M = 2,
        Q10_C = 2,
        Q10_A = 1.2,
        Q10_O = 1.5,
        a_O = 5,
        b_O = 0.67,
        epsilon_SDA = 0.2,
        K_g = 2,
        h_g = 1
    )
}

.tpo2_get_config <- function(params) {
    cfg <- other_params(params)$tpo2
    if (is.null(cfg) || !isTRUE(cfg$enabled)) {
        stop("The tpo2 community extension is not enabled for this params object.")
    }
    cfg
}

.tpo2_make_dimnames <- function(params) {
    dimnames(getMaxIntakeRate(params))
}

.tpo2_check_species_size_array <- function(params, x, name) {
    if (!is.array(x) || !identical(dim(x), dim(getMaxIntakeRate(params)))) {
        stop(name, " must be a species x size array matching getMaxIntakeRate(params).")
    }
    invisible(TRUE)
}

.tpo2_as_species_size_array <- function(params, value, name) {
    target_dim <- dim(getMaxIntakeRate(params))
    target_dimnames <- .tpo2_make_dimnames(params)
    if (length(value) == 1) {
        out <- array(value, dim = target_dim, dimnames = target_dimnames)
        return(out)
    }
    if (is.array(value) && identical(dim(value), target_dim)) {
        if (is.null(dimnames(value))) {
            dimnames(value) <- target_dimnames
        }
        return(value)
    }
    stop(name, " must be a scalar or a species x size array.")
}

.tpo2_eval_forcing <- function(x, t, name) {
    if (is.function(x)) {
        return(x(t))
    }
    if (is.list(x) && !is.null(x$fun) && is.function(x$fun)) {
        return(x$fun(t))
    }
    if (is.list(x) && !is.null(x$vec) && !is.null(x$times)) {
        return(stats::approx(x = x$times, y = x$vec, xout = t, rule = 2)$y)
    }
    stop("forcing must provide ", name, "_fun(t) or ", name, "_vec together with times.")
}

.tpo2_make_forcing <- function(forcing) {
    if (!is.list(forcing)) {
        stop("forcing must be a list.")
    }
    out <- list()
    if (!is.null(forcing$T_fun)) {
        out$T <- list(fun = forcing$T_fun)
    } else if (!is.null(forcing$T_vec)) {
        if (is.null(forcing$times)) {
            stop("forcing$times must be supplied with forcing$T_vec.")
        }
        out$T <- list(vec = forcing$T_vec, times = forcing$times)
    } else {
        stop("forcing must provide T_fun or T_vec.")
    }
    if (!is.null(forcing$pO2_fun)) {
        out$pO2 <- list(fun = forcing$pO2_fun)
    } else if (!is.null(forcing$pO2_vec)) {
        if (is.null(forcing$times)) {
            stop("forcing$times must be supplied with forcing$pO2_vec.")
        }
        out$pO2 <- list(vec = forcing$pO2_vec, times = forcing$times)
    } else {
        stop("forcing must provide pO2_fun or pO2_vec.")
    }
    out
}

.tpo2_get_state <- function(params, n, n_pp, n_other, t,
                            encounter = NULL, feeding_level = NULL) {
    cfg <- .tpo2_get_config(params)
    if (is.null(encounter)) {
        encounter <- getEncounter(params, n = n, n_pp = n_pp, n_other = n_other,
                                  t = t)
    }
    if (is.null(feeding_level)) {
        feeding_level <- encounter / (encounter + .tpo2_h(params, t))
    }
    h_t <- .tpo2_h(params, t)
    U <- .tpo2_u_from_f(feeding_level)
    M_A_star <- .tpo2_ma_max(params, t) * U
    g_info <- .tpo2_g(params, t = t, f = feeding_level, h_t = h_t,
                      M_A_star = M_A_star)
    list(
        encounter = encounter,
        f = feeding_level,
        U = U,
        h_t = h_t,
        M_M = g_info$M_M,
        M_A_star = M_A_star,
        G_O = g_info$G_O,
        D_approx = g_info$D_approx,
        pO2_int_approx = g_info$pO2_int_approx,
        g = g_info$g,
        M_A = M_A_star * g_info$g,
        C = h_t * feeding_level * g_info$g,
        SDA = cfg$epsilon_SDA * h_t * feeding_level * g_info$g,
        f_real = feeding_level * g_info$g
    )
}

#' Enable the temperature-oxygen community extension
#'
#' Enables a lightweight temperature-oxygen (tpo2) extension for an existing
#' `MizerParams` object. The extension keeps the standard mizer encounter logic,
#' but replaces feeding-level, predation-rate and energy calculations with
#' temperature- and oxygen-aware rate functions. Standard mizer behaviour remains
#' unchanged unless this function is called.
#'
#' The extension stores all of its settings in `other_params(params)$tpo2`,
#' including user-editable defaults for temperature and oxygen parameters,
#' forcing data, and the lagged oxygen downregulation scalar `g_prev`.
#'
#' @param params A `MizerParams` object.
#' @param forcing A list supplying either `T_fun(t)` or `T_vec` with `times`,
#'   and either `pO2_fun(t)` or `pO2_vec` with `times`.
#' @param tpo2_pars Named list of parameter overrides. Supported entries are
#'   `T_ref`, `Q10_M`, `Q10_C`, `Q10_A`, `Q10_O`, `a_O`, `b_O`,
#'   `epsilon_SDA`, `K_g` and `h_g`.
#' @param g_init Optional initial lagged oxygen downregulation scalar. Must be
#'   `NULL`, a scalar, or a species x size array. If `NULL`, it is initialised
#'   to 1 everywhere.
#'
#' @return A modified `MizerParams` object with `other_params(params)$tpo2`
#'   populated and custom rate functions registered.
#' @export
#' @examples
#' params <- newCommunityParams()
#' forcing <- list(
#'     times = c(0, 5),
#'     T_vec = c(10, 10),
#'     pO2_vec = c(20, 20)
#' )
#' params_tpo2 <- enable_tpo2_community(params, forcing = forcing)
enable_tpo2_community <- function(params, forcing, tpo2_pars = list(),
                                  g_init = NULL) {
    params <- validParams(params)
    cfg <- utils::modifyList(.tpo2_defaults(), tpo2_pars)
    cfg$forcing <- .tpo2_make_forcing(forcing)
    if (is.null(g_init)) {
        cfg$g_prev <- array(1, dim = dim(getMaxIntakeRate(params)),
                            dimnames = .tpo2_make_dimnames(params))
    } else {
        cfg$g_prev <- .tpo2_as_species_size_array(params, g_init, "g_init")
    }
    other <- other_params(params)
    other$tpo2 <- cfg
    other_params(params) <- other
    params <- setRateFunction(params, "FeedingLevel", "tpo2FeedingLevel")
    params <- setRateFunction(params, "PredRate", "tpo2PredRate")
    params <- setRateFunction(params, "EReproAndGrowth", "tpo2EReproAndGrowth")
    params
}

#' Get temperature and oxygen forcing
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#'
#' @return A list with scalar entries `T_t` and `pO2_env_t`.
#' @keywords internal
.tpo2_get_env <- function(params, t) {
    cfg <- .tpo2_get_config(params)
    list(
        T_t = .tpo2_eval_forcing(cfg$forcing$T, t, "T"),
        pO2_env_t = .tpo2_eval_forcing(cfg$forcing$pO2, t, "pO2")
    )
}

#' Q10 temperature multiplier
#'
#' @param T_t Temperature.
#' @param T_ref Reference temperature.
#' @param Q10 Q10 coefficient.
#'
#' @return Numeric Q10 multiplier.
#' @keywords internal
.tpo2_theta_q10 <- function(T_t, T_ref, Q10) {
    Q10^((T_t - T_ref) / 10)
}

#' Temperature-adjusted maximum intake rate
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#'
#' @return A species x size array.
#' @keywords internal
.tpo2_h <- function(params, t) {
    cfg <- .tpo2_get_config(params)
    env <- .tpo2_get_env(params, t)
    getMaxIntakeRate(params) *
        .tpo2_theta_q10(env$T_t, cfg$T_ref, cfg$Q10_C)
}

#' Temperature-adjusted standard metabolism
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#'
#' @return A species x size array.
#' @keywords internal
.tpo2_mm <- function(params, t) {
    cfg <- .tpo2_get_config(params)
    env <- .tpo2_get_env(params, t)
    sweep(outer(params@species_params$p, params@w, function(x, y) y^x),
          1, params@species_params$ks, "*") *
        .tpo2_theta_q10(env$T_t, cfg$T_ref, cfg$Q10_M)
}

#' Temperature-adjusted active-metabolism asymptote
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#'
#' @return A species x size array.
#' @keywords internal
.tpo2_ma_max <- function(params, t) {
    cfg <- .tpo2_get_config(params)
    env <- .tpo2_get_env(params, t)
    outer(params@species_params$k, params@w) *
        .tpo2_theta_q10(env$T_t, cfg$T_ref, cfg$Q10_A)
}

#' Oxygen conductance term
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#'
#' @return A species x size array.
#' @keywords internal
.tpo2_go <- function(params, t) {
    cfg <- .tpo2_get_config(params)
    env <- .tpo2_get_env(params, t)
    cfg$a_O * outer(rep(1, nrow(params@species_params)), params@w^cfg$b_O) *
        .tpo2_theta_q10(env$T_t, cfg$T_ref, cfg$Q10_O)
}

#' Intended activity scalar from feeding level
#'
#' @param f Feeding level.
#'
#' @return Feeding-level-shaped intended activity scalar clipped to [0, 1].
#' @keywords internal
.tpo2_u_from_f <- function(f) {
    pmin(pmax(4 * f * (1 - f), 0), 1)
}

#' Oxygen downregulation helper
#'
#' @param params A `MizerParams` object with the tpo2 extension enabled.
#' @param t Numeric time.
#' @param f Raw feeding level.
#' @param h_t Temperature-adjusted maximum intake rate.
#' @param M_A_star Intended active metabolism.
#'
#' @return A list with `g`, `D_approx`, `pO2_int_approx`, `M_M` and `G_O`.
#' @keywords internal
.tpo2_g <- function(params, t, f, h_t, M_A_star) {
    cfg <- .tpo2_get_config(params)
    .tpo2_check_species_size_array(params, f, "f")
    .tpo2_check_species_size_array(params, h_t, "h_t")
    .tpo2_check_species_size_array(params, M_A_star, "M_A_star")
    env <- .tpo2_get_env(params, t)
    M_M <- .tpo2_mm(params, t)
    G_O <- .tpo2_go(params, t)
    D_approx <- M_M + cfg$g_prev * (M_A_star + cfg$epsilon_SDA * h_t * f)
    pO2_int_approx <- env$pO2_env_t - D_approx / G_O
    pO2_pos <- pmax(pO2_int_approx, 0)
    g <- pO2_pos^cfg$h_g / (pO2_pos^cfg$h_g + cfg$K_g^cfg$h_g)
    list(
        g = g,
        D_approx = D_approx,
        pO2_int_approx = pO2_int_approx,
        M_M = M_M,
        G_O = G_O
    )
}

#' Feeding level with temperature-adjusted intake rate
#'
#' @inheritParams mizerFeedingLevel
#' @return A species x size array with the raw feeding level.
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2FeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    h_t <- .tpo2_h(params, t)
    encounter / (encounter + h_t)
}

#' Predation rate with oxygen-downregulated realised feeding level
#'
#' @inheritParams mizerPredRate
#' @return A predator species x prey size array.
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2PredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
    h_t <- .tpo2_h(params, t)
    M_A_star <- .tpo2_ma_max(params, t) * .tpo2_u_from_f(feeding_level)
    g_info <- .tpo2_g(params, t = t, f = feeding_level, h_t = h_t,
                      M_A_star = M_A_star)
    f_real <- feeding_level * g_info$g
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)

    if (!is.null(comment(params@pred_kernel))) {
        n_total_in_size_bins <- sweep(n, 2, params@dw, '*', check.margin = FALSE)
        pred_rate <- sweep(params@pred_kernel, c(1, 2),
                           (1 - f_real) * params@search_vol *
                               n_total_in_size_bins,
                           "*", check.margin = FALSE)
        pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
        return(pred_rate)
    }

    idx_sp <- (no_w_full - no_w + 1):no_w_full
    Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
    Q[, idx_sp] <- sweep((1 - f_real) * params@search_vol * n, 2,
                         params@dw, "*")

    pred_rate <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_p) *
        mvfft(base::t(Q)), inverse = TRUE))) / no_w_full
    pred_rate[pred_rate < 1e-18] <- 0
    pred_rate * params@ft_mask
}

#' Energy available for reproduction and growth under the tpo2 extension
#'
#' @inheritParams mizerEReproAndGrowth
#' @return A species x size array with energy available for reproduction and growth.
#' @export
#' @family mizer rate functions
#' @keywords internal
tpo2EReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter,
                                feeding_level, ...) {
    cfg <- .tpo2_get_config(params)
    h_t <- .tpo2_h(params, t)
    U <- .tpo2_u_from_f(feeding_level)
    M_A_star <- .tpo2_ma_max(params, t) * U
    g_info <- .tpo2_g(params, t = t, f = feeding_level, h_t = h_t,
                      M_A_star = M_A_star)
    M_A <- M_A_star * g_info$g
    C <- h_t * feeding_level * g_info$g
    SDA <- cfg$epsilon_SDA * C
    C - SDA - g_info$M_M - M_A
}

#' Project a tpo2-enabled model while updating lagged oxygen limitation
#'
#' Runs [project()] one saved interval at a time, recomputes the end-of-step
#' oxygen downregulation scalar `g`, writes it into `other_params(params)$tpo2`
#' for the next interval, and returns the combined simulation together with the
#' final parameter object and the saved `g` history.
#'
#' @param params A tpo2-enabled `MizerParams` object.
#' @param t_max Length of projection in years.
#' @param dt Time step passed to [project()].
#' @param effort Fishing effort passed to [project()].
#' @param ... Additional arguments passed to [project()].
#'
#' @return A named list with entries `sim`, `params_final` and `saved_g`.
#' @export
#' @examples
#' params <- newCommunityParams()
#' forcing <- list(times = c(0, 1), T_vec = c(10, 10), pO2_vec = c(20, 20))
#' params <- enable_tpo2_community(params, forcing)
#' out <- project_tpo2(params, t_max = 0.2, dt = 0.1, effort = 0)
project_tpo2 <- function(params, t_max, dt = 0.1, effort = 0, ...) {
    params <- validParams(params)
    cfg <- .tpo2_get_config(params)
    times <- seq(0, t_max, by = dt)
    if (tail(times, 1) < t_max) {
        times <- c(times, t_max)
    }
    sim_combined <- NULL
    saved_g <- array(NA,
                     dim = c(length(times), dim(cfg$g_prev)),
                     dimnames = c(list(time = as.character(times)),
                                  dimnames(cfg$g_prev)))
    saved_g[1, , ] <- cfg$g_prev

    current_params <- params
    for (i in 2:length(times)) {
        interval <- times[i] - times[i - 1]
        sim_step <- project(current_params, t_start = times[i - 1],
                            t_max = interval, dt = dt,
                            t_save = interval, effort = effort,
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
        n_end <- finalN(sim_step)
        n_pp_end <- finalNResource(sim_step)
        n_other_end <- finalNOther(sim_step)
        end_time <- times[i]
        encounter <- getEncounter(current_params, n = n_end, n_pp = n_pp_end,
                                  n_other = n_other_end, t = end_time)
        f <- tpo2FeedingLevel(current_params, n = n_end, n_pp = n_pp_end,
                              n_other = n_other_end, t = end_time,
                              encounter = encounter)
        h_t <- .tpo2_h(current_params, end_time)
        M_A_star <- .tpo2_ma_max(current_params, end_time) * .tpo2_u_from_f(f)
        g_new <- .tpo2_g(current_params, t = end_time, f = f, h_t = h_t,
                         M_A_star = M_A_star)$g
        other <- other_params(current_params)
        other$tpo2$g_prev <- g_new
        other_params(current_params) <- other
        current_params <- setInitialValues(current_params, sim_step)
        other <- other_params(current_params)
        other$tpo2$g_prev <- g_new
        other_params(current_params) <- other
        saved_g[i, , ] <- g_new
    }
    list(sim = sim_combined, params_final = current_params, saved_g = saved_g)
}

#' Diagnostics for the tpo2 community extension
#'
#' Computes a named set of species x size diagnostic arrays for the current
#' model state.
#'
#' @param params A tpo2-enabled `MizerParams` object.
#' @param n Species abundance array. Defaults to [initialN()].
#' @param n_pp Resource abundance vector. Defaults to [initialNResource()].
#' @param n_other Other ecosystem components. Defaults to [initialNOther()].
#' @param t Numeric time.
#'
#' @return A named list containing `encounter`, `f`, `U`, `h_t`, `M_M`,
#'   `M_A_star`, `G_O`, `D_approx`, `pO2_int_approx`, `g`, `M_A`, `C`, `SDA`
#'   and `f_real`.
#' @export
#' @examples
#' params <- newCommunityParams()
#' forcing <- list(times = c(0, 1), T_vec = c(10, 10), pO2_vec = c(20, 20))
#' params <- enable_tpo2_community(params, forcing)
#' diag <- getTPO2Diagnostics(params)
#' names(diag)
getTPO2Diagnostics <- function(params, n = initialN(params),
                               n_pp = initialNResource(params),
                               n_other = initialNOther(params), t = 0) {
    params <- validParams(params)
    encounter <- getEncounter(params, n = n, n_pp = n_pp, n_other = n_other,
                              t = t)
    feeding_level <- tpo2FeedingLevel(params, n = n, n_pp = n_pp,
                                      n_other = n_other, t = t,
                                      encounter = encounter)
    .tpo2_get_state(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                    encounter = encounter, feeding_level = feeding_level)
}
