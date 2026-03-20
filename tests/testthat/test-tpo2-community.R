make_tpo2_forcing <- function(T = 10, pO2 = 20) {
    list(times = c(0, 1), T_vec = c(T, T), pO2_vec = c(pO2, pO2))
}

test_that("standard mizer is unchanged when tpo2 is disabled", {
    params <- newCommunityParams()
    encounter <- getEncounter(params)
    feeding_level <- getFeedingLevel(params)
    expect_equal(feeding_level, mizerFeedingLevel(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0,
        encounter = encounter
    ))
    expect_equal(getPredRate(params), mizerPredRate(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0,
        feeding_level = feeding_level
    ))
    expect_equal(getEReproAndGrowth(params), mizerEReproAndGrowth(
        params,
        n = initialN(params),
        n_pp = initialNResource(params),
        n_other = initialNOther(params),
        t = 0,
        encounter = encounter,
        feeding_level = feeding_level
    ))
})

test_that("high oxygen at reference temperature stays close in feeding but changes energy via U", {
    params <- newCommunityParams()
    params_tpo2 <- enable_tpo2_community(params, make_tpo2_forcing())
    expect_equal(getFeedingLevel(params_tpo2), getFeedingLevel(params), tolerance = 1e-12)
    e_std <- getEReproAndGrowth(params)
    e_tpo2 <- getEReproAndGrowth(params_tpo2)
    expect_false(isTRUE(all.equal(e_std, e_tpo2, tolerance = 1e-8)))
})

test_that("lower oxygen lowers g, ingestion, active metabolism and predation mortality", {
    params <- newCommunityParams()
    hi <- enable_tpo2_community(params, make_tpo2_forcing(pO2 = 100))
    lo <- enable_tpo2_community(params, make_tpo2_forcing(pO2 = 0.1))
    dhi <- getTPO2Diagnostics(hi)
    dlo <- getTPO2Diagnostics(lo)
    expect_true(all(dlo$g <= dhi$g))
    expect_true(all(dlo$C <= dhi$C))
    expect_true(all(dlo$M_A <= dhi$M_A))
    pm_hi <- getPredMort(hi)
    pm_lo <- getPredMort(lo)
    expect_true(all(pm_lo <= pm_hi))
})

test_that("U has the required shape", {
    expect_equal(.tpo2_u_from_f(0), 0)
    expect_equal(.tpo2_u_from_f(1), 0)
    expect_equal(.tpo2_u_from_f(0.5), 1)
})

test_that("project_tpo2 updates g_prev across steps", {
    params <- enable_tpo2_community(newCommunityParams(), make_tpo2_forcing(pO2 = 0.5))
    g0 <- other_params(params)$tpo2$g_prev
    out <- project_tpo2(params, t_max = 0.2, dt = 0.1, effort = 0)
    g1 <- other_params(out$params_final)$tpo2$g_prev
    expect_false(isTRUE(all.equal(g0, g1)))
    expect_identical(dim(out$saved_g)[2:3], dim(g0))
})

test_that("diagnostics keep species x size dimensions", {
    params <- enable_tpo2_community(newCommunityParams(), make_tpo2_forcing())
    diag <- getTPO2Diagnostics(params)
    expected_dim <- dim(getMaxIntakeRate(params))
    for (nm in names(diag)) {
        expect_identical(dim(diag[[nm]]), expected_dim)
    }
})
