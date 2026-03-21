% Generated manually.
\name{tpo2DVMFeedingLevel}
\alias{tpo2DVMFeedingLevel}
\alias{tpo2DVMPredRate}
\alias{tpo2DVMMort}
\alias{tpo2DVMEReproAndGrowth}
\title{Effective DVM rate overrides}
\usage{
tpo2DVMFeedingLevel(params, n, n_pp, n_other, t, encounter, ...)

tpo2DVMPredRate(params, n, n_pp, n_other, t, feeding_level, ...)

tpo2DVMMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...)

tpo2DVMEReproAndGrowth(params, n, n_pp, n_other, t, encounter, feeding_level, ...)
}
\arguments{
  \item{params}{A \code{MizerParams} object with the DVM wrapper enabled.}
  \item{n, n_pp, n_other, t, encounter, feeding_level, f_mort, pred_mort, ...}{Standard mizer rate-function arguments.}
}
\value{
The effective vertically integrated rate array stored in \code{other_params(params)$tpo2_dvm$effective}.
}
\description{
Internal rate functions used by the DVM wrapper to make the integrated \code{project()} step consume the effective feeding-level, predation-rate, mortality, and energy arrays computed from local day and night profiles.
}
\keyword{internal}
