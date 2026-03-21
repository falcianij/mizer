% Generated manually.
\name{enable_tpo2_dvm}
\alias{enable_tpo2_dvm}
\alias{project_tpo2_dvm}
\alias{getTPO2DVMProfiles}
\title{Temperature-oxygen diel vertical migration wrapper}
\usage{
enable_tpo2_dvm(params, dvm_forcing, dvm_pars = list(), init = list())

project_tpo2_dvm(params, t_max, dt = 1/365, effort = 0, ...)

getTPO2DVMProfiles(params, time_index = NULL)
}
\arguments{
  \item{params}{A \code{MizerParams} object.}
  \item{dvm_forcing}{A named list with either explicit \code{z_mid}, optional \code{z_top}/\code{z_bot}, and phase-specific \code{T}, \code{pO2}, \code{L}, and \code{n_pp} forcings, or a \code{profile_data} data frame with \code{depth_mid_m}, \code{depth_top_m}, \code{depth_bot_m}, \code{T_mean}, \code{pO2_mean}, \code{zmeso_mean}, \code{zmicro_mean}, \code{I_day_rel}, and \code{I_night_rel} so the wrapper can derive a depth-resolved resource spectrum.}
  \item{dvm_pars}{Named list overriding DVM and T-pO2 defaults.}
  \item{init}{Named list with optional initial \code{p_day0}, \code{p_night0}, \code{g_prev_day0}, and \code{g_prev_night0}.}
  \item{t_max}{Projection length in years.}
  \item{dt}{Integrated time step passed to \code{project()}.}
  \item{effort}{Fishing effort passed to \code{project()}.}
  \item{...}{Additional arguments passed to \code{project()}.}
  \item{time_index}{Optional saved-time indices to extract.}
}
\value{
\code{enable_tpo2_dvm()} returns a modified \code{MizerParams} object. \code{project_tpo2_dvm()} returns a list with \code{sim}, \code{params_final}, and \code{history}. \code{getTPO2DVMProfiles()} returns the saved DVM profile history.
}
\description{
Utilities for enabling and projecting the temperature-oxygen-light diel vertical migration wrapper and retrieving the saved day/night profile history.
}
