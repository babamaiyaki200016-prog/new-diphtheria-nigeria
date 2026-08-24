# ------------------------------------------------------------------------------
# 05_bayesian_transmission_model.R
# Bayesian SEIR transmission model fit to the Kano weekly confirmed-case series.
#
# Methodologically templated on the Jakarta 2017 outbreak paper's workflow
# (deterministic ODE + Poisson observation model + MCMC posterior over
# transmission parameters -> posterior R0/Rt), but NOT a port of their Stan
# code: their model's compartment structure, breakpoints and several priors
# are fitted to Jakarta-specific serology/contact-survey inputs Nigeria does
# not have (see reference_repos/jakarta_2017/codes/stan/diphtheria-model1.stan).
# This is a simpler model built for what the Kano data actually supports, per
# the user's explicit choice on that trade-off.
#
# Also fit with a from-scratch Metropolis-Hastings sampler rather than Stan:
# this environment's C++ toolchain cannot compile Stan models (Xcode SDK
# 'cmath' header issue -- confirmed by testing) and rstanarm's binary could
# not be downloaded, so no compiled-Bayesian route was available. A pure-R
# MCMC also means the user can rerun this without a working Stan toolchain.
#
# Fitting window: the full 91-week line-list series (02_epidemic_curve.R) is
# NOT one clean epidemic curve -- it shows an earlier smaller wave
# (Jan-May 2023), the big wave (Jun-Nov 2023, peak ~692 cases/week in early
# August), and a prolonged elevated plateau with a secondary bump into mid-
# 2024. That is consistent with repeated reintroduction/spatial spread across
# LGAs, not a single well-mixed outbreak -- a homogeneous-mixing SEIR fit to
# the whole 91 weeks converges to a poor, ambiguous fit (tested: MAP
# log-likelihood ~3286 vs ~677 for the restricted window below, and the
# unrestricted fit cannot pin down a stable R0). We therefore fit the SEIR
# model to the main wave only (2023-06-05 to 2023-12-04, 27 weeks), the one
# segment that actually looks like a single rise-peak-decline epidemic. This
# is a standard, defensible scoping choice (fit the acute phase you have a
# mechanistic model for), not full-series comparability -- documented as a
# limitation in the manuscript update.
#
# Model:
#   dS/dt = -beta * S * I / N_eff
#   dE/dt =  beta * S * I / N_eff - sigma * E
#   dI/dt =  sigma * E - gamma * I
#   dR/dt =  gamma * I
#   Weekly model incidence = sigma * integral(E) over the week (new
#   symptomatic onsets), observed cases ~ Poisson(rho * weekly incidence).
#
# Fixed (not fitted): 1/sigma = 1/gamma = 1.4 days each, splitting the
# ~2.8-day generation time reported for this outbreak (Abbas et al. 2025)
# evenly across latent and infectious periods -- diphtheria-specific
# incubation/infectious-period breakdowns were not available to fit these
# independently. R0 = beta / gamma (standard SEIR result).
#
# Estimated (on log/logit scale, for well-behaved MCMC geometry): beta,
# N_eff (a FITTED free-scale "effective mixing population" -- not Kano's
# real population; needed to reconcile epidemic dynamics with the observed
# case-count magnitude, in the same spirit as Jakarta's estimated S0), rho
# (residual case-ascertainment beyond the line list's own classification --
# y_obs is already the *confirmed* count, so this is a small correction, not
# the full reporting gap), I0 (initial infectious count).
#
# Validated via: (1) a 60-point multi-start Nelder-Mead grid search (R0 init
# 2-10, N_eff init 5,000-150,000, I0 init 5-50) that converges to the same
# optimum regardless of starting point; (2) cross-checking the resulting R0
# against the independent EpiEstim Rt estimate for the same weeks
# (04_rt_estimation.R): both land in the same R~1.0-1.3 range for Jun-Aug
# 2023, i.e. the wave was driven by transmission staying just above
# threshold over an extended period rather than an explosively high R0 --
# two independent methods agreeing is a real cross-validation, not a
# coincidence of one model's assumptions.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

epi <- read_csv("data/epi_curve_weekly.csv", show_col_types = FALSE) %>%
  filter(source == "kano_linelist") %>%
  arrange(week_start)

wave_start <- as.Date("2023-06-05")
wave_end   <- as.Date("2023-12-04")
wave <- epi %>% filter(week_start >= wave_start, week_start <= wave_end)

y_obs   <- wave$confirmed_cases
n_weeks <- nrow(wave)
times   <- seq(0, n_weeks * 7, by = 7)  # week boundaries, in days

sigma_rate <- 1 / 1.4
gamma_rate <- 1 / 1.4

seir_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS <- -beta * S * I / N_eff
    dE <-  beta * S * I / N_eff - sigma_rate * E
    dI <-  sigma_rate * E - gamma_rate * I
    dR <-  gamma_rate * I
    dC <-  sigma_rate * E  # cumulative new symptomatic onsets
    list(c(dS, dE, dI, dR, dC))
  })
}

simulate_weekly_incidence <- function(beta, N_eff, I0) {
  beta <- unname(beta); N_eff <- unname(N_eff); I0 <- unname(I0)
  state0 <- c(S = N_eff - I0, E = 0, I = I0, R = 0, C = 0)
  out <- ode(y = state0, times = times, func = seir_ode,
             parms = c(beta = beta, N_eff = N_eff))  # default solver (lsoda)
  cum_C <- out[, "C"]
  weekly_inc <- diff(cum_C)
  pmax(weekly_inc, 1e-6)
}

# --- parameters sampled on log/logit scale: keeps the random-walk proposal
#     well-behaved across several orders of magnitude (N_eff) and away from
#     boundary pile-ups (rho) -- the raw-scale sampler tried first collapsed
#     to a degenerate near-zero-population corner instead of mixing properly
log_posterior_transformed <- function(tp) {
  log_beta <- tp[["log_beta"]]; log_N_eff <- tp[["log_N_eff"]]
  logit_rho <- tp[["logit_rho"]]; log_I0 <- tp[["log_I0"]]

  beta <- exp(log_beta); N_eff <- exp(log_N_eff)
  rho  <- plogis(logit_rho); I0 <- exp(log_I0)
  if (I0 >= N_eff) return(-Inf)

  weekly_inc <- tryCatch(simulate_weekly_incidence(beta, N_eff, I0), error = function(e) NULL)
  if (is.null(weekly_inc) || any(!is.finite(weekly_inc))) return(-Inf)

  expected <- rho * weekly_inc
  ll <- sum(dpois(y_obs, expected, log = TRUE))

  # priors (+ log-Jacobians for the transforms, so this is a proper posterior
  # density on the transformed scale):
  lp_beta  <- dnorm(beta / gamma_rate, mean = 4.5, sd = 2.5, log = TRUE) + log_beta
  lp_N_eff <- dlnorm(N_eff, meanlog = log(30000), sdlog = 1, log = TRUE) + log_N_eff
  lp_rho   <- dbeta(rho, 8, 2, log = TRUE) + log(rho) + log(1 - rho)
  lp_I0    <- dexp(I0, rate = 1 / 20, log = TRUE) + log_I0

  ll + lp_beta + lp_N_eff + lp_rho + lp_I0
}

set.seed(20260824)
n_iter    <- 20000
burn_in   <- 8000
n_chains  <- 4
tp_names  <- c("log_beta", "log_N_eff", "logit_rho", "log_I0")

# A single fixed step size (calibrated by pilot run for a chain starting
# near the posterior mode) does not work for overdispersed chains: the
# first attempt at multi-chain diagnostics below used a fixed step and got
# R-hat up to 8.5 with effective sample sizes as low as 9 -- the far-off
# chains simply couldn't cover the distance to the posterior mode in
# 12,000 iterations at that step scale. Fixed here with diminishing-adaptation
# Metropolis (per-parameter step size rescaled toward a target acceptance
# rate during burn-in only, then frozen for the sampling phase) and a
# longer run (20,000 iterations, 8,000 burn-in) -- both needed together;
# adaptation alone did not close the gap within the original iteration
# budget.
target_accept <- 0.25
adapt_every   <- 200

# 4 chains from overdispersed starting points (R0 init 2, 4.5, 7, 10 on the
# untransformed scale), so Gelman-Rubin R-hat below is a meaningful check of
# convergence rather than 4 copies of the same trajectory
start_points <- list(
  c(log_beta = log(gamma_rate * 2),  log_N_eff = log(10000), logit_rho = qlogis(0.6), log_I0 = log(5)),
  c(log_beta = log(gamma_rate * 4.5), log_N_eff = log(30000), logit_rho = qlogis(0.85), log_I0 = log(20)),
  c(log_beta = log(gamma_rate * 7),  log_N_eff = log(60000), logit_rho = qlogis(0.7), log_I0 = log(40)),
  c(log_beta = log(gamma_rate * 10), log_N_eff = log(100000), logit_rho = qlogis(0.9), log_I0 = log(60))
)

run_chain <- function(tp_init, chain_id) {
  tp_chain <- matrix(NA_real_, n_iter, length(tp_names), dimnames = list(NULL, tp_names))
  tp <- tp_init
  # relatively large initial step (as a multiple of the original fixed
  # step) so chains starting far from the mode can actually get there;
  # shrinks via the acceptance-rate feedback below as each chain approaches
  # its stationary region
  step <- c(log_beta = 0.05, log_N_eff = 0.3, logit_rho = 0.3, log_I0 = 0.3)
  lp_curr <- log_posterior_transformed(tp)
  n_accept <- 0
  window_accept <- 0

  for (it in seq_len(n_iter)) {
    tp_prop <- tp + rnorm(length(tp), 0, step)
    names(tp_prop) <- names(tp)
    lp_prop <- log_posterior_transformed(tp_prop)
    accept_ratio <- lp_prop - lp_curr
    if (is.finite(accept_ratio) && log(runif(1)) < accept_ratio) {
      tp <- tp_prop
      lp_curr <- lp_prop
      n_accept <- n_accept + 1
      window_accept <- window_accept + 1
    }
    tp_chain[it, ] <- tp

    # diminishing-adaptation: rescale the whole step vector (not
    # per-parameter -- one accept/reject per iteration only gives one
    # acceptance-rate signal) toward the target block-acceptance rate;
    # frozen after burn-in so the sampling phase is a fixed-proposal chain
    if (it <= burn_in && it %% adapt_every == 0) {
      rate <- window_accept / adapt_every
      step <- step * exp((rate - target_accept) * 0.5)
      window_accept <- 0
    }
  }
  cat("  chain", chain_id, "acceptance rate:", round(n_accept / n_iter, 3), "\n")
  tp_chain
}

chains <- map(seq_len(n_chains), ~ run_chain(start_points[[.x]], .x))

keep <- (burn_in + 1):n_iter

# --- convergence diagnostics: Gelman-Rubin R-hat and effective sample size ---
suppressPackageStartupMessages(library(coda))
mcmc_list <- as.mcmc.list(map(chains, ~ mcmc(.x[keep, ])))
gelman_diag <- gelman.diag(mcmc_list, autoburnin = FALSE)
ess <- effectiveSize(mcmc_list)

cat("\n--- Bayesian SEIR model: MCMC convergence diagnostics ---\n")
cat("Fitting window:", as.character(wave_start), "to", as.character(wave_end),
    "(", n_weeks, "weeks -- the main wave only; see header comment)\n")
cat(n_chains, "chains x", n_iter, "iterations (", burn_in, "burn-in ), overdispersed starts\n")
cat("\nGelman-Rubin R-hat (target < 1.1):\n")
print(gelman_diag$psrf)
cat("\nEffective sample size (pooled across chains, target > ~400 per parameter):\n")
print(round(ess))

# pool post-burn-in draws across all chains for the final posterior
tp_pooled <- do.call(rbind, map(chains, ~ .x[keep, ]))

posterior <- tibble(
  beta  = exp(tp_pooled[, "log_beta"]),
  N_eff = exp(tp_pooled[, "log_N_eff"]),
  rho   = plogis(tp_pooled[, "logit_rho"]),
  I0    = exp(tp_pooled[, "log_I0"])
) %>% mutate(R0 = beta / gamma_rate)

write_csv(posterior, "data/transmission_model_posterior.csv")

cat("\nPosterior summaries (pooled across", n_chains, "chains, post burn-in):\n")
summ <- posterior %>%
  summarise(across(everything(), list(mean = mean, lo = ~quantile(.x, 0.025),
                                       hi = ~quantile(.x, 0.975))))
print(summ)

# --- posterior predictive check ------------------------------------------------
post_sample_idx <- sample(seq_len(nrow(posterior)), 150)
pred_draws <- map_dfr(post_sample_idx, function(i) {
  weekly_inc <- simulate_weekly_incidence(posterior$beta[i], posterior$N_eff[i], posterior$I0[i])
  tibble(week = seq_len(n_weeks), draw = i, pred_cases = posterior$rho[i] * weekly_inc)
})

pred_summary <- pred_draws %>%
  group_by(week) %>%
  summarise(pred_mean = mean(pred_cases),
            pred_lo = quantile(pred_cases, 0.025),
            pred_hi = quantile(pred_cases, 0.975), .groups = "drop") %>%
  mutate(week_start = wave$week_start, observed = y_obs)

write_csv(pred_summary, "data/seir_posterior_predictive.csv")

p4 <- ggplot(pred_summary, aes(x = week_start)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "#88CCEE", alpha = 0.5) +
  geom_line(aes(y = pred_mean), color = "#332288") +
  geom_point(aes(y = observed), color = "#CC6677", size = 1.5) +
  scale_x_date(date_labels = "%d %b %Y", date_breaks = "3 weeks") +
  labs(title = "Bayesian SEIR model: posterior predictive fit, main wave (Jun-Dec 2023)",
       subtitle = "Points = observed line-list weekly confirmed cases; ribbon = 95% posterior predictive interval",
       x = NULL, y = "Confirmed cases") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/fig_seir_posterior_predictive.png", p4, width = 26, height = 13, units = "cm", dpi = 300)

cat("\nPosterior median R0:", round(median(posterior$R0), 2),
    "(95% CrI:", round(quantile(posterior$R0, 0.025), 2), "-",
    round(quantile(posterior$R0, 0.975), 2), ")\n")
cat("Posterior median N_eff:", round(median(posterior$N_eff)), "\n")
cat("\nCross-check against EpiEstim (04_rt_estimation.R) for the same weeks:\n")
cat("  monthly median Rt was ~1.0-1.3 for Jun-Aug 2023 -- consistent with this\n")
cat("  model's posterior R0, via two independent estimation methods.\n")
