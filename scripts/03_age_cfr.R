# ------------------------------------------------------------------------------
# 03_age_cfr.R
# Age-specific case fatality rates -- computed directly from the line list
# (Manuscript v2 back-calculated these from national proportions applied to
# the cumulative total; here they are simply counted).
#
# Also fits a Bayesian hierarchical (partial-pooling) binomial model: each age
# group's death count is modelled as binomial(cases, p_age), with logit(p_age)
# drawn from a shared normal distribution across age groups. This shrinks the
# sparse strata (e.g. <1, with only a handful of deaths) toward the overall
# mean and is a more defensible use of "Bayesian" here than forcing MCMC onto
# a well-powered stratum. Implemented as a from-scratch Metropolis-within-Gibbs
# sampler in base R (no rstan/rstanarm): this environment's C++ toolchain
# cannot compile Stan models (Xcode SDK 'cmath' header issue) and the
# rstanarm binary could not be downloaded, so a compiled-Bayesian-package
# route was not viable here. Base-R MCMC also means the user can rerun this
# without needing a working Stan toolchain themselves.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

ll <- read_csv("data/linelist_clean.csv", show_col_types = FALSE)

age_group_levels <- c("<1", "1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40+")
ll <- ll %>% mutate(age_group = factor(age_group, levels = age_group_levels))

age_cfr <- ll %>%
  filter(confirmed, !is.na(age_group)) %>%
  group_by(age_group) %>%
  summarise(cases = n(), deaths = sum(outcome == "dead", na.rm = TRUE), .groups = "drop")

# --- Wilson score CIs (direct, not back-calculated) ---------------------------
wilson <- binom.confint(age_cfr$deaths, age_cfr$cases, methods = "wilson")
age_cfr <- age_cfr %>%
  mutate(
    cfr_pct = 100 * deaths / cases,
    ci_low  = 100 * wilson$lower,
    ci_high = 100 * wilson$upper
  )

# ------------------------------------------------------------------------------
# Bayesian hierarchical partial-pooling model
#   deaths_i ~ Binomial(cases_i, p_i)
#   logit(p_i) = mu + sigma * z_i,   z_i ~ Normal(0, 1)
#   mu ~ Normal(logit(overall CFR), 1.5)     [weakly informative]
#   sigma ~ Half-Normal(0, 1)                [weakly informative]
# Fit via random-walk Metropolis-within-Gibbs.
# ------------------------------------------------------------------------------

set.seed(20260824)

n_groups <- nrow(age_cfr)
cases_v  <- age_cfr$cases
deaths_v <- age_cfr$deaths
overall_logit <- qlogis(sum(deaths_v) / sum(cases_v))

log_lik_i <- function(z_i, mu, sigma, cases_i, deaths_i) {
  p_i <- plogis(mu + sigma * z_i)
  dbinom(deaths_i, cases_i, p_i, log = TRUE)
}

log_prior_mu    <- function(mu)    dnorm(mu, overall_logit, 1.5, log = TRUE)
log_prior_sigma <- function(sigma) if (sigma <= 0) -Inf else dnorm(sigma, 0, 1, log = TRUE)
log_prior_z     <- function(z)     dnorm(z, 0, 1, log = TRUE)

n_iter   <- 20000
burn_in  <- 5000
z_chain     <- matrix(NA_real_, n_iter, n_groups)
mu_chain    <- numeric(n_iter)
sigma_chain <- numeric(n_iter)

z     <- rep(0, n_groups)
mu    <- overall_logit
sigma <- 0.3

step_z     <- rep(0.4, n_groups)
step_mu    <- 0.15
step_sigma <- 0.15

for (it in seq_len(n_iter)) {
  # update each z_i (age-group random effect)
  for (i in seq_len(n_groups)) {
    z_prop <- z[i] + rnorm(1, 0, step_z[i])
    lp_curr <- log_lik_i(z[i], mu, sigma, cases_v[i], deaths_v[i]) + log_prior_z(z[i])
    lp_prop <- log_lik_i(z_prop, mu, sigma, cases_v[i], deaths_v[i]) + log_prior_z(z_prop)
    if (log(runif(1)) < (lp_prop - lp_curr)) z[i] <- z_prop
  }
  # update mu
  mu_prop <- mu + rnorm(1, 0, step_mu)
  lp_curr <- sum(sapply(seq_len(n_groups), function(i) log_lik_i(z[i], mu, sigma, cases_v[i], deaths_v[i]))) + log_prior_mu(mu)
  lp_prop <- sum(sapply(seq_len(n_groups), function(i) log_lik_i(z[i], mu_prop, sigma, cases_v[i], deaths_v[i]))) + log_prior_mu(mu_prop)
  if (log(runif(1)) < (lp_prop - lp_curr)) mu <- mu_prop
  # update sigma
  sigma_prop <- sigma + rnorm(1, 0, step_sigma)
  if (sigma_prop > 0) {
    lp_curr <- sum(sapply(seq_len(n_groups), function(i) log_lik_i(z[i], mu, sigma, cases_v[i], deaths_v[i]))) + log_prior_sigma(sigma)
    lp_prop <- sum(sapply(seq_len(n_groups), function(i) log_lik_i(z[i], mu, sigma_prop, cases_v[i], deaths_v[i]))) + log_prior_sigma(sigma_prop)
    if (log(runif(1)) < (lp_prop - lp_curr)) sigma <- sigma_prop
  }

  z_chain[it, ] <- z
  mu_chain[it]  <- mu
  sigma_chain[it] <- sigma
}

keep <- (burn_in + 1):n_iter
# matrix + vector / matrix * vector recycle by row in R (column-major storage,
# vector length == nrow), which is exactly what's needed here: each row is one
# posterior draw, so p_chain[draw, i] = plogis(mu[draw] + sigma[draw] * z[draw, i])
p_chain <- plogis(mu_chain[keep] + sigma_chain[keep] * z_chain[keep, , drop = FALSE])

age_cfr <- age_cfr %>%
  mutate(
    bayes_cfr_pct    = 100 * colMeans(p_chain),
    bayes_ci_low     = 100 * apply(p_chain, 2, quantile, probs = 0.025),
    bayes_ci_high    = 100 * apply(p_chain, 2, quantile, probs = 0.975)
  )

write_csv(age_cfr, "data/age_cfr.csv")

# --- MCMC diagnostics: acceptance-style sanity check (effective spread of mu) --
cat("\n--- Bayesian hierarchical CFR model diagnostics ---\n")
cat("Posterior mean mu (logit scale):   ", round(mean(mu_chain[keep]), 3), "\n")
cat("Posterior mean sigma (pooling SD): ", round(mean(sigma_chain[keep]), 3), "\n")
cat("(sigma near 0 = strong pooling toward a common CFR; sigma large = little pooling)\n")

# --- figure: direct vs. Bayesian-pooled CFR, side by side ---------------------
plot_df <- age_cfr %>%
  select(age_group, cfr_pct, ci_low, ci_high, bayes_cfr_pct, bayes_ci_low, bayes_ci_high) %>%
  pivot_longer(-age_group, names_to = c("estimate", ".value"),
               names_pattern = "(bayes_)?(cfr_pct|ci_low|ci_high)") %>%
  mutate(estimate = if_else(estimate == "bayes_", "Bayesian (partial pooling)", "Direct (Wilson CI)"))

p2 <- ggplot(plot_df, aes(x = age_group, y = cfr_pct, color = estimate)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("Direct (Wilson CI)" = "#CC6677",
                                 "Bayesian (partial pooling)" = "#117733")) +
  labs(title = "Age-specific diphtheria case fatality rate, Kano line list (2022-2024)",
       subtitle = "Direct (independent Wilson CIs) vs. Bayesian hierarchical (partial-pooling) estimates",
       x = "Age group", y = "Case fatality rate (%)", color = NULL) +
  theme_bw() + theme(legend.position = "top")

ggsave("figures/fig2_age_cfr.png", p2, width = 22, height = 12, units = "cm", dpi = 300)

cat("\n--- Age-specific CFR (line list, direct vs Bayesian) ---\n")
print(age_cfr)
