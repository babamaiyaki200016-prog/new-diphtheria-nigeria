# ------------------------------------------------------------------------------
# 06_vaccination_scenarios.R
# Posterior-predictive vaccination-intervention scenarios, replacing Manuscript
# v2's Figure 5 (a single deterministic SIR run per scenario, absolute
# doses/day against an assumed N = 1,000,000). Here each scenario is
# simulated across the full posterior from 05, so the comparison carries
# credible intervals rather than one point estimate per curve.
#
# Vaccination is expressed as a DAILY COVERAGE RATE (fraction of remaining
# susceptibles vaccinated per day) rather than an absolute daily dose count:
# 05's N_eff is a fitted free-scale parameter (~60,000; see that script's
# header), not Kano's real population, so an absolute "500 doses/day" has no
# stable real-world meaning here -- it would be a huge daily fraction of a
# 60,000-person pool but a tiny one against Kano's actual population. A
# coverage rate is scale-invariant and translates directly to a real-world
# campaign size once applied to an actual target population in the
# manuscript text.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

posterior <- read_csv("data/transmission_model_posterior.csv", show_col_types = FALSE)

sigma_rate <- 1 / 1.4
gamma_rate <- 1 / 1.4
sim_days   <- 200

seir_vacc_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dS <- -beta * S * I / N_eff - vacc_cov_rate * S
    dE <-  beta * S * I / N_eff - sigma_rate * E
    dI <-  sigma_rate * E - gamma_rate * I
    dR <-  gamma_rate * I + vacc_cov_rate * S
    dC <-  sigma_rate * E
    list(c(dS, dE, dI, dR, dC))
  })
}

scenarios <- tibble(
  scenario      = c("No vaccination", "Moderate (0.2%/day coverage)", "Intensive (1%/day coverage)"),
  vacc_cov_rate = c(0, 0.002, 0.01)
)

set.seed(20260824)
post_sample <- posterior %>% slice_sample(n = 150) %>% mutate(draw_id = row_number())

run_scenario <- function(beta, N_eff, I0, vacc_cov_rate) {
  beta <- unname(beta); N_eff <- unname(N_eff); I0 <- unname(I0)
  vacc_cov_rate <- unname(vacc_cov_rate)
  state0 <- c(S = N_eff - I0, E = 0, I = I0, R = 0, C = 0)
  out <- ode(y = state0, times = seq(0, sim_days, by = 1), func = seir_vacc_ode,
             parms = c(beta = beta, N_eff = N_eff, vacc_cov_rate = vacc_cov_rate))
  as_tibble(as.data.frame(out))
}

sim_results <- post_sample %>%
  cross_join(scenarios) %>%
  pmap_dfr(function(beta, N_eff, I0, draw_id, scenario, vacc_cov_rate, ...) {
    run_scenario(beta, N_eff, I0, vacc_cov_rate) %>%
      mutate(draw_id = draw_id, scenario = scenario)
  })

write_csv(sim_results %>% filter(draw_id <= 20), "data/vaccination_scenarios_sample_draws.csv")

scenario_summary <- sim_results %>%
  group_by(scenario, time) %>%
  summarise(I_mean = mean(I), I_lo = quantile(I, 0.025), I_hi = quantile(I, 0.975), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario))

write_csv(scenario_summary, "data/vaccination_scenarios_summary.csv")

peak_summary <- sim_results %>%
  group_by(scenario, draw_id) %>%
  summarise(peak_I = max(I), .groups = "drop") %>%
  group_by(scenario) %>%
  summarise(peak_mean = mean(peak_I), peak_lo = quantile(peak_I, 0.025),
            peak_hi = quantile(peak_I, 0.975), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
  arrange(scenario)

write_csv(peak_summary, "data/vaccination_scenarios_peak_summary.csv")

p5 <- ggplot(scenario_summary, aes(x = time, y = I_mean, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = I_lo, ymax = I_hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("No vaccination" = "#CC6677",
                                 "Moderate (0.2%/day coverage)" = "#DDCC77",
                                 "Intensive (1%/day coverage)" = "#117733")) +
  scale_fill_manual(values = c("No vaccination" = "#CC6677",
                                "Moderate (0.2%/day coverage)" = "#DDCC77",
                                "Intensive (1%/day coverage)" = "#117733")) +
  labs(title = "Posterior-predictive effect of vaccination on epidemic trajectory",
       subtitle = "Bayesian SEIR fitted to Kano line-list data (main wave); ribbons = 95% posterior predictive interval",
       x = "Days from wave onset", y = "Infectious individuals (model-scale population)",
       color = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "top")

ggsave("figures/fig5_vaccination_scenarios.png", p5, width = 24, height = 12, units = "cm", dpi = 300)

cat("\n--- Vaccination scenario peak infection summary (model-scale population, N_eff ~",
    round(median(posterior$N_eff)), ") ---\n")
print(peak_summary)
cat("\nRelative peak reduction vs. no vaccination:\n")
no_vacc_peak <- peak_summary$peak_mean[peak_summary$scenario == "No vaccination"]
peak_summary %>%
  mutate(pct_reduction = 100 * (1 - peak_mean / no_vacc_peak)) %>%
  select(scenario, peak_mean, pct_reduction) %>%
  print()
