# ------------------------------------------------------------------------------
# 07_policy_impact.R
# Counterfactual "cases averted" for policy interventions, replacing
# Manuscript v2's Table 4 / Figure 4 bar chart -- which the manuscript itself
# labelled as illustrative/placeholder numbers "for demonstration only".
#
# Methodologically parallels the Jakarta paper's calculate_diff() scenario
# comparisons (no-ORI, no-contact-tracing, delayed-ORI vs. observed), but run
# here across the full MCMC posterior from 05, so each "cases averted"
# estimate carries a credible interval rather than being a single number from
# one parameter set.
#
# Of the manuscript's 5 original interventions (Table 4), only two act on
# *transmission* and can be quantified this way:
#   - Reactive vaccination (fast vs. delayed, moderate vs. intensive coverage)
#   - Contact tracing and prophylaxis (modelled as a reduction in beta from
#     faster case isolation/post-exposure prophylaxis)
# The other three -- routine immunisation intensification (a years-long,
# pre-outbreak intervention, not a within-wave counterfactual), secure
# DAT/antibiotic supply (acts on case fatality, not incidence -- already
# addressed by the age-CFR analysis in 03_age_cfr.R), and community
# engagement (acts diffusely on both vaccine uptake and care-seeking, without
# a clean single model parameter) -- are not mechanistically quantifiable
# from this model and are reported qualitatively only, as the manuscript
# itself does for several of them.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

posterior <- read_csv("data/transmission_model_posterior.csv", show_col_types = FALSE)

sigma_rate <- 1 / 1.4
gamma_rate <- 1 / 1.4
sim_days   <- 200

seir_intervention_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    beta_t <- if (t >= delay_days) beta * (1 - beta_reduction) else beta
    v_t    <- if (t >= delay_days) vacc_cov_rate else 0
    dS <- -beta_t * S * I / N_eff - v_t * S
    dE <-  beta_t * S * I / N_eff - sigma_rate * E
    dI <-  sigma_rate * E - gamma_rate * I
    dR <-  gamma_rate * I + v_t * S
    dC <-  sigma_rate * E
    list(c(dS, dE, dI, dR, dC))
  })
}

scenarios <- tribble(
  ~scenario,                                  ~vacc_cov_rate, ~beta_reduction, ~delay_days,
  "Baseline (no intervention)",                0,              0,               0,
  "Reactive vaccination, immediate, moderate", 0.002,          0,               0,
  "Reactive vaccination, immediate, intensive",0.01,           0,               0,
  "Reactive vaccination, delayed 4 weeks",     0.01,           0,               28,
  "Contact tracing & prophylaxis (-25% beta)", 0,              0.25,            0
)

run_scenario <- function(beta, N_eff, I0, vacc_cov_rate, beta_reduction, delay_days) {
  beta <- unname(beta); N_eff <- unname(N_eff); I0 <- unname(I0)
  state0 <- c(S = N_eff - I0, E = 0, I = I0, R = 0, C = 0)
  out <- ode(y = state0, times = seq(0, sim_days, by = 1), func = seir_intervention_ode,
             parms = c(beta = beta, N_eff = N_eff, vacc_cov_rate = vacc_cov_rate,
                       beta_reduction = beta_reduction, delay_days = delay_days))
  out[nrow(out), "C"]  # cumulative symptomatic onsets at end of simulation
}

set.seed(20260824)
post_sample <- posterior %>% slice_sample(n = 300) %>% mutate(draw_id = row_number())

cat("Running", nrow(post_sample), "posterior draws x", nrow(scenarios), "scenarios...\n")

results <- post_sample %>%
  cross_join(scenarios) %>%
  pmap_dfr(function(beta, N_eff, I0, rho, draw_id, scenario, vacc_cov_rate, beta_reduction, delay_days, ...) {
    cum_C <- run_scenario(beta, N_eff, I0, vacc_cov_rate, beta_reduction, delay_days)
    tibble(draw_id = draw_id, scenario = scenario, rho = rho, cum_cases = rho * cum_C)
  })

baseline_cases <- results %>% filter(scenario == "Baseline (no intervention)") %>%
  select(draw_id, baseline_cases = cum_cases)

impact <- results %>%
  left_join(baseline_cases, by = "draw_id") %>%
  mutate(cases_averted = baseline_cases - cum_cases,
         pct_averted = 100 * cases_averted / baseline_cases)

impact_summary <- impact %>%
  group_by(scenario) %>%
  summarise(
    cum_cases_mean = mean(cum_cases), cum_cases_lo = quantile(cum_cases, 0.025), cum_cases_hi = quantile(cum_cases, 0.975),
    cases_averted_mean = mean(cases_averted), cases_averted_lo = quantile(cases_averted, 0.025), cases_averted_hi = quantile(cases_averted, 0.975),
    pct_averted_mean = mean(pct_averted), pct_averted_lo = quantile(pct_averted, 0.025), pct_averted_hi = quantile(pct_averted, 0.975),
    .groups = "drop"
  ) %>%
  mutate(scenario = factor(scenario, levels = scenarios$scenario)) %>%
  arrange(scenario)

write_csv(impact_summary, "data/policy_impact_summary.csv")

cat("\n--- Policy impact: cases averted vs. baseline (model-scale population) ---\n")
print(impact_summary %>% select(scenario, cum_cases_mean, cases_averted_mean, pct_averted_mean))

plot_df <- impact_summary %>% filter(scenario != "Baseline (no intervention)")

p6 <- ggplot(plot_df, aes(x = reorder(scenario, pct_averted_mean), y = pct_averted_mean)) +
  geom_col(fill = "#117733") +
  geom_errorbar(aes(ymin = pct_averted_lo, ymax = pct_averted_hi), width = 0.2) +
  coord_flip() +
  labs(title = "Modelled transmission-reduction interventions: cases averted vs. no intervention",
       subtitle = "Bayesian SEIR posterior (main wave fit); error bars = 95% credible interval.\nRoutine immunisation, DAT/antibiotic supply, and community engagement are not\nmodelled here (see script header) and remain qualitative recommendations.",
       x = NULL, y = "Cases averted (%, relative to no intervention)") +
  theme_bw()

ggsave("figures/fig4_policy_impact.png", p6, width = 24, height = 12, units = "cm", dpi = 300)

cat("\nNote: these percentages are relative to this fitted model's own no-intervention\n")
cat("baseline, not literal historical case counts, and cover only the 27-week main-wave\n")
cat("fitting window (05_bayesian_transmission_model.R) -- see that script's header for\n")
cat("why the model is not fit to the full multi-wave series.\n")
