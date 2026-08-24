# ------------------------------------------------------------------------------
# 04_rt_estimation.R
# Time-varying effective reproduction number (Rt) via EpiEstim.
#
# Manuscript v2 attempted this on only 11 weekly points and got an
# acknowledged-implausible result (median Rt 6.55, 95% CI 5.85-7.30) which the
# authors discarded. With ~21,700 dated cases from the line list we can build
# a daily incidence series over the main wave (2022-05 to 2024-05) and get a
# materially more stable estimate -- for the period where there's enough
# incidence to support it.
#
# 2022-06 to 2022-12 is too sparse for meaningful 7-day-window Rt estimates:
# a 14-day rolling case count stays in the 0-8 range for that entire stretch
# (checked directly), so EpiEstim swings wildly between near-zero and >20
# with enormous credible intervals whenever a case follows a gap -- an
# artefact of dividing small integers by other small integers, not a real
# transmission signal (EpiEstim itself warns "estimating R too early in the
# epidemic" for this stretch). Case density only becomes adequate for stable
# estimation once the first (smaller) wave gets underway in January 2023
# (14-day count reaches 17, then 60, in the first two weeks of that month).
# We therefore restrict Rt estimation to 2023-01-01 onward and say so
# explicitly, rather than showing a headline "stable, coherent" Rt curve
# that's actually dominated by six months of noise -- the same standard the
# manuscript already holds its own original Rt attempt to.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

ll <- read_csv("data/linelist_clean.csv", show_col_types = FALSE)

rt_start_date <- as.Date("2023-01-01")

daily_incidence <- ll %>%
  filter(confirmed) %>%
  count(doo, name = "I") %>%
  complete(doo = seq(min(doo), max(doo), by = "day"), fill = list(I = 0)) %>%
  arrange(doo) %>%
  filter(doo >= rt_start_date) %>%
  rename(dates = doo)

# Serial interval: manuscript's Methods notes recent Kano-specific evidence of
# a ~2.8-day generation time (Abbas et al. 2025), versus the 15+/-5 day
# literature default it fell back on. We use the shorter, outbreak-specific
# estimate here since it's the better-supported figure for this outbreak.
si_mean <- 2.8
si_sd   <- 1.0  # not reported in Abbas et al.; kept narrow but non-degenerate

config <- make_config(list(mean_si = si_mean, std_si = si_sd,
                            t_start = seq(8, nrow(daily_incidence) - 6),
                            t_end   = seq(14, nrow(daily_incidence))))

rt_est <- estimate_R(incid = daily_incidence, method = "parametric_si", config = config)

rt_out <- rt_est$R %>%
  mutate(
    window_end_date = daily_incidence$dates[t_end],
    window_start_date = daily_incidence$dates[t_start]
  ) %>%
  select(window_start_date, window_end_date,
         Rt_mean = `Mean(R)`, Rt_median = `Median(R)`,
         Rt_lower = `Quantile.0.025(R)`, Rt_upper = `Quantile.0.975(R)`)

write_csv(rt_out, "data/rt_estimates.csv")

p3 <- ggplot(rt_out, aes(x = window_end_date, y = Rt_median)) +
  geom_ribbon(aes(ymin = Rt_lower, ymax = Rt_upper), fill = "#88CCEE", alpha = 0.4) +
  geom_line(color = "#332288", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  labs(
    title = "Time-varying effective reproduction number (Rt), Kano diphtheria outbreak",
    subtitle = str_wrap(paste0(
      "EpiEstim, 7-day sliding windows, serial interval ", si_mean, " ± ", si_sd,
      " days (Abbas et al. 2025); n = ", sum(daily_incidence$I),
      " confirmed cases with onset dates, ", format(rt_start_date, "%b %Y"),
      " onward (earlier data too sparse for stable estimation -- see script header)"
    ), width = 100),
    x = NULL, y = "Rt (median, 95% CrI)"
  ) +
  theme_bw()

ggsave("figures/fig3_rt_estimate.png", p3, width = 24, height = 12, units = "cm", dpi = 300)

cat("\n--- Rt estimation summary ---\n")
cat("Daily series length:", nrow(daily_incidence), "days,",
    as.character(min(daily_incidence$dates)), "to", as.character(max(daily_incidence$dates)), "\n")
cat("Number of 7-day windows estimated:", nrow(rt_out), "\n")
cat("Rt range across windows: ", round(min(rt_out$Rt_median, na.rm = TRUE), 2), "-",
    round(max(rt_out$Rt_median, na.rm = TRUE), 2), "\n")
cat("Rt at peak transmission (max median):\n")
print(rt_out %>% filter(Rt_median == max(Rt_median, na.rm = TRUE)))
cat("\nCompare to Manuscript v2's discarded estimate: median Rt 6.55 (95% CI 5.85-7.30),\n")
cat("based on only 7 usable data points -- flagged by the authors themselves as implausible.\n")
