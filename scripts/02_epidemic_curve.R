# ------------------------------------------------------------------------------
# 02_epidemic_curve.R
# Build the real (case-level) weekly epidemic curve from the Kano line list,
# 2022-05-14 to 2024-05-30, and extend it with the manuscript's sparse
# situation-report points for the 2024-08 -> 2026-01 gap the line list does
# not cover. This replaces Manuscript v2's Table 1 / Figure 1, which relied
# on only 11 weekly points for the *entire* outbreak.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

ll <- read_csv("data/linelist_clean.csv", show_col_types = FALSE)

# --- weekly incidence directly from the line list (confirmed cases only) -----
weekly_linelist <- ll %>%
  filter(confirmed) %>%
  mutate(week_start = floor_date(doo, unit = "week", week_start = 1)) %>%
  count(week_start, name = "confirmed_cases") %>%
  mutate(source = "kano_linelist") %>%
  arrange(week_start)

weekly_deaths_linelist <- ll %>%
  filter(confirmed, outcome == "dead") %>%
  mutate(week_start = floor_date(doo, unit = "week", week_start = 1)) %>%
  count(week_start, name = "deaths") %>%
  arrange(week_start)

weekly_linelist <- weekly_linelist %>%
  left_join(weekly_deaths_linelist, by = "week_start") %>%
  mutate(deaths = replace_na(deaths, 0))

# --- gap-filling manuscript points, restricted to weeks after the line list's
#     last observed week (2024-05-30) so we do not double count -----------------
manuscript_pts <- read_csv("data/manuscript_table1_weekly.csv", show_col_types = FALSE) %>%
  filter(week_start > max(weekly_linelist$week_start)) %>%
  transmute(week_start, confirmed_cases, deaths = NA_real_, source)

epi_curve <- bind_rows(weekly_linelist, manuscript_pts) %>%
  arrange(week_start)

write_csv(epi_curve, "data/epi_curve_weekly.csv")

# --- figure --------------------------------------------------------------------
gap_start <- max(weekly_linelist$week_start)
gap_end   <- min(manuscript_pts$week_start)

p1 <- ggplot(epi_curve, aes(x = week_start, y = confirmed_cases, fill = source)) +
  geom_col(width = 6) +
  annotate("rect", xmin = gap_start, xmax = gap_end, ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "grey30") +
  scale_fill_manual(
    values = c(kano_linelist = "#117733", manuscript_v2_table1 = "#CC6677"),
    labels = c(kano_linelist = "Kano line list (case-level)",
               manuscript_v2_table1 = "NCDC situation report points (Table 1)")
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "4 months") +
  labs(
    title = "Weekly confirmed diphtheria cases, Nigeria, 2022-2026",
    subtitle = "Shaded band: no case-level or sitrep data available (2024-06 to 2023-11 gap in reporting)",
    x = NULL, y = "Confirmed cases", fill = NULL
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

ggsave("figures/fig1_epidemic_curve.png", p1, width = 24, height = 12, units = "cm", dpi = 300)

cat("\n--- Epidemic curve summary ---\n")
cat("Line-list weeks:      ", nrow(weekly_linelist), "(", as.character(min(weekly_linelist$week_start)),
    "to", as.character(max(weekly_linelist$week_start)), ")\n")
cat("Manuscript gap-fill weeks:", nrow(manuscript_pts), "\n")
cat("Total confirmed, line-list period:", sum(weekly_linelist$confirmed_cases), "\n")
cat("Total deaths, line-list period:   ", sum(weekly_linelist$deaths), "\n")
cat("Peak week (line list):", as.character(weekly_linelist$week_start[which.max(weekly_linelist$confirmed_cases)]),
    "with", max(weekly_linelist$confirmed_cases), "confirmed cases\n")