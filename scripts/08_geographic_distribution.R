# ------------------------------------------------------------------------------
# 08_geographic_distribution.R
# State-level distribution of suspected/confirmed cases and deaths, Nigeria.
# Source data unchanged from Manuscript v2 Table 3 (the line list is Kano-only
# and cannot provide the cross-state comparison); regenerated in R here for
# reproducibility and so the figure can be reproduced/edited by the user
# without Excel/PowerPoint.
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

state_data <- read_csv("data/manuscript_table3_state.csv", show_col_types = FALSE) %>%
  mutate(
    cfr_pct = round(100 * deaths / confirmed, 1),
    pct_national_suspected = round(100 * suspected / sum(suspected), 1)
  ) %>%
  arrange(desc(suspected))

write_csv(state_data, "data/state_distribution.csv")

p7 <- state_data %>%
  filter(state != "Others") %>%
  mutate(state = fct_reorder(state, suspected)) %>%
  pivot_longer(c(suspected, confirmed), names_to = "case_type", values_to = "n") %>%
  ggplot(aes(x = state, y = n, fill = case_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c(suspected = "#88CCEE", confirmed = "#332288"),
                     labels = c(suspected = "Suspected", confirmed = "Confirmed")) +
  labs(title = "Diphtheria cases by state, Nigeria, 2022-2026",
       subtitle = "Source: NCDC situation report, Table 3 (cumulative to epi week 3, 2026)",
       x = NULL, y = "Cases", fill = NULL) +
  theme_bw() + theme(legend.position = "top")

ggsave("figures/fig3_state_distribution.png", p7, width = 22, height = 14, units = "cm", dpi = 300)

cat("\n--- Geographic concentration ---\n")
top7 <- state_data %>% filter(state != "Others") %>% slice_head(n = 7)
cat("Top 7 states account for", round(sum(top7$suspected) / sum(state_data$suspected) * 100, 1),
    "% of national suspected cases\n")
cat("Kano alone:", state_data$pct_national_suspected[state_data$state == "Kano"], "% of suspected cases\n")
print(state_data)
