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

# --- bubble map: state centroids sized by confirmed cases -------------------
suppressPackageStartupMessages({
  library(sf)
  library(rnaturalearth)
})

nigeria_states <- ne_states(country = "Nigeria", returnclass = "sf")

# name reconciliation between Table 3's labels and Natural Earth's
name_map <- c(
  "Kano" = "Kano", "Yobe" = "Yobe", "Borno" = "Borno", "Katsina" = "Katsina",
  "Bauchi" = "Bauchi", "Kaduna" = "Kaduna", "Jigawa" = "Jigawa",
  "Sokoto" = "Sokoto", "Gombe" = "Gombe", "Lagos" = "Lagos", "Plateau" = "Plateau"
)

centroids <- nigeria_states %>%
  filter(name %in% name_map) %>%
  st_centroid() %>%
  transmute(state = name, geometry) %>%
  st_as_sf()
coords <- st_coordinates(centroids)
centroids <- centroids %>% mutate(lon = coords[, 1], lat = coords[, 2]) %>% st_drop_geometry()

map_data <- state_data %>%
  filter(state != "Others") %>%
  inner_join(centroids, by = "state")

if (nrow(map_data) != 11) {
  warning("Expected 11 named states to match Natural Earth boundaries, got ", nrow(map_data))
}

p8 <- ggplot() +
  geom_sf(data = nigeria_states, fill = "grey90", color = "grey50", linewidth = 0.3) +
  geom_point(data = map_data, aes(x = lon, y = lat, size = confirmed),
             color = "#CC3333", alpha = 0.75) +
  scale_size_area(max_size = 26, breaks = c(0, 5000, 10000, 15000, 20000, 25000),
                   labels = scales::comma, name = "Confirmed\ncases") +
  labs(title = "Geographic distribution of confirmed diphtheria cases, Nigeria, 2022-2026",
       subtitle = "Bubble size represents cumulative confirmed cases by state (Table 3)",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave("figures/fig3_map.png", p8, width = 22, height = 18, units = "cm", dpi = 300)

cat("\n--- Geographic concentration ---\n")
top7 <- state_data %>% filter(state != "Others") %>% slice_head(n = 7)
cat("Top 7 states account for", round(sum(top7$suspected) / sum(state_data$suspected) * 100, 1),
    "% of national suspected cases\n")
cat("Kano alone:", state_data$pct_national_suspected[state_data$state == "Kano"], "% of suspected cases\n")
print(state_data)
