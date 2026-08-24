# ------------------------------------------------------------------------------
# 01_clean_linelist.R
# Read and clean the Kano diphtheria line list (case-level surveillance data).
#
# The raw workbook's "Final Classification" / "Epiweek" / "Age_Months" /
# "Age Group" columns are Excel formulas that reference other workbook copies
# and read as errors once exported; we recompute them directly in R using the
# same 3-pathway NCDC case-definition logic described in the manuscript's
# Methods (lab-confirmed / epi-linked / clinically compatible / discarded /
# pending / unknown).
#
# Output: data/linelist_clean.csv (one row per case, Kano State only,
#         2022-05-14 to 2024-05-30).
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")

raw <- read_excel("data/kano_linelist_raw.xlsx", sheet = "Line List",
                   col_names = TRUE, guess_max = 25000)

# The sheet has ~230 trailing blank columns left over from Excel formatting;
# keep only the named, populated columns.
raw <- raw[, !is.na(names(raw)) & names(raw) != ""]

ll <- raw %>%
  transmute(
    row_id            = `...1`,
    epid_number       = `EPID Number`,
    state              = `State of residence`,
    lga                = `LGA of residence`,
    ward               = Ward,
    settlement         = Settlement,
    facility           = `Name of health facility`,
    dob                = as.Date(`Date of birth (day/mnth/yr)`),
    age_years          = as.numeric(`Age (yrs)`),
    age_months_raw     = as.numeric(`Age (mnth)`),
    sex                = str_to_lower(`Sex (m/f)`),
    date_seen          = as.Date(`Date seen at HF (day/mnth/yr)`),
    date_onset         = as.Date(`Date of onset (day/mnth/yr)`),
    pseudomembrane     = str_to_lower(`Pseudo-membrane`),
    epid_linkage       = str_to_lower(`Epid linkage`),
    hospitalised       = str_to_lower(Hospitalization),
    outcome_raw        = Outcome,
    status_raw         = Status,
    severity           = Severity,
    sample_collected   = str_to_lower(`Sample collected (yes/no)`),
    culture_result     = str_to_lower(`Culture result for C diphtheriea`),
    toxin_result       = str_to_lower(`Final Results (Toxin Test)`),
    vacc_status_raw    = str_to_lower(`Vaccination status`),
    n_doses            = suppressWarnings(as.numeric(`No of doses`)),
    antibiotics        = `Antibiotics administered`,
    dat_administered   = str_to_lower(`DAT administered`),
    contact_history    = str_to_lower(`History of Contact`)
  )

# --- outcome / vital status: collapse inconsistent capitalisation ------------
ll <- ll %>%
  mutate(
    outcome = case_when(
      str_to_lower(str_trim(outcome_raw)) %in% c("dead", "died") ~ "dead",
      str_to_lower(str_trim(outcome_raw)) == "alive"              ~ "alive",
      TRUE                                                        ~ NA_character_
    )
  )

# --- vaccination status: collapse to 4 categories -----------------------------
ll <- ll %>%
  mutate(
    vacc_status = case_when(
      str_detect(vacc_status_raw, "not")   ~ "not vaccinated",
      str_detect(vacc_status_raw, "part")  ~ "partially vaccinated",
      str_detect(vacc_status_raw, "full")  ~ "fully vaccinated",
      str_detect(vacc_status_raw, "unk")   ~ "unknown",
      TRUE ~ NA_character_
    )
  )

# --- age in completed months (fills gaps the way the workbook formula did) ---
ll <- ll %>%
  mutate(
    age_months = if_else(is.na(age_months_raw) | age_months_raw == 0,
                          age_years * 12, age_months_raw + age_years * 12),
    age_group = cut(
      age_months,
      breaks = c(-Inf, 12, 5*12, 10*12, 15*12, 20*12, 30*12, 40*12, Inf),
      labels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40+"),
      right = FALSE
    )
  )

# --- date of onset (DOO), falling back to date seen where onset is missing,
#     mirroring the workbook's own IF(ISBLANK(...)) logic ---------------------
ll <- ll %>%
  mutate(doo = coalesce(date_onset, date_seen))

# --- Final Classification: NCDC 3-pathway case definition ---------------------
# Reconstructed from the original workbook's nested-IF formula (columns V, AA,
# AJ, AP, AU), which gates on the *clinical* case definition (pseudomembrane
# present, column V) rather than "Epid linkage" (column AA) as its primary
# criterion -- "Epid linkage" is in fact blank/"no" for every single row in
# this file, so that column carries no real information here and the
# "2_Epid linked" branch (which the original formula nests *after* the
# clinically-compatible branch, so it can never actually fire) is effectively
# dead in both the source workbook and this reconstruction.
#
# 1_Lab confirmed:    toxin test positive
# 4_Discarded:        culture or toxin test negative
# 3_Clin compatible:  pseudomembrane present + (culture positive, toxin not
#                      done), OR pseudomembrane present + no sample collected
# 2_Epid linked:      pseudomembrane present + epi-linked + no sample
#                      collected (unreachable here -- epid_linkage is always "no")
# 5_Pending:          sample collected + culture/toxin result pending
# 6_Unknown:          no pseudomembrane, or pseudomembrane unknown
ll <- ll %>%
  mutate(
    final_classification = case_when(
      toxin_result == "positive" ~ "1_Lab confirmed",
      culture_result == "negative" | toxin_result == "negative" ~ "4_Discarded",
      pseudomembrane == "yes" & culture_result == "positive" & toxin_result == "not done" ~ "3_Clin compatible",
      pseudomembrane == "yes" & sample_collected == "no" ~ "3_Clin compatible",
      pseudomembrane == "yes" & epid_linkage == "yes" & sample_collected == "no" ~ "2_Epid linked",
      sample_collected == "yes" & (culture_result == "pending" | toxin_result == "pending") ~ "5_Pending",
      pseudomembrane %in% c("no", "unknown") ~ "6_Unknown",
      TRUE ~ NA_character_
    ),
    confirmed = final_classification %in% c("1_Lab confirmed", "2_Epid linked", "3_Clin compatible")
  )

# --- epi week / year, from date of onset --------------------------------------
ll <- ll %>%
  mutate(
    epiweek = epiweek(doo),
    epiyear = epiyear(doo)
  )

# --- drop rows with no usable onset/seen date and no case number (blank pad rows)
ll <- ll %>% filter(!is.na(epid_number) | !is.na(doo))

# dob is >99% missing and unused downstream (age comes from age_years /
# age_months); dropping it avoids readr's column-type guesser mis-inferring
# logical from a mostly-NA sample and then choking on the rare real date.
ll <- ll %>% select(-dob)

write_csv(ll, "data/linelist_clean.csv")

# --- console summary for sanity-checking against the manuscript ---------------
cat("\n--- Kano line list: cleaning summary ---\n")
cat("Total rows:              ", nrow(ll), "\n")
cat("Rows with onset date:    ", sum(!is.na(ll$date_onset)), "\n")
cat("Date range (DOO):        ", as.character(min(ll$doo, na.rm = TRUE)), "to",
    as.character(max(ll$doo, na.rm = TRUE)), "\n")
cat("Confirmed cases:         ", sum(ll$confirmed, na.rm = TRUE), "\n")
cat("Deaths among confirmed:  ", sum(ll$confirmed & ll$outcome == "dead", na.rm = TRUE), "\n")
cat("\nFinal classification breakdown:\n")
print(count(ll, final_classification, sort = TRUE))
cat("\nCompare to manuscript's reported Kano row (Table 3, cumulative to wk3 2026):\n")
cat("  31,171 suspected / 24,687 confirmed / 1,246 deaths\n")
cat("  (line list covers 2022-05 to 2024-05 only, i.e. the early/main wave --\n")
cat("   expect this to be a substantial share of, not equal to, those totals)\n")