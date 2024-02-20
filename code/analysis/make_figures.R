library(tidyverse)
library(patchwork)
library(survival)
library(viridis)

# Tables S1-S3 are found in the analysis files

if (Sys.info()["sysname"] == "Darwin") {
  root_dir <- "/Volumes/lss_jsimmeri/LBD/"
} else {
  root_dir <- "/Shared/lss_jsimmeri/LBD/"
}

# Table 1
model_data <- read_rds(glue::glue("{root_dir}/treated_model_data.rds"))
outcomes_2rx <- model_data |>
  filter(two_rx) |>
  select(enrolid, drug, develops_pd, survival_time)

tz_tam_matches <- read_rds(glue::glue("/{root_dir}/matches_tz_tam_2rx.rds")) |>
  mutate(enrolid = as.numeric(enrolid))
tz_5ari_matches <- read_rds(glue::glue("/{root_dir}/matches_tz_5ari_2rx.rds")) |>
  mutate(enrolid = as.numeric(enrolid))
tam_5ari_matches <- read_rds(glue::glue("/{root_dir}/matches_tam_5ari_2rx.rds")) |>
  mutate(enrolid = as.numeric(enrolid))

matched_outcomes_tz_tam <- inner_join(outcomes_2rx, tz_tam_matches, by = c("enrolid"))
matched_outcomes_tz_5ari <- inner_join(outcomes_2rx, tz_5ari_matches, by = c("enrolid"))
matched_outcomes_tam_5ari <- inner_join(outcomes_2rx, tam_5ari_matches, by = c("enrolid"))

outcomes_2rx |>
  group_by(drug) |>
  summarize(
    n = n(),
    cases = sum(develops_pd),
    person_years = sum(survival_time / 365)
  ) |>
  mutate(
    incidence = round(cases / person_years * 10000, 2)
  ) |>
  knitr::kable()

matched_outcomes_tz_tam |>
  group_by(drug) |>
  summarize(
    n = n(),
    cases = sum(develops_pd),
    person_years = sum(survival_time / 365)
  ) |>
  mutate(
    incidence = round(cases / person_years * 10000, 2)
  ) |>
  knitr::kable()

matched_outcomes_tz_5ari |>
  group_by(drug) |>
  summarize(
    n = n(),
    cases = sum(develops_pd),
    person_years = sum(survival_time / 365)
  ) |>
  mutate(
    incidence = round(cases / person_years * 10000, 2)
  ) |>
  knitr::kable()

matched_outcomes_tam_5ari |>
  group_by(drug) |>
  summarize(
    n = n(),
    cases = sum(develops_pd),
    person_years = sum(survival_time / 365)
  ) |>
  mutate(
    incidence = round(cases / person_years * 10000, 2)
  ) |>
  knitr::kable()

# Figure 1
bind_rows(
  survfit(Surv(survival_time, develops_pd) ~ treatment,
          data = outcomes_2rx |>
                   mutate(
                     treatment = drug == "tz/dz/az"
                   ) |>
                   filter(
                     drug %in% c("tz/dz/az", "tamsulosin")
                   )) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=TRUE", "Tz/Dz/Az", "Tamsulosin"),
      comparison = "Tz/Dz/Az vs Tamsulosin",
      matching = "Unadjusted"
    ),
  survfit(Surv(survival_time, develops_pd) ~ treatment,
          data = outcomes_2rx |>
                   mutate(
                     treatment = drug == "tz/dz/az"
                   ) |>
                   filter(
                     drug %in% c("tz/dz/az", "5ari")
                   )) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=TRUE", "Tz/Dz/Az", "5ARI"),
      comparison = "Tz/Dz/Az vs 5ARI",
      matching = "Unadjusted"
    ),
  survfit(Surv(survival_time, develops_pd) ~ treatment,
          data = outcomes_2rx |>
                   mutate(
                     treatment = drug == "tamsulosin"
                   ) |>
                   filter(
                     drug %in% c("tamsulosin", "5ari")
                   )) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=TRUE", "Tamsulosin", "5ARI"),
      comparison = "Tamsulosin vs 5ARI",
      matching = "Unadjusted"
    ),
  survfit(Surv(survival_time, develops_pd) ~ treatment, 
          data = matched_outcomes_tz_tam) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=1", "Tz/Dz/Az", "Tamsulosin"),
      comparison = "Tz/Dz/Az vs Tamsulosin",
      matching = "Propensity Score Matched"
    ),
  survfit(Surv(survival_time, develops_pd) ~ treatment, 
          data = matched_outcomes_tz_5ari) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=1", "Tz/Dz/Az", "5ARI"),
      comparison = "Tz/Dz/Az vs 5ARI",
      matching = "Propensity Score Matched"
    ),
  survfit(Surv(survival_time, develops_pd) ~ treatment, 
          data = matched_outcomes_tam_5ari) |>
    broom::tidy() |>
    mutate(
      drug = ifelse(strata == "treatment=1", "Tamsulosin", "5ARI"),
      comparison = "Tamsulosin vs 5ARI",
      matching = "Propensity Score Matched"
    )
) |>
  mutate(
    matching = forcats::fct_relevel(matching, "Unadjusted"),
    comparison = forcats::fct_relevel(comparison, "Tz/Dz/Az versus Tamsulosin", "Tz/Dz/Az versus 5ARI")
  ) |>
  filter(time <= (365 * 10)) |>
  filter(matching == "Propensity Score Matched") |>
  ggplot(aes(x = time / 365, y = estimate * 100,
             ymin = conf.low * 100, ymax = conf.high * 100, 
             color = drug, fill = drug)) + 
  geom_step() + 
  geom_ribbon(alpha = 0.33, linetype = 0) + 
  facet_grid(cols = vars(forcats::fct_rev(comparison))) +
  theme_minimal() + 
  labs(x = "Years of Follow-Up", y = "Percent Remaining DLB Free",
       color = "", fill = "") + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10), minor_breaks = c(1, 3, 5, 7, 9)) +
  theme(legend.position = c(0.125, 0.25), 
        legend.background = element_rect(fill = NA, linetype = 0)) +
  scale_color_manual(values = c("#00A9E0", "#00AF66", "#FF7000")) +
  scale_fill_manual(values = c("#00A9E0", "#00AF66", "#FF7000"))
ggsave("~/projects/dlb-tz/fig_1.svg", width = 6, height = 3)
