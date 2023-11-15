# accepts ptx_model_data which contains the data used for propensity score
# matching in match.Rmd. If you want the post-matching data, use
# `filter()` or `semi_join()` first to restrict to the rows included in the
# matched data
cohens_d_no_ci <- function(data, variable, by, conf.level, ...) {
  d <- cohens_d(data[[variable]] ~ data[[by]], ci = conf.level) |>
    as_tibble() |>
    select(
      estimate = Cohens_d
    )
  return(d)
}

cohens_w <- function(data, variable, by, conf.level, ...) {
  counts <- data |>
    group_by(!!sym(variable), !!sym(by)) |>
    count() |>
    ungroup() |>
    pivot_wider(
      names_from = !!sym(variable),
      values_from = "n"
    )

  xsq_test <- chisq.test(as.matrix(counts[, -1]))
  w <- effectsize(xsq_test, type = "cohens_w", ci = NA) |>
    as_tibble() |>
    rename(estimate = Cohens_w)

  return(w)
}

describe_cohort <- function(ptx_model_data) {
  demo_table <- ptx_model_data %>%
    mutate(start_year = as.numeric(as.character(start_year))) %>%
    select(drug, age, start_year) %>%
    tbl_summary(by = "drug",
                label = list(
                  age ~ "Age at Index Date",
                  start_year ~ "Medication Index Date"
                )
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                    all_continuous() ~ function(x) style_number(x, digits = 2),
                    all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  duration_table <- ptx_model_data %>%
    mutate(lookback_time = index_date - start_date) %>%
    select(drug, lookback_time, follow_up_time) %>%
    mutate(
      lookback_time = lookback_time / 365,
      follow_up_time = follow_up_time / 365
    ) %>%
    tbl_summary(
      by = "drug",
      label = list(
        lookback_time ~ "Lookback Time (Years)",
        follow_up_time ~ "Follow Up Time (Years)"
      )
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                    all_continuous() ~ function(x) style_number(x, digits = 2),
                    all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  rates_table <- ptx_model_data %>%
    select(drug, in_rate, out_rate, mean_ndx, out_dx_rate) %>%
    mutate(mean_ndx_out_rate = mean_ndx * out_rate) %>%
    tbl_summary(
      by = "drug",
      label = list(
        in_rate ~ "Inpatient Events Per Year",
        out_rate ~ "Outpatient Events Per Year",
        mean_ndx ~ "Mean Number of DX Codes Per Outpatient Event",
        out_dx_rate ~ "Outpatient DX Code Incidence Per Year"
      )
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                    all_continuous() ~ function(x) style_number(x, digits = 2),
                    all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  elix_table <- ptx_model_data %>%
    select(drug,
      cm_elix_Alcohol, cm_elix_Anemia, cm_elix_BloodLoss,
      cm_elix_CHF, cm_elix_Coagulopathy, cm_elix_Depression,
      cm_elix_DM, cm_elix_DMcx, cm_elix_Drugs,
      cm_elix_FluidsLytes, cm_elix_HIV, cm_elix_HTN, cm_elix_HTNcx,
      cm_elix_Hypothyroid, cm_elix_Liver,
      cm_elix_Lymphoma, cm_elix_Mets, cm_elix_NeuroOther,
      cm_elix_Obesity, cm_elix_Paralysis, cm_elix_PHTN,
      cm_elix_Psychoses, cm_elix_PUD, cm_elix_Pulmonary,
      cm_elix_PVD, cm_elix_Renal, cm_elix_Rheumatic,
      cm_elix_Tumor, cm_elix_Valvular, cm_elix_WeightLoss
    ) %>%
    tbl_summary(by = "drug",
                label = list(
                  cm_elix_Alcohol ~ "Alcohol Abuse",
                  cm_elix_Anemia ~ "Anemia",
                  cm_elix_BloodLoss ~ "Blood Loss",
                  cm_elix_CHF ~ "Heart Failure",
                  cm_elix_Coagulopathy ~ "Coagulopathy",
                  cm_elix_Depression ~ "Depression",
                  cm_elix_DM ~ "Diabetes Without Complications",
                  cm_elix_DMcx ~ "Diabetes With Complications",
                  cm_elix_Drugs ~ "Drug Abuse",
                  cm_elix_FluidsLytes ~ "Fluid or Electrolyte Disorders",
                  cm_elix_HIV ~ "HIV",
                  cm_elix_HTN ~ "HTN Without Complications",
                  cm_elix_HTNcx ~ "HTN With Complications",
                  cm_elix_Hypothyroid ~ "Hypothyroidism",
                  cm_elix_Liver ~ "Liver Disease",
                  cm_elix_Lymphoma ~ "Lymphoma",
                  cm_elix_Mets ~ "Metastatic Cancer",
                  cm_elix_NeuroOther ~ "Other Neurologic Disorders",
                  cm_elix_Obesity ~ "Obesity",
                  cm_elix_Paralysis ~ "Paralysis",
                  cm_elix_PHTN ~ "Pulmonary Hypertension",
                  cm_elix_Psychoses ~ "Psychoses",
                  cm_elix_PUD ~ "Peptic Ulcer Disease",
                  cm_elix_Pulmonary ~ "COPD",
                  cm_elix_PVD ~ "Peripheral Vascular Disease",
                  cm_elix_Renal ~ "Renal Failure",
                  cm_elix_Rheumatic ~ "Rheumatoid Arthritis",
                  cm_elix_Tumor ~ "Solid Tumor",
                  cm_elix_Valvular ~ "Valvular Disease",
                  cm_elix_WeightLoss ~ "Weight Loss"
              ),
              statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
              digits = list(all_dichotomous() ~ c(0, 1))

    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                     all_continuous() ~ function(x) style_number(x, digits = 2),
                     all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  hypothension_table <- ptx_model_data %>%
    select(drug,
      cm_orthostatic, cm_other_hypo
    ) %>%
    tbl_summary(by = "drug",
                label = list(
                  cm_orthostatic ~ "Orthostatic Hypotension",
                  cm_other_hypo ~ "Other Hypotension"
                ),
                statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
                digits = list(all_dichotomous() ~ c(0, 1))
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                     all_continuous() ~ function(x) style_number(x, digits = 2),
                     all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  psa_table <- ptx_model_data %>%
    mutate(psa_measured = as.numeric(psa_measured),
           cm_abnormal_psa = as.numeric(cm_abnormal_psa)
           ) %>%
    select(drug,
      psa_measured, cm_abnormal_psa
    ) %>%
    tbl_summary(by = "drug",
                label = list(
                  psa_measured ~ "PSA Measurement Taken",
                  cm_abnormal_psa ~ "Diagnosis of Abnormal PSA"
                ),
                statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
                digits = list(all_dichotomous() ~ c(0, 1))
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                     all_continuous() ~ function(x) style_number(x, digits = 2),
                     all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  luts_table <- ptx_model_data %>%
    mutate(cm_slow_stream = as.numeric(cm_slow_stream),
           uroflow = as.numeric(uroflow),
           cystometrogram = as.numeric(cystometrogram)) %>%
    select(drug,
      cm_slow_stream, uroflow, cystometrogram, cm_bph
    ) %>%
    tbl_summary(by = "drug",
                label = list(
                  cm_slow_stream ~ "Slow Urinary Stream Diagnosis",
                  uroflow ~ "Uroflow Study Performed",
                  cystometrogram ~ "Cystometrogram Performed",
                  cm_bph ~ "Diagnosis of BPH"
                ),
                statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
                digits = list(all_dichotomous() ~ c(0, 1))
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                     all_continuous() ~ function(x) style_number(x, digits = 2),
                     all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  other_table <- ptx_model_data %>%
    mutate(cm_anxiety = as.numeric(cm_anxiety),
           cm_ed = as.numeric(cm_ed)) %>%
    select(drug, cm_anxiety, cm_ed) %>%
    tbl_summary(by = "drug",
                label = list(
                  cm_anxiety ~ "Diagnosis of Anxiety",
                  cm_ed ~ "Diagnosis of Erectile Dysfunction"
                ),
                statistic = list(all_dichotomous() ~ "{n} ({p}%)"),
                digits = list(all_dichotomous() ~ c(0, 1))
    ) %>%
    add_difference(test = list(all_continuous() ~ "cohens_d_no_ci",
                               all_categorical() ~ "cohens_w"),
                   estimate_fun = list(
                     all_continuous() ~ function(x) style_number(x, digits = 2),
                     all_categorical() ~ function(x) style_number(x, digits = 2)
                   )
                 )

  final_table <- tbl_stack(
      list(
        demo_table,
        duration_table,
        rates_table,
        elix_table,
        hypothension_table,
        psa_table,
        luts_table,
        other_table),
      group_header = c(
        "Age and Medication Start Time",
        "Duration of Enrollment Time",
        "Rates of Inpatient and Outpatient Encounters",
        "Elixhauser/AHRQ Comorbidity Flags",
        "Hypotension",
        "Prostate Specific Antigen",
        "Bladder Function and Urinary Flow",
        "Other Comorbidities"
      )
  ) %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    as_flex_table()

  return(final_table)
}

if (Sys.info()["sysname"] == "Darwin") {
  shared_dir <- "/Volumes"
} else {
  shared_dir <- "/Shared"
}
model_data <- read_rds(glue::glue("{shared_dir}/lss_jsimmeri/LBD/treated_model_data.rds"))

