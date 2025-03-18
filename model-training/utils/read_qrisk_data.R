
read_qrisk_data <- function(params) {

  # Read baseline data ------------------------------------------------------

  df_baseline <- data.table::fread(paste0(params$data_path, "/preprocessed-data/baseline_characteristics.tsv")) %>% 
    dplyr::mutate(eid = as.character(eid))  %>%
    dplyr::select(
      eid, age, sex, body_mass_index_bmi, current_tobacco_smoking, smoking_status,
      systolic_blood_pressure, diastolic_blood_pressure, ethnic_background,
      blood_pressure_medication, cholesterol_lowering_medication,
      number_of_cigarettes_currently_smoked_daily_current_cigarette_smokers,
      family_history_of_heart_disease,
      medications, standing_height, weight
    ) %>% 
    dplyr::mutate(
      sex_male = ifelse(sex == "Male", 1, 0),
      number_of_cigarattes_smoked = as.numeric(number_of_cigarettes_currently_smoked_daily_current_cigarette_smokers)
    )
  

  # Read clinical biochemistry ----------------------------------------------

  
  df_biochemistry <- data.table::fread(paste0(params$data_path, "/preprocessed-data/blood_biochemistry.tsv")) %>% 
    dplyr::mutate(eid = as.character(eid)) %>% 
    dplyr::select(eid, total_cholesterol = cholesterol, hdl_cholesterol = hdl_cholesterol)
  
  
  
  # Prevalent conditions ----------------------------------------------------
  
  preval_conditions <- 
    tibble::tribble(
      ~icd10, ~name,
      "I48", "atrial_fibrillation",
      c("M05", "M06"), "rheumatoid_arthritis",
      "N18", "chronic_kidney_disease",
      "E10", "type_1_diabetes",
      c("E11", "E14"), "type_2_diabetes",
      "G43", "migraine",
      "M32", "systemic_lupus_erythematosus",
      c("F20", "F29", "F31"), "severe_mental_illness",
      "N48", "erectile_dysfunction"
    )
  
  df_preval <-
    purrr::map2(
      .x = preval_conditions$icd10,
      .y = preval_conditions$name,
      .f = function(icd10_codes, name) {
        df <- purrr::map(
          .x = icd10_codes,
          .f = ~readr::read_tsv(paste0(params$data_path, "/endpoints/first-occurrences/", .x, ".tsv")) %>% 
            dplyr::select(eid, starts_with("PREVAL")) 
        ) %>% 
          purrr::reduce(full_join, by = "eid") %>% 
          dplyr::mutate(
            PREVALENT_ENDPOINT = rowSums(select(., starts_with("PREVAL_")), na.rm = FALSE),
            # Recode to 0 and 1 (at least one ICD10 code for a prevalent disease)
            PREVALENT_ENDPOINT = dplyr::if_else(PREVALENT_ENDPOINT == 0, 0, 1, missing = NA_real_)
          ) %>% 
          dplyr::select(eid, PREVALENT_ENDPOINT) %>% 
          dplyr::rename(!!paste0("prevalent_", name) := PREVALENT_ENDPOINT)
        
        return(df)
      }
    ) %>% 
    purrr::reduce(full_join, by = "eid") %>% 
    dplyr::mutate(eid = as.character(eid))
  

  # Combine -----------------------------------------------------------------

  df_qrisk_data <- df_baseline %>% 
    dplyr::left_join(df_preval, by = "eid") %>% 
    dplyr::left_join(df_biochemistry, by = "eid") %>% 
    # Exclude individuals already on cholesterol lowering treatment
    dplyr::filter(cholesterol_lowering_medication == 0) %>% 
    dplyr::transmute(
      eid,
      age, 
      sex_male,
      height = standing_height,
      weight = weight,
      body_mass_index = body_mass_index_bmi,
      ethnic_background = dplyr::case_when(
        .data$ethnic_background %in% c("White", "British", "Irish", NA, "Any other white background") ~ "White",
        .data$ethnic_background == "Indian" ~ "Indian",
        .data$ethnic_background == "Pakistani" ~ "Pakistani",
        .data$ethnic_background == "Bangladeshi" ~ "Bangladeshi",
        .data$ethnic_background == "Any other Asian background" ~ "Other Asian",
        .data$ethnic_background == "Caribbean" ~ "Caribbean",
        .data$ethnic_background == "African" ~ "African",
        .data$ethnic_background == "Chinese" ~ "Chinese",
        TRUE ~ "Other"
      ),
      systolic_blood_pressure,
      smoking_status = dplyr::case_when(
        .data$smoking_status == "Never" ~ "Never",
        .data$smoking_status == "Previous" ~ "Former",
        .data$smoking_status == "Current" & number_of_cigarattes_smoked< 10 ~ "Light smoker",
        .data$smoking_status == "Current" & between(number_of_cigarattes_smoked, 10, 19) ~ "Moderate smoker",
        .data$smoking_status == "Current" & number_of_cigarattes_smoked >= 20 ~ "Heavy smoker",
        # If number of cigarettes smoked daily is missing, set to moderate smoker (median)
        .data$smoking_status == "Current" & is.na(number_of_cigarattes_smoked) ~ "Moderate smoker",
      ),
      total_cholesterol_to_hdl_cholesterol = total_cholesterol/hdl_cholesterol,
      family_history_of_heart_disease,
      type_1_diabetes = prevalent_type_1_diabetes,
      type_2_diabetes = prevalent_type_2_diabetes,
      systemic_lupus_erythematosus = prevalent_systemic_lupus_erythematosus,
      treated_hypertension = blood_pressure_medication,
      rheumatoid_arthritis = prevalent_rheumatoid_arthritis,
      atrial_fibrillation = prevalent_atrial_fibrillation,
      chronic_kidney_disease = prevalent_chronic_kidney_disease,
      migraine = prevalent_migraine,
      corticosteroid_use = ifelse(
        str_detect(medications, "prednisolone|betamethasone|cortisone|depo-medrone|dexamethasone|deflazacort|efcortesol|hydrocortisone|methylprednisolone|triamcinolone"), 1, 0
      ),
      atypical_antipsychotic_use = ifelse(
        str_detect(medications, "misulpride|aripiprazole|clozapine|lurasidone|olanzapine|paliperidone|quetiapine|risperidone|sertindole|zotepine"), 1, 0
      ),
      severe_mental_illness = prevalent_severe_mental_illness,
      erectile_dysfunction = ifelse(
        sex_male == 1 & (str_detect(medications, "alprostadil|phosphodiesterase|papaverine|phentolamine") | prevalent_erectile_dysfunction == 1), 1, 0
      )
    ) 
  

    
  # Log-transform log-normally distributed values ---------------------------
  
  df_qrisk_data_preprocessed <- df_qrisk_data %>% 
    dplyr::mutate_at(
      .vars = vars(total_cholesterol_to_hdl_cholesterol), 
      .funs = ~log1p(.)
    ) %>% 
    dplyr::mutate_at(
      .vars = vars(body_mass_index, systolic_blood_pressure, total_cholesterol_to_hdl_cholesterol, height, weight), 
      .funs = ~winsorize(., sd_factor = 5)
    ) %>%
    # Transformed variables included in QRISK3
    dplyr::mutate(
      age_1 = (age / 10)^-2,
      dbmi = body_mass_index / 10,
      bmi_1 = dbmi^-2,
      bmi_2 = dbmi^-2 * log(dbmi)
    ) %>% 
    select(-dbmi)
  
  
  # Add interaction terms ---------------------------------------------------
  
  options(na.action='na.pass')
  
  df_final <-
    stats::model.matrix(
      ~., 
      df_qrisk_data_preprocessed %>% select(-eid),
      na.action = na.pass
    ) %>% 
    as_tibble() %>% 
    dplyr::select(-"(Intercept)") %>% 
    dplyr::mutate(eid = df_qrisk_data_preprocessed$eid) 
  
  return(df_final)
  
}



