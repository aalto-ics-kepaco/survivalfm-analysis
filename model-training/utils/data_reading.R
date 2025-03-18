create_endpoints <- function(model_list, params) {

  df_endpoints <- purrr::map(
    .x = unique(model_list$endpoint_var),
    .f = function(endpoint_var) {

      message(endpoint_var)

      icd10_codes <- model_list %>%
        dplyr::filter(endpoint_var == !!endpoint_var) %>%
        pull(included_icd10_codes) %>%
        strsplit(., ";", ) %>% 
        unlist() %>% 
        unique()

      df_endpoints <- purrr::map(
        .x = icd10_codes,
        .f = ~if (str_detect(.x, "^C|^D")) {
          file_path <- paste0(params$data_path, "/endpoints/cancer-register/", .x, ".tsv")
          if (file.exists(file_path)) {data.table::fread(file_path)} else {NULL}
        } else {
          file_path <-paste0(params$data_path, "/endpoints/first-occurrences/", .x, ".tsv")
          if (file.exists(file_path)) {data.table::fread(file_path)} else {NULL}
        }
      ) %>%
        purrr::compact() 

      if (length(df_endpoints) > 0) {
        df_endpoints <- df_endpoints %>%
          purrr::reduce(full_join, by = "eid") %>% 
          dplyr::mutate_all(~as.numeric(.)) %>% 
          dplyr::mutate(eid = as.character(eid)) 
      } else {
        return()
      }

      # Read endpoints
      df_endpoint <- df_endpoints %>%
        dplyr::mutate(
          # Any prevalent endpoint excludes incident one
          INCIDENT_ENDPOINT =  rowSums(select(., starts_with("INCIDENT_")), na.rm = FALSE), 
          # Recode to 0 and 1 (at least one ICD10 code for an incident disease)
          INCIDENT_ENDPOINT = dplyr::if_else(INCIDENT_ENDPOINT == 0, 0, 1, missing = NA_real_),
          PREVALENT_ENDPOINT = rowSums(select(., starts_with("PREVAL_")), na.rm = FALSE),
          # Recode to 0 and 1 (at least one ICD10 code for a prevalent disease)
          PREVALENT_ENDPOINT = dplyr::if_else(PREVALENT_ENDPOINT == 0, 0, 1, missing = NA_real_),
          ENDPOINT_AGEDIFF = do.call(pmin, select(., ends_with("_AGEDIFF")))
        ) %>%
        dplyr::select(
          eid,
          INCIDENT_ENDPOINT,
          PREVALENT_ENDPOINT,
          ENDPOINT_AGEDIFF
        )

      # Censoring
      df_endpoint <- df_endpoint %>%
        dplyr::mutate(
          INCIDENT_ENDPOINT = ifelse(INCIDENT_ENDPOINT == 1 &  ENDPOINT_AGEDIFF <= params$survival_time, 1, 0),
          ENDPOINT_AGEDIFF  = ifelse(ENDPOINT_AGEDIFF  >= params$survival_time, params$survival_time, ENDPOINT_AGEDIFF)
        )

      df_final <- df_endpoint %>%
        dplyr::select(
          eid,
          !!paste0("INCIDENT_", endpoint_var) := "INCIDENT_ENDPOINT",
          !!paste0(endpoint_var, "_AGEDIFF") := "ENDPOINT_AGEDIFF"
        )

      return(df_final)
    }
  ) %>%
    purrr::reduce(inner_join, by = "eid")

}


winsorize <- function(x, sd_factor) {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- sd(x, na.rm = TRUE)
  lower_bound <- x_mean - sd_factor * x_sd
  upper_bound <- x_mean + sd_factor * x_sd
  pmin(pmax(x, lower_bound), upper_bound)
}


read_ukb_data <- function(params, model_row) {
  
  # Baseline characteristics ------------------------------------------------
  
  message("Reading baseline data...")
  
  preval_conditions <- 
    tibble::tribble(
      ~icd10, ~name,
      c("E11", "E14"), "type_2_diabetes"
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
            PREVALENT_ENDPOINT = dplyr::if_else(PREVALENT_ENDPOINT == 0, 0, 1, missing = NA_real_)
          ) %>% 
          dplyr::select(eid, PREVALENT_ENDPOINT) %>% 
          dplyr::rename(!!paste0("prevalent_", name) := PREVALENT_ENDPOINT)
        
        return(df)
      }
    ) %>% 
    purrr::reduce(full_join, by = "eid") %>% 
    dplyr::mutate(eid = as.character(eid))
  
  
  df_cholesterols <- data.table::fread(paste0(params$data_path, "/preprocessed-data/blood_biochemistry.tsv")) %>% 
    select(eid, cholesterol, ldl_direct, hdl_cholesterol) %>% 
    dplyr::mutate(eid = as.character(eid)) 
  
  df_cholesterols_preprocessed <- df_cholesterols %>%
    dplyr::mutate_at(vars(-eid), ~log1p(.) %>% winsorize(., sd_factor = 4)) 
  
  df_clinical <- data.table::fread(paste0(params$data_path, "/preprocessed-data/baseline_characteristics.tsv")) %>% 
    dplyr::mutate(eid = as.character(eid)) 
  
  df_clinical_preprocessed <-  dplyr::inner_join(df_clinical, df_cholesterols_preprocessed, by = "eid") %>%
    dplyr::left_join(df_preval, by = "eid") %>% 
    dplyr::select(
      eid, age, sex, body_mass_index_bmi, smoking_status,
      systolic_blood_pressure, diastolic_blood_pressure, ethnic_background,
      blood_pressure_medication,  cholesterol_lowering_medication,
      prevalent_diabetes = prevalent_type_2_diabetes,
      total_cholesterol = cholesterol, 
      hdl_cholesterol = hdl_cholesterol,
      ldl_cholesterol = ldl_direct
    ) 
  
  continuous_vars <- c("body_mass_index_bmi", "systolic_blood_pressure", "diastolic_blood_pressure")
  
  df_standard <- df_clinical_preprocessed %>%
    mutate_at(vars(continuous_vars), ~winsorize(., sd_factor = 4)) %>%
    dplyr::mutate(
      sex_male = ifelse(sex == "Male", 1, 0),
      current_smoker = ifelse(smoking_status %in% c("Current"), 1, 0),
      former_smoker = ifelse(smoking_status %in% c("Previous"), 1, 0),
      ethnic_background_white = ifelse(ethnic_background %in% c("British", "Irish", "Any other white background"), 1, 0),
      ethnic_background_mixed = ifelse(ethnic_background %in% c("White and Black Caribbean", "White and Black African", "White and Asian", "Any other mixed background"), 1, 0),
      ethnic_background_asian = ifelse(ethnic_background %in% c("Indian", "Pakistani", "Bangladeshi", "Any other Asian background"), 1, 0),
      ethnic_background_black = ifelse(ethnic_background %in% c("Caribbean", "African", "Any other Black background"), 1, 0),
      ethnic_background_chinese = ifelse(ethnic_background %in% c("Chinese"), 1, 0)
    ) %>%
    dplyr::select(-sex, -smoking_status, -ethnic_background)
  
  
  
  if (model_row$input_name == "Standard risk factors") {
    df_final <- df_standard 
    return(df_final)
  } 
  
  if (model_row$input_name == "Polygenic risk scores" | model_row$input_name == "All") {
    
    df_prs_preprocessed <- data.table::fread(paste0(params$data_path, "/preprocessed-data/polygenic_risk_scores.tsv")) %>%
      select(eid, starts_with("standard_"))  %>% 
      dplyr::mutate(eid = as.character(eid)) %>%
      mutate_at(vars(-eid), ~winsorize(., sd_factor = 4)) %>%
      dplyr::rename_with(~paste0("polygenic_risk_scores__", .), -eid) %>% 
      tidyr::drop_na()
    
    rm(df_prs)
    
    df_final <- df_prs_preprocessed %>% 
      dplyr::left_join(df_standard, by = "eid") 
    
    if (model_row$input_name == "Polygenic risk scores") {
      return(df_final)
    }
    
  } 
  
  if (model_row$input_name == "Metabolomics biomarkers" | model_row$input_name == "All") {
    
    message("Reading metabolomics data...")
    
    df_metabolomics <- data.table::fread(paste0(params$data_path, "/preprocessed-data/metabolomics.tsv")) %>%
      dplyr::mutate(eid = as.character(eid))
    
    df_metabolomics_meta <- data.table::fread(paste0(params$data_path, "/helper-files/metabolite_encodings.tsv")) %>%
      dplyr::filter(!units %in% c("percent", "ratio"))
    
    df_metabolomics_preprocessed <- df_metabolomics %>%
      dplyr::select(eid, any_of(df_metabolomics_meta$col_name)) %>%
      mutate_at(vars(-eid), ~log1p(.) %>% winsorize(., sd_factor = 4))
    
    df_final <- df_metabolomics_preprocessed %>% 
      dplyr::left_join(df_standard %>% select(-ldl_cholesterol, -hdl_cholesterol, -total_cholesterol), by = "eid") 
    
    if (model_row$input_name == "Metabolomics biomarkers") {
      return(df_final)
    }
  } 
  
  if (model_row$input_name == "Clinical biochemistry and blood counts" | model_row$input_name == "All") {
    
    message("Reading blood count data...")
    
    # Blood count -------------------------------------------------------------
    
    df_blood_counts <- data.table::fread(paste0(params$data_path, "/preprocessed-data/blood_count.tsv")) %>%
      dplyr::mutate(eid = as.character(eid)) %>% 
      # These are mostly 0
      dplyr::select(-contains("nucleated_red_blood_cell")) 
    
    df_blood_counts_preprocessed  <- df_blood_counts %>%
      dplyr::mutate_at(vars(-eid), ~log1p(.) %>% winsorize(., sd_factor = 4)) 
    
    # Blood biochemistry ------------------------------------------------------
    
    message("Reading blood biochemistry data...")
    
    df_full_blood_biochemistry <- data.table::fread(paste0(params$data_path, "/preprocessed-data/blood_biochemistry.tsv")) %>% 
      dplyr::mutate(eid = as.character(eid)) 
    
    df_full_blood_biochemistry_preprocessed <- df_full_blood_biochemistry  %>%
      dplyr::mutate(eid = as.character(eid)) %>% 
      dplyr::mutate_at(vars(-eid), ~log1p(.) %>% winsorize(., sd_factor = 4)) %>%
      dplyr::select(where(~sum(is.na(.)) / length(.) <= 0.2)) 
    
    df_final <- df_full_blood_biochemistry_preprocessed  %>% 
      dplyr::inner_join(df_blood_counts_preprocessed, by = "eid") %>% 
      dplyr::inner_join(df_standard %>% dplyr::select(-hdl_cholesterol, -ldl_cholesterol, -total_cholesterol), by = "eid") 
    
    if (model_row$input_name == "Clinical biochemistry and blood counts") {
      return(df_final)
    }
  }
  
  if (model_row$input_name == "All") {
    df_final <-  df_prs_preprocessed %>% 
      dplyr::inner_join(df_metabolomics_preprocessed,  by = "eid") %>% 
      dplyr::inner_join(
        df_full_blood_biochemistry_preprocessed  %>% 
          # Remove repeated biomarkers measured with different platforms
          dplyr::select(-hdl_cholesterol, -ldl_direct, -cholesterol, -apolipoprotein_b, -glucose, -creatinine, -albumin, -apolipoprotein_a), 
        by = "eid"
      ) %>% 
      dplyr::inner_join(df_blood_counts_preprocessed, by = "eid") %>% 
      dplyr::inner_join(df_standard %>% dplyr::select(-hdl_cholesterol, -ldl_cholesterol, -total_cholesterol), by = "eid") 
  }
  
}



znormalize_data <- function(input_df,
                            znorm_data,
                            .keep = TRUE) {

  if (!("eid" %in% names(input_df))) {
    stop("Input data must have 'eid' column.")
  }

  df_znormed <- input_df %>%
    dplyr::select(eid, tidyselect::any_of(znorm_data$variable)) %>%
    tidyr::gather(key = "variable", value = "value", -eid) %>%
    dplyr::inner_join(znorm_data, by = "variable") %>%
    dplyr::mutate(value_znorm = (value - mean) / sd) %>%
    dplyr::select(eid, variable, value_znorm) %>%
    tidyr::spread(key = variable, value = value_znorm)


  if (.keep == TRUE) {
    df_znormed <- input_df %>%
      dplyr::select(-tidyselect::any_of(znorm_data$variable)) %>%
      dplyr::inner_join(df_znormed, by = "eid")
  }

  # Keep the original order of columns and rows
  df_out <- df_znormed %>%
    dplyr::select(tidyselect::any_of(names(input_df))) %>%
    dplyr::arrange(match(eid, input_df$eid))


  return(df_out)

}




