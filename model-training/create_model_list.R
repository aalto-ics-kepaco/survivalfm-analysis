library(tidyverse)

setwd("~/survivalfm-analysis/model-training")

df_selected_endpoints <-
  tibble::tribble(
    ~endpoint_var, ~endpoint_name, ~included_icd10_codes,
    "ALZ", "Alzheimer's disease", "F00;G30",
    "VASC", "Vascular and other dementia", "F01;F02;F03",
    "LIVER", "Liver disease", "K70;K71;K72;K73;K74;K75;K76;K77",
    "LUNG_CANC", "Lung cancer", "C33;C34",
    "DIAB", "Type 2 diabetes", "E11;E14",
    "STR", "Stroke", "I60;I61;I63;I64",
    "MI", "Myocardial infarction", "I21;I22",
    "CKD", "Chronic kidney disease", "N18",
    "OST", "Osteoporosis", "M81;M80"
  )

# Create model list ------------------------------------------------------

model_list <-
  rbind(
    df_selected_endpoints %>% mutate(input_name = "Standard risk factors"),
    df_selected_endpoints %>% mutate(input_name = "Clinical biochemistry and blood counts"),
    df_selected_endpoints %>% mutate(input_name = "Metabolomics biomarkers"),
    df_selected_endpoints %>% mutate(input_name = "Polygenic risk scores")
  ) %>% 
  dplyr::mutate(input_name = factor(input_name, levels = unique(.data$input_name)))

model_list <-
  rbind(
    model_list %>% mutate(training_method = "Standard Cox regression"),
    model_list %>% mutate(training_method = "survivalFM")
  ) %>% 
  dplyr::bind_rows(
    model_list %>% filter(input_name == "Standard risk factors") %>% mutate(training_method = "Standard Cox regression with interaction terms")
  )

model_list <- model_list %>%
  dplyr::mutate(training_method = factor(training_method, levels = unique(.data$training_method))) %>%
  arrange(training_method, input_name, endpoint_var) %>% 
  dplyr::bind_rows(
    dplyr::bind_rows(
      df_selected_endpoints %>% mutate(input_name = "All") %>% dplyr::arrange(endpoint_var) %>% mutate(training_method = "Standard Cox regression"),
      df_selected_endpoints %>% mutate(input_name = "All") %>% dplyr::arrange(endpoint_var) %>% mutate(training_method = "survivalFM")
    ) 
  ) %>% 
  dplyr::mutate(idx = row_number()) 

readr::write_tsv(model_list, "model_list.tsv")

