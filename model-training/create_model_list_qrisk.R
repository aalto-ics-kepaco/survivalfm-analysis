library(tidyverse)

setwd("~/survivalfm-analysis/model-training")

df_selected_endpoints <-
  tibble::tribble(
    ~endpoint_var, ~endpoint_name, ~included_icd10_codes,
    "CVD", "Cardiovascular disease", "G45;I63;I64;I20;I21;I22;I23;I24;I25",
  )

# Create model slate ------------------------------------------------------

model_list <-
  rbind(
    df_selected_endpoints %>%
      mutate(display_name = "Linear terms (standard Cox)") %>% 
      mutate(training_method = "Standard Cox regression"),
    df_selected_endpoints %>%
      mutate(display_name = "Linear terms + QRISK3 interaction terms with age (standard Cox)") %>% 
      mutate(training_method = "Standard Cox regression"),
    df_selected_endpoints %>%
      mutate(display_name = "Linear terms + factorized interaction terms (survivalFM)") %>% 
      mutate(training_method = "survivalFM")
  ) %>% 
  dplyr::mutate(input_name = "QRISK3") %>% 
  dplyr::mutate(display_name = factor(display_name, levels = unique(.data$display_name))) %>%
  dplyr::mutate(training_method = factor(training_method, levels = unique(.data$training_method))) 


readr::write_tsv(model_list, "model_list_qrisk.tsv")
