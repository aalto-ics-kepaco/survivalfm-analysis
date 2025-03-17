library(tidyverse)
library(future)
library(furrr)

setwd("~/survivalfm-analysis/compute-bootstrapped-performance-metrics")

n_bootstrap <- 1000


# Read result files -------------------------------------------------------

result_files <- list.files("../model-training/result-slates", full.names = T)

result_slate  <- purrr::map_dfr(
  .x = result_files,
  .f =  ~readr::read_csv(.x, col_types = readr::cols())
)  %>% 
  dplyr::filter(input_name != "All") %>% 
  dplyr::filter(training_method %in% c("Standard Cox regression", "survivalFM"))


# Function to compute performance metrics ---------------------------------

metrics <- function(df) {
  
  # C-index
  cindex_cox <- survival::concordance(survival::Surv(times, status) ~ lp_cox, data = df, reverse = TRUE)$concordance
  cindex_survivalfm <- survival::concordance(survival::Surv(times, status) ~ lp_survivalfm, data = df, reverse = TRUE)$concordance
  delta_cindex <- cindex_survivalfm - cindex_cox
  
  # NRI
  nri <- nricens::nricens(
    time = df$times,
    event = df$status,
    p.std = df$abs_risk_cox,
    p.new = df$abs_risk_survivalfm,
    cut = 0,
    t0 = 10,
    updown = "diff",
    niter = 1
  )$nri
  
  # Royston's R-squared
  fit_cox <- survival::coxph(
    as.formula("survival::Surv(time = times, event = status) ~ lp_cox"),
    data = df
  )
  
  fit_survivalfm <- survival::coxph(
    as.formula("survival::Surv(time = times, event = status) ~ lp_survivalfm"),
    data = df
  )
  rsq_cox <- survMisc::rsq(fit_cox, sigD = NULL)
  rsq_survivalfm <- survMisc::rsq(fit_survivalfm, sigD = NULL)
  
  delta_rsq <- rsq_survivalfm$mev - rsq_cox$mev
  
  return(
    c(
      cindex_cox = cindex_cox, 
      cindex_survivalfm = cindex_survivalfm, 
      rsq_cox = rsq_cox$mev,
      rsq_survivalfm = rsq_survivalfm$mev,
      delta_cindex = delta_cindex, 
      delta_rsq = delta_rsq, 
      nri_overall = nri["NRI", "Estimate"], 
      nri_events = nri["NRI+", "Estimate"], 
      nri_non_events = nri["NRI-", "Estimate"]
    )
  )
}


# Run bootstrap -----------------------------------------------------------

# Set up caching directory
cache_path <- "cache/performance_metrics__test__survivalfm_vs_cox"
dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)

start_time <- Sys.time()



# Run the combined function with bootstrapping
df_combined <- 
  purrr::map2_dfr(
    .x = result_list %>% select(input_name, endpoint_name) %>% distinct() %>% pull(input_name),
    .y = result_slate %>% select(input_name, endpoint_name) %>% distinct() %>% pull(endpoint_name),
    .f = function(input_name, endpoint_name) {
      
      message(paste0(input_name, "-", endpoint_name))
      
      # Read data ---------------------------------------------------------------
      
      # Standard Cox regression
      cox_path <- result_slate %>% 
        filter(input_name == !!input_name, endpoint_name == !!endpoint_name, training_method == "Standard Cox regression") %>% 
        pull(predictions_path) %>% 
        str_remove("/m/cs/scratch/ukbiobank-kepaco/tmp-heli/")
      
      df_cox <- readRDS(paste0("../model-training/", cox_path))
      
      # survivalFM
      survivalfm_path <- result_slate %>% 
        filter(input_name == !!input_name, endpoint_name == !!endpoint_name, training_method == "survivalFM") %>% 
        pull(predictions_path) %>% 
        str_remove("/m/cs/scratch/ukbiobank-kepaco/tmp-heli/")
      
      df_survivalfm <- readRDS(paste0("../model-training/", survivalfm_path))
      
      df <-
        dplyr::full_join(
          df_survivalfm %>% select(eid, status, times, lp_survivalfm = lp, abs_risk_survivalfm = abs_risk),
          df_cox %>% select(eid, status, times, lp_cox = lp, abs_risk_cox = abs_risk),
          by = c("eid", "status", "times")
        )
      
      stopifnot(nrow(df) == nrow(df_survivalfm))
      stopifnot(nrow(df) == nrow(df_cox))
      
      rm(df_survivalfm, df_cox)
      
      
      # Bootsrapping ------------------------------------------------------------
      
      filename <- paste0(cache_path, "/bootstrap__", str_remove_all(input_name, "\\n"), "__", endpoint_name, "__", n_bootstrap, ".rds")
      
      if (file.exists(filename)) {
        # Read from cache if already computed
        boot_res <- readRDS(filename)
        message("Already exists, skipping...")
      } else {
        set.seed(123)
        # Run the bootstrap 
        boot_res <-  boot::censboot(
          data = df, 
          statistic = metrics, 
          R = n_bootstrap, 
          index = match(c("times", "status"), names(df)),
          parallel = "multicore",
          ncpus = 3
        )
        saveRDS(boot_res, filename)
      }
      
      
      # Confidence intervals 
      ci_cindex_cox <- boot::boot.ci(boot_res, type = "perc", index = 1) # C-index
      ci_cindex_survivalfm <- boot::boot.ci(boot_res, type = "perc", index = 2) # C-index
      ci_rsq_cox <- boot::boot.ci(boot_res, type = "perc", index = 3) # R-sq
      ci_rsq_survivalfm <- boot::boot.ci(boot_res, type = "perc", index = 4) # R-sq
      ci_delta_cindex <- boot::boot.ci(boot_res, type = "perc", index = 5) # Delta cindex
      ci_delta_rsq <- boot::boot.ci(boot_res, type = "perc", index = 6) # Delta R-sq
      ci_nri_overall <- boot::boot.ci(boot_res, type = "perc", index = 7) # NRI overall
      ci_nri_events <- boot::boot.ci(boot_res, type = "perc", index = 8) # NRI events
      ci_nri_nonevents <- boot::boot.ci(boot_res, type = "perc", index = 9) # NRI non-events
      
      
      tibble(
        input = input_name,
        endpoint = endpoint_name,
        cindex_cox = boot_res$t0[1],
        cindex_cox_ci_lower = ci_cindex_cox$perc[4],
        cindex_cox_ci_upper = ci_cindex_cox$perc[5],
        cindex_survivalfm = boot_res$t0[2],
        cindex_survivalfm_ci_lower= ci_cindex_survivalfm$perc[4],
        cindex_survivalfm_ci_upper = ci_cindex_survivalfm$perc[5],
        rsq_cox = boot_res$t0[3],
        rsq_cox_ci_lower = ci_rsq_cox$perc[4],
        rsq_cox_ci_upper = ci_rsq_cox$perc[5],
        rsq_survivalfm = boot_res$t0[4],
        rsq_survivalfm_ci_lower = ci_rsq_survivalfm$perc[4],
        rsq_survivalfm_ci_upper = ci_rsq_survivalfm$perc[5],
        delta_cindex = boot_res$t0[5],
        delta_cindex_ci_lower = ci_delta_cindex$perc[4],
        delta_cindex_ci_upper = ci_delta_cindex$perc[5],
        delta_rsq = boot_res$t0[6],
        delta_rsq_ci_lower = ci_delta_rsq$perc[4],
        delta_rsq_ci_upper = ci_delta_rsq$perc[5],
        nri_overall = boot_res$t0[7],
        nri_overall_ci_lower = ci_nri_overall$perc[4],
        nri_overall_ci_upper = ci_nri_overall$perc[5],
        nri_events = boot_res$t0[8],
        nri_events_ci_lower = ci_nri_events$perc[4],
        nri_events_ci_upper = ci_nri_events$perc[5],
        nri_nonevents = boot_res$t0[9],
        nri_nonevents_ci_lower = ci_nri_nonevents$perc[4],
        nri_nonevents_ci_upper = ci_nri_nonevents$perc[5],
        n = nrow(df),
        n_events = sum(df$status)
      )
    }
  )



# Save results to CSV
df_combined %>%
  readr::write_csv("../performance-metrics/performance_metrics__test__survivalfm_vs_cox.csv")

