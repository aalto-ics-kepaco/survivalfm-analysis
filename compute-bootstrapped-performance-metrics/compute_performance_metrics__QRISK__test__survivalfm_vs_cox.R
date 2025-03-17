library(tidyverse)
library(future)
library(furrr)

setwd("~/survivalfm-analysis/compute-bootstrapped-performance-metrics")

n_bootstrap <- 1000


# Read result files -------------------------------------------------------

result_files <- list.files("../model-training/result-lists-qrisk", full.names = T)

result_list  <- purrr::map_dfr(
  .x = result_files,
  .f =  ~readr::read_csv(.x, col_types = readr::cols())
)  %>% 
  dplyr::filter(display_name %in% c("Linear terms (standard Cox)" , "Linear terms + factorized interaction terms (survivalFM)"))


# Function to compute performance metrics ---------------------------------

metrics <- function(df) {
  
  # C-index
  cindex_std <- survival::concordance(survival::Surv(times, status) ~ lp_std, data = df, reverse = TRUE)$concordance
  cindex_new <- survival::concordance(survival::Surv(times, status) ~ lp_new, data = df, reverse = TRUE)$concordance
  delta_cindex <- cindex_new - cindex_std
  
  # NRI
  nri <- nricens::nricens(
    time = df$times,
    event = df$status,
    p.std = df$abs_risk_std,
    p.new = df$abs_risk_new,
    cut = 0.10,
    t0 = 10,
    updown = "category",
    niter = 1
  )$nri
  
  # Royston's R-squared
  fit_std <- survival::coxph(
    as.formula("survival::Surv(time = times, event = status) ~ lp_std"),
    data = df
  )
  
  fit_new <- survival::coxph(
    as.formula("survival::Surv(time = times, event = status) ~ lp_new"),
    data = df
  )
  rsq_std <- survMisc::rsq(fit_std, sigD = NULL)
  rsq_new <- survMisc::rsq(fit_new, sigD = NULL)
  
  delta_rsq <- rsq_new$mev - rsq_std$mev
  
  return(
    c(
      cindex_std = cindex_std, 
      cindex_new = cindex_new, 
      rsq_std = rsq_std$mev,
      rsq_new = rsq_new$mev,
      delta_cindex = delta_cindex, 
      delta_rsq = delta_rsq, 
      nri_overall = nri["NRI", "Estimate"], 
      nri_events = nri["NRI+", "Estimate"], 
      nri_non_events = nri["NRI-", "Estimate"]
    )
  )
}


# Read data ---------------------------------------------------------------

# Standard Cox regression
cox_path <- result_list %>% 
  filter(display_name == "Linear terms (standard Cox)") %>% 
  pull(predictions_path) %>% 
  str_remove("/m/cs/scratch/ukbiobank-kepaco/tmp-heli/")

df_std <- readRDS(paste0("../model-training/", cox_path))

# survivalFM
survivalfm_path <- result_list %>% 
  filter(display_name == "Linear terms + factorized interaction terms (survivalFM)") %>% 
  pull(predictions_path) %>% 
  str_remove("/m/cs/scratch/ukbiobank-kepaco/tmp-heli/")

df_new <- readRDS(paste0("../model-training/", survivalfm_path))

df <-
  dplyr::full_join(
    df_new %>% select(eid, status, times, lp_new = lp, abs_risk_new = abs_risk),
    df_std %>% select(eid, status, times, lp_std = lp, abs_risk_std = abs_risk),
    by = c("eid", "status", "times")
  )

stopifnot(nrow(df) == nrow(df_new))
stopifnot(nrow(df) == nrow(df_std))

rm(df_new, df_std)


# Bootsrapping ------------------------------------------------------------


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

# Confidence intervals 
ci_cindex_std <- boot::boot.ci(boot_res, type = "perc", index = 1) # C-index
ci_cindex_new <- boot::boot.ci(boot_res, type = "perc", index = 2) # C-index
ci_rsq_std <- boot::boot.ci(boot_res, type = "perc", index = 3) # R-sq
ci_rsq_new <- boot::boot.ci(boot_res, type = "perc", index = 4) # R-sq
ci_delta_cindex <- boot::boot.ci(boot_res, type = "perc", index = 5) # Delta cindex
ci_delta_rsq <- boot::boot.ci(boot_res, type = "perc", index = 6) # Delta R-sq
ci_nri_overall <- boot::boot.ci(boot_res, type = "perc", index = 7) # NRI overall
ci_nri_events <- boot::boot.ci(boot_res, type = "perc", index = 8) # NRI events
ci_nri_nonevents <- boot::boot.ci(boot_res, type = "perc", index = 9) # NRI non-events


df_combined <- tibble(
  std_model = "Cox regression",
  new_model = "survivalFM",
  cindex_std = boot_res$t0[1],
  cindex_std_ci_lower = ci_cindex_std$perc[4],
  cindex_std_ci_upper = ci_cindex_std$perc[5],
  cindex_new = boot_res$t0[2],
  cindex_new_ci_lower= ci_cindex_new$perc[4],
  cindex_new_ci_upper = ci_cindex_new$perc[5],
  rsq_std = boot_res$t0[3],
  rsq_std_ci_lower = ci_rsq_std$perc[4],
  rsq_std_ci_upper = ci_rsq_std$perc[5],
  rsq_new = boot_res$t0[4],
  rsq_new_ci_lower = ci_rsq_new$perc[4],
  rsq_new_ci_upper = ci_rsq_new$perc[5],
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



# Save results to CSV
df_combined %>%
  readr::write_csv("../performance-metrics/performance_metrics__QRISK__test__survivalfm_vs_cox.csv")

