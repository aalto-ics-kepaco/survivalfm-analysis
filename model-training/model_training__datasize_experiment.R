on_cluster = T
setwd("~/survivalfm-analysis/model-training")

# Command line arguments --------------------------------------------------

if (on_cluster) {
  
  args <- commandArgs(TRUE)
  
  if (length(args) != 1) {
    stop(
      "Incorrect number of arguments",
      call. = FALSE
    )
  }
  array_index <- as.integer(args[1L])
  
} else {
  
  Sys.setenv('SLURM_CPUS_PER_TASK' = 5)
  
  array_index <- ..
  
}

library(tidyverse)
library(foreach)
library(survival)
library(glmnet)
library(survivalfm)
library(doParallel)
library(survMisc)

message("Libraries loaded. ")


# Setup -------------------------------------------------------------------


purrr::walk(list.files("utils", full.names = T), ~source(.x))


params <-
  list(
    seed = 12,
    survival_time = 10,
    nfolds_inner = 10,
    list_path = "model_list.tsv",
    result_lists = "result-lists-datasize-experiment",
    data_path = if (on_cluster) {...} else {...},  # Location of the data
    outputs_path = if (on_cluster) {...} else {...}, # Location to store results
    lambda_range = 10^seq(log10(1e-4), log10(1), length.out = 9) %>% round(5),
    k = 10
  )

dir.create(params$result_lists, showWarnings = F)
dir.create(params$outputs_path, showWarnings = F)

purrr::walk(
  .x = c(
    file.path(params$outputs_path, "cv-results"),
    file.path(params$outputs_path, "final-fits"),
    file.path(params$outputs_path, "predictions"),
    file.path(params$outputs_path, "cv-folds"),
    file.path(params$outputs_path, "shuffled-indices"),
    recursive = T
  ),
  .f = ~dir.create(.x)
)


# Set seed for reproducibility
set.seed(params$seed)

# Read model list --------------------------------------------------------

message(paste0("list path: ", params$list_path))

model_list <- readr::read_tsv(params$list_path)

model_list <- model_list %>%
  rowwise() %>%
  dplyr::mutate(
    identifier  = paste(
      c(stringr::str_replace_all(endpoint_name, " ", "_"),
        stringr::str_replace_all(input_name, " ", "_"),
        stringr::str_replace_all(training_method, " ", "_")
      ) %>% na.omit() %>% str_replace_all("\\+", "_"),
      collapse = "__")
  ) %>%
  ungroup()

stopifnot(nrow(model_list) == n_distinct(model_list$identifier))

model_row <- model_list %>% dplyr::filter(idx == !!array_index)


message(paste0("Training model ", array_index, "/", nrow(model_list), ": ", model_row$endpoint_name, " - ", model_row$input_name, "\nTraining method: ", model_row$training_method))



# Read data ---------------------------------------------------------------

df_input <- read_ukb_data(params, model_row) 

df_endpoint <- create_endpoints(model_list, params) %>% 
  dplyr::select(
    eid,
    incident_event = paste0("INCIDENT_", model_row$endpoint_var),
    event_agediff = paste0(model_row$endpoint_var, "_AGEDIFF")
  )

df_current <- df_input %>%
  dplyr::left_join(df_endpoint, by = "eid") %>%
  # Remove prevalent events
  dplyr::filter(event_agediff > 0) %>% 
  drop_na(incident_event, event_agediff) 


# Extract assessment center locations -------------------------------------

df_centers <- data.table::fread(paste0(params$data_path, "/preprocessed-data/baseline_characteristics.tsv")) %>% 
  dplyr::select(eid, uk_biobank_assessment_centre)

scotland_eids <- df_centers %>% 
  dplyr::filter(uk_biobank_assessment_centre %in% c("Glasgow", "Edinburgh")) %>% 
  dplyr::pull(eid)

# Scotland for testing
df_test <- df_current %>% filter(eid %in% scotland_eids) 
# England and Wales for training
df_train <- df_current %>% filter(!eid %in% scotland_eids)

stopifnot(length(intersect(df_test$eid, df_train$eid)) == 0)

rm(df_input)
rm(df_endpoint)
rm(df_centers)
rm(model_list)

gc()


# Numbers of individuals in training --------------------------------------

if (model_row$input_name == "Metabolomics biomarkers") {
  data_sizes <- c(10000, 25000, 50000, 100000, 150000, 200000)
} else if (model_row$input_name == "Standard risk factors") {
  data_sizes <- c(10000, 25000, 50000, 100000, 200000, 300000, 400000)
} else {
  data_sizes <- c(10000, 25000, 50000, 100000, 200000, 300000, 400000)
}


for (rep in seq(1:params$n_repeats)) {
  
  message(paste0("REPEAT: ", rep))
  
  for (n_idx in 1:length(data_sizes)) {
    
    n <- data_sizes[n_idx]
    
    message(paste0("N train: ", n))
    
    model_row_current <- model_row
    model_row_current$identifier <- paste0(model_row$identifier, "__n", n, "__rep", rep)
    
    if (file.exists(file.path(params$result_lists, paste0(model_row_current$identifier, ".csv")))) {
      message("Alrady exists, skipping...")
      next
    }
    
    # Subset the training data ------------------------------------------------
    
    shuffled_idx_path <- file.path(params$outputs_path, "shuffled-indices", paste0(model_row$input_name, "__", model_row$endpoint_name, "__rep", rep, ".rds"))
    
    if (file.exists(shuffled_idx_path)){
      shuffled_indices <- readRDS(shuffled_idx_path)
      print("Read shuffled indices from cache. ")
    } else {
      shuffled_indices <- sample(nrow(df_train))
      saveRDS(shuffled_indices, shuffled_idx_path)
    }
    
    df_train_sub <- df_train[shuffled_indices, ] %>% slice(1:n)
    
    # Create cross-validation folds -------------------------------------------

    cv_fold_path <- file.path(params$outputs_path, paste0("cv-folds/", model_row$endpoint_name, "__", model_row$input_name, "__n", n, "__rep", rep, ".rds")) %>%
      stringr::str_replace_all(" ", "_") %>%
      stringr::str_to_lower()

    if (file.exists(cv_fold_path)){
      folds <- readRDS(cv_fold_path)
      print("Read CV folds from cache. ")
    } else {
      folds <- caret::createFolds(df_train_sub$incident_event, k = params$nfolds_inner, list = TRUE, returnTrain = FALSE)
      saveRDS(folds, cv_fold_path)
      print("Created new CV folds. ")
    }

    # Cross-validation to optimize parameters ---------------------------------
    
    message("Cross-validation..")
    
    model_row_current$cv_checkpoints <- file.path(params$outputs_path, "cv-checkpoints", paste0(model_row_current$identifier))
    dir.create(model_row_current$cv_checkpoints, recursive = T)
    
    cv_result_path <- file.path(params$outputs_path, "cv-results", paste0(model_row_current$identifier, ".rds"))
    
    
    cv_results <- cross_validation(df_train = df_train_sub, folds = folds, params, model_row = model_row_current)
    saveRDS(cv_results, cv_result_path)
    

    
    # Select optimal regularization parameters --------------------------------
    
    all_params <- c("lambda1", "lambda2")
    missing_vars <- all_params[!all_params %in% names(cv_results)]
    
    for (var_name in missing_vars) {
      cv_results[[var_name]] <- NA
    }
    
    cv_results_summary <- cv_results %>%
      group_by(lambda1, lambda2) %>%
      summarize(mean_cindex = mean(cindex_val)) %>% 
      ungroup()
    
    best_params <- cv_results_summary[which.max(cv_results_summary$mean_cindex)[1],]
    
    message("Best parameters: ")
    message(paste0(names(best_params), ": ",  best_params, collapse = "\n"))

    model_row_current$best_lambda1 = best_params$lambda1
    model_row_current$best_lambda2 = best_params$lambda2


    # Preprocess data (scale) -------------------------------------------------
    
    is_binary <- function(x) {all(x %in% c(0, 1, NA))}
    binary_vars <- df_current %>% select(where(is_binary)) %>% names() 
    
    message("Z-normalizing data...")
    
    znorm_data <- df_train_sub %>%
      dplyr::select(-eid, -any_of(binary_vars), -incident_event, -event_agediff) %>%
      tidyr::gather(key = "variable", value = "value") %>%
      dplyr::group_by(.data$variable) %>%
      dplyr::summarise(
        mean = mean(.data$value, na.rm = T),
        sd = stats::sd(.data$value, na.rm = T),
        .groups = "keep"
      ) %>%
      dplyr::ungroup()
    
    df_train_preprocessed <- df_train_sub %>%
      znormalize_data(znorm_data =  znorm_data) 
    
    df_test_preprocessed <- df_test %>%
      znormalize_data(znorm_data =  znorm_data) 
    
    
    X_train <- df_train_preprocessed %>% tibble::column_to_rownames("eid") %>% dplyr::select(-incident_event, -event_agediff)
    status_train <- df_train_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(incident_event)
    times_train  <- df_train_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(event_agediff)

    # Impute only training data
    sink("NUL") # Sink messages from impute.knn
    X_train  <- impute::impute.knn(data = as.matrix(X_train), k = 10)$data %>% as.data.frame()
    sink()

    # Remove all samples with missing values from test
    df_test_preprocessed <- tidyr::drop_na(df_test_preprocessed)
    
    X_test <- df_test_preprocessed %>% tibble::column_to_rownames("eid") %>% dplyr::select(-incident_event, -event_agediff)
    status_test <- df_test_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(incident_event)
    times_test <- df_test_preprocessed  %>% tibble::column_to_rownames("eid") %>% pull(event_agediff)
    

    # Train final model -------------------------------------------------------

    model_row_current$final_fit_path <- file.path(params$outputs_path, "final-fits", paste0(model_row_current$identifier, ".rds"))
    
    message("\nTraining final model..")
    
    start_time <- Sys.time()
    
    final_fit <- fit_model(
      X = X_train,
      times = times_train,
      status = status_train,
      model_row =  model_row,
      model_params = list(
        lambda1 = best_params$lambda1,
        lambda2 = best_params$lambda2,
        k = params$k
      )
    )
    
    end_time <- Sys.time()
    training_time <- difftime(end_time, start_time, units = "mins")
    final_fit$training_time <- training_time
    saveRDS(final_fit, model_row_current$final_fit_path)
    
    
    # Evaluate performance ----------------------------------------------------
    
    message("Evaluating performance...")
    
    model_row_current$predictions_path <- file.path(params$outputs_path, "predictions", paste0(model_row_current$identifier, ".rds"))
    
    df_pred_test <- predict_by_model(X_test, model_row_current, final_fit) 
    
    baseline_survival <- survivalfm:::basesurv(
      time = times_train,
      event = status_train,
      lp = df_pred_train$lp,
      times.eval = params$survival_time
    )$base_surv

    df_pred_test <- df_pred_test %>%
      dplyr::mutate(
        abs_risk = (1 - baseline_survival^exp(.data$lp)),
        status = status_test,
        times = times_test
      )
    
    saveRDS(df_pred_test, model_row_current$predictions_path)
    
    cindex_test <- survival::concordance(survival::Surv(times, status) ~ lp, data = df_pred_test, reverse = T)$concordance
    cindex_test_se <- survival::concordance(survival::Surv(times, status) ~ lp, data = df_pred_test, reverse = T)$var %>% sqrt()
    
    fit_test <- survival::coxph(as.formula("survival::Surv(times, status) ~ lp"), data = df_pred_test)
    rsq_test <- survMisc::rsq(fit_test, sigD = NULL)$mev
    
    message(paste0("FINAL C-INDEX (test): ", cindex_test))
    message(paste0("FINAL R-sq (test): ", rsq_test))
    
    rm(X_test, status_test, times_test, df_pred_test)

    df_pred_train <- predict_by_model(X_train, model_row_current, final_fit)

    df_pred_train <- df_pred_train %>%
      dplyr::mutate(
        abs_risk = (1 - baseline_survival^exp(.data$lp))*100,
        status = status_train,
        times = times_train
      )

    cindex_train <- survival::concordance(survival::Surv(times, status) ~ lp, data = df_pred_train, reverse = T)$concordance
    cindex_train_se <- survival::concordance(survival::Surv(times, status) ~ lp, data = df_pred_train, reverse = T)$var %>% sqrt()

    fit_train <- survival::coxph(as.formula("survival::Surv(times, status) ~ lp"), data = df_pred_train)
    rsq_train <- survMisc::rsq(fit_train, sigD = NULL)$mev

    rm(X_train, times_train, df_pred_train)
    
    model_row_current <- model_row_current %>% 
      dplyr::mutate(
        rep = rep,
        n_train = n,
        n_events_train = sum(status_train),
        rsq_test = rsq_test,
        rsq_train = rsq_train,
        cindex_test = cindex_test,
        cindex_test_se = cindex_test_se,
        cindex_train = cindex_train,
        cindex_train_se = cindex_train_se,
        cindex_val = best_params$mean_cindex,
        final_model_training_time = as.numeric(final_fit$training_time)
      )
    
    readr::write_csv(model_row_current, file.path(params$result_lists, paste0(model_row_current$identifier, ".csv")))
    
  }
}
