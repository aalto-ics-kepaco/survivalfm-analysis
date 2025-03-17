
cross_validation <- function(df_train, folds, params, model_row) {
  
  message("Model selection...")
  
  cv_results <- tibble()
  
  for (fold_idx in seq_along(folds)) {
    
    message(paste0("Fold: ", fold_idx))
    fold <- folds[[fold_idx]]
    
    df_train_cv <- df_train[-fold,]
    df_val_cv <- df_train[fold,]
    
    stopifnot(intersect(df_train_cv$eid, df_val_cv$eid) == 0)
    
  
    # Compute znorm data ------------------------------------------------------

    is_binary <- function(x) {all(x %in% c(0, 1, NA))}
    binary_vars <- df_train %>% select(where(is_binary)) %>% names() 
    
    message("Z-normalizing data...")
    
    # Scaling factors from training data
    znorm_data <- df_train_cv %>%
      dplyr::select(-eid, -any_of(binary_vars), -incident_event, -event_agediff) %>%
      tidyr::gather(key = "variable", value = "value") %>%
      dplyr::group_by(.data$variable) %>%
      dplyr::summarise(
        mean = mean(.data$value, na.rm = T),
        sd = stats::sd(.data$value, na.rm = T),
        .groups = "keep"
      ) %>%
      dplyr::ungroup()

    # Preprocess training data ------------------------------------------------

    # Z-normalize
    df_train_cv_preprocessed <- df_train_cv %>%
      znormalize_data(znorm_data =  znorm_data) 
    
    X_train_cv <- df_train_cv_preprocessed %>% tibble::column_to_rownames("eid") %>% dplyr::select(-incident_event, -event_agediff)
    status_train_cv <- df_train_cv_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(incident_event)
    times_train_cv <- df_train_cv_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(event_agediff)
    
    # Impute only training data
    sink("NUL") # Sink messages from impute.knn
    X_train_cv  <- impute::impute.knn(data = as.matrix(X_train_cv), k = 10)$data %>% as.data.frame()
    sink()

    # Preprocess validation data ----------------------------------------------

    # No imputation, remove all samples with missing values
    df_val_cv <- drop_na(df_val_cv)
    
    # Z-normalize
    df_val_cv_preprocessed <- df_val_cv %>%
      znormalize_data(znorm_data =  znorm_data) 
    
    X_val_cv <- df_val_cv_preprocessed %>% tibble::column_to_rownames("eid") %>% dplyr::select(-incident_event, -event_agediff)
    status_val_cv <- df_val_cv_preprocessed %>% tibble::column_to_rownames("eid") %>% pull(incident_event)
    times_val_cv <- df_val_cv_preprocessed  %>% tibble::column_to_rownames("eid") %>% pull(event_agediff)
    
    rm(df_train_cv, df_val_cv, df_train_cv_preprocessed, df_val_cv_preprocessed)
    
    stopifnot(length(intersect(rownames(X_val_cv), rownames(X_train_cv))) == 0)
    

    # Add interactions --------------------------------------------------------

    if (model_row$training_method == "Standard Cox regression with interaction terms") {
      # Add all possible interaction terms
      X_train_cv <- add_interaction_terms(X_train_cv)
      X_val_cv <- add_interaction_terms(X_val_cv)
    } 
    
    if (model_row$input_name == "QRISK3") {
      if (model_row$display_name == "Linear terms + QRISK3 interaction terms with age (standard Cox)") {
          # Add age interaction terms from QRISK3
          X_train_cv <- add_qrisk_age_interaction_terms(X_train_cv)
          X_val_cv <- add_qrisk_age_interaction_terms(X_val_cv)
        }
    }
    
    # Parameter grid ----------------------------------------------------------

    if (model_row$training_method == "survivalFM" | model_row$training_method == "Standard Cox regression with interaction terms") {
      param_combinations <- expand.grid(lambda1 = params$lambda_range, lambda2 = params$lambda_range) %>% 
        arrange(lambda1, lambda2)
    } else {
      param_combinations <- expand.grid(lambda1 = params$lambda_range) %>% 
        arrange(desc(lambda1))
    }
    
    
    # Register parallel backend -----------------------------------------------
    
    n_cores <- max(Sys.getenv('SLURM_CPUS_PER_TASK') %>% as.numeric(), 1, na.rm = T)
    message(paste0("Using ", n_cores, " cores"))
    
    cluster <- parallel::makeForkCluster(n_cores, type = "FORK")
    doParallel::registerDoParallel(cl = cluster)
    print(paste0("Number of workers available: ", foreach::getDoParWorkers()))
    
     
    results <- foreach(param_idx = 1:nrow(param_combinations), 
                       .combine = 'rbind', 
                       .inorder = FALSE,
                       .packages = c("survival", "survivalfm", "glmnet")
                       ) %dopar% {
      
      message(paste0("Processing parameter combination ", param_idx, "/", nrow(param_combinations)))
      param_comb <- param_combinations %>% dplyr::slice(!!param_idx)
      checkpoint_path <- paste0(model_row$cv_checkpoints, "/", paste0(names(param_comb), param_comb, collapse = "_"), "__fold", fold_idx, ".rds")
 
      if (file.exists(checkpoint_path)) {
        fit <- readRDS(checkpoint_path)
        message("Read cv checkpoint from files!")
      } else {
        start_time <- Sys.time()
        fit <- fit_model(
          X = X_train_cv,
          times = times_train_cv,
          status = status_train_cv,
          model_row = model_row,
          model_params = list(
            lambda1 = param_comb$lambda1,
            lambda2 = param_comb$lambda2,
            k = params$k
          )
        )
        end_time <- Sys.time()
        fit$model_training_time <- difftime(end_time, start_time, units = "mins")
        saveRDS(fit, checkpoint_path)
      }
      
      df_pred <- tibble(
        lp =  predict_by_model(X = X_val_cv, model_row = model_row, fit = fit)$lp,
        times = times_val_cv,
        status = status_val_cv
      )
      
      result <- param_comb %>%
        mutate(
          cindex_val = survival::concordance(survival::Surv(times, status) ~ lp, data = df_pred, reverse = T)$concordance,
          training_time = fit$model_training_time,
          checkpoint_path = checkpoint_path
        )
      
      rm(fit)
      gc()
      return(result)
    }
    
    parallel::stopCluster(cl = cluster)
    rm(list = "cluster")
    
    stopifnot(nrow(results) == nrow(param_combinations))
    
    cv_results <- cv_results %>% bind_rows(results %>% mutate(fold = fold_idx))
    
    rm(X_train_cv, X_val_cv)
    gc()
    
  }
  
  return(cv_results)
  
}

