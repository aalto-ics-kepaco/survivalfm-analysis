
add_interaction_terms <- function(X) {
  
  interaction_variables <- colnames(X)
  
  interaction_terms <- c()
  for (i in 1:(length(interaction_variables) - 1)) {
    for (j in (i + 1):length(interaction_variables)) {
      interaction_terms <- c(interaction_terms, paste0("`", interaction_variables[i], "`", "*", "`", interaction_variables[j], "`"))
    }
  }
  
  interaction_terms_formula <- paste(interaction_terms, collapse = " + ")
  
  formula <- as.formula(paste("~ . +", interaction_terms_formula))

  options(na.action = 'na.pass')
  
  X_out <-
    stats::model.matrix(formula, X, na.action = na.pass) %>% 
    as.data.frame() %>% 
    dplyr::select(-"(Intercept)")
  
  return(X_out)
}



add_qrisk_age_interaction_terms <- function(X) {
  
  # Add interaction terms
  interaction_variables <- colnames(X)[str_detect(colnames(X), "smoking|atrial|corti|migraine|kidney|lupus|hypertension|diab|body|family|systolic")]
  
  # Create the interaction terms with age
  interaction_terms_age_1 <- paste(paste0("`", interaction_variables, "`"), "age",  sep="*", collapse=" + ")
  interaction_terms_age_2 <- paste(paste0("`", interaction_variables, "`"), "age_1",  sep="*", collapse=" + ")
  formula <- as.formula(paste("~ . +", interaction_terms_age_1, "+", interaction_terms_age_2))
  
  # Add interaction terms in the data matrix
  options(na.action = 'na.pass')
  
  X_out <-
    stats::model.matrix(formula, X, na.action = na.pass) %>% 
    as.data.frame() %>% 
    dplyr::select(-"(Intercept)")
  
  return(X_out)
}


fit_model <- function(X, times, status, model_row, model_params) {

 if (model_row$training_method == "survivalFM") {

    print(paste0("Method: survivalFM, k = ", model_params$k))

    stopifnot(!is.null(model_params$lambda1))
    stopifnot(!is.null(model_params$lambda2))
    stopifnot(!is.null(model_params$k))
    
    
    fit <- survivalfm:::survivalfm(
      x = as.matrix(X),
      y = survival::Surv(times, status),
      rank = model_params$k,
      lambda1 = model_params$lambda1,
      lambda2 =   model_params$lambda2,
      maxiter = 2000
    )
    
 } else if (model_row$training_method == "Standard Cox regression") {
   
   print("Method: Cox regression")
   stopifnot(!is.null(model_params$lambda1))
   
   fit <-  glmnet::glmnet(
     x = as.matrix(X),
     y = survival::Surv(times, status),
     family = "cox",
     lambda = model_params$lambda1,
     alpha = 0, # L2 (ridge)
     standardize = FALSE # Already standardized
   )
   
 } else if (model_row$training_method == "Standard Cox regression with interaction terms") {
   
   print("Method: Cox regression with interaction terms")
   stopifnot(!is.null(model_params$lambda1))
   stopifnot(!is.null(model_params$lambda2))
   
   # Allow differential shrinkage for interaction effects
   penalty_factor <- rep(1, ncol(X))
   penalty_factor[str_detect(names(X), ":")] <- model_params$lambda2/model_params$lambda1
   
   fit <-  glmnet::glmnet(
     x = as.matrix(X),
     y = survival::Surv(times, status),
     family = "cox",
     lambda = model_params$lambda1,
     penalty.factor = penalty_factor,
     alpha = 0, # L2 (ridge)
     standardize = FALSE # Already standardized
   )
   
 } 
  return(fit)
}


predict_by_model <- function(X, model_row, fit) {
  
  if (model_row$training_method == "survivalFM") {
    df_pred <- tibble(
      eid = rownames(X),
      lp = survivalfm:::predict.survivalfm(fit, as.matrix(X), type = "link") %>% as.numeric()
    )
  } else if (str_detect(model_row$training_method, "Standard Cox regression")) {
    df_pred <- tibble(
      eid = rownames(X),
      lp =  stats::predict(fit, as.matrix(X), type = "link") %>% as.numeric()
    )
  } 
  
  return(df_pred)
}

