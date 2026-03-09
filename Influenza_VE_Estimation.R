# ------------------------------------------------------------------------------
# Influenza VE Estimation: Bayesian hierarchical meta-analysis with JAGS
# ------------------------------------------------------------------------------

# Step 0: Set working directory to script location (RStudio or Rscript --file=)
tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  } else {
    script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(script_path) > 0) {
      path <- sub("^--file=", "", script_path)
      if (file.exists(path)) setwd(dirname(normalizePath(path)))
    }
  }
}, error = function(e) invisible(NULL))

# Optional: load showtext for fonts (Windows); skip if unavailable
suppressWarnings({
  if (requireNamespace("showtext", quietly = TRUE)) {
    tryCatch({
      library(showtext)
      font_add(family = "Times New Roman", regular = "C:/Windows/Fonts/times.ttf",
               bold = "C:/Windows/Fonts/timesbd.ttf", italic = "C:/Windows/Fonts/timesi.ttf",
               bolditalic = "C:/Windows/Fonts/timesbi.ttf")
      showtext_auto()
    }, error = function(e) invisible(NULL))
  }
})

# Step 1: Check and load required R packages
required_pkgs <- c("HDInterval", "ggplot2", "matrixStats", "dplyr", "magrittr", "readxl", "stringr", "tidyr", "metafor", "tidyverse", "corrplot")
missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing R packages: ", paste(missing, collapse = ", "), "\n  Install with: install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))")
}
library(HDInterval)
library(ggplot2)
library(matrixStats)
library(dplyr)
library(magrittr)
library(readxl)
library(stringr)
library(tidyr)
library(metafor)

# Lazy-load rjags/coda only when running JAGS (requires JAGS installed on system)
ensure_jags_loaded <- function() {
  if (!requireNamespace("rjags", quietly = TRUE)) stop("Install rjags (requires JAGS): install.packages(\"rjags\")")
  if (!requireNamespace("coda", quietly = TRUE)) stop("Install coda: install.packages(\"coda\")")
  suppressPackageStartupMessages({ library(coda); library(rjags) })
}

# Step 2: Data preparation
# Input file: model_input.xlsx (required sheets: country, ve2, global)
input_file <- "model_input.xlsx"
# 2a. Country-level covariates (sheet "country"): average 2022--2024, log-transform GDP/PopDens
df <- read_excel(input_file, sheet = "country")

mean_exclude_zero <- function(x) {
  if (all(x == 0, na.rm = TRUE)) {
    return(0)
  } else {
    x_nonzero <- x[x != 0]
    if (length(x_nonzero) == 0) return(0)
    return(mean(x_nonzero, na.rm = TRUE))
  }
}

mean_simple <- function(x) {
  mean(x, na.rm = TRUE)
}

df_new <- df %>%
  mutate(
    InfVax6m = apply(select(., InfVax6m_2022, InfVax6m_2023, InfVax6m_2024), 1, mean_exclude_zero),
    IVC65p   = apply(select(., IVC65p_2022, IVC65p_2023, IVC65p_2024), 1, mean_exclude_zero),
    PLGP     = apply(select(., PLGP_2022, PLGP_2023, PLGP_2024), 1, mean_exclude_zero),
    perGDP  = apply(select(., perGDP_2022, perGDP_2023, perGDP_2024), 1, mean_simple),
    p14p     = apply(select(., p14p_2022, p14p_2023, p14p_2024), 1, mean_simple),
    p65p     = apply(select(., p65p_2022, p65p_2023, p65p_2024), 1, mean_simple),
    PopDens  = apply(select(., PopDens_2022, PopDens_2023), 1, mean_simple),
    re       = apply(select(., re_2022, re_2023, re_2024), 1, mean_simple)
  )

data2 <- df_new %>%
  select(
    `Country Name`,
    `Country Code`,
    WB_region,
    BBS_pop_2022,
    AvgSAT,
    GDP,
    perGDP,
    InfVax6m,
    IVC65p,
    PLGP,
    p14p,
    p65p,
    PopDens,
    re,
    PM_2020,
    INT_ARR_2022
  )


data2 <- data2 %>%
  mutate(
    `Country Name` = str_trim(`Country Name`),
    `Country Code` = str_trim(`Country Code`)
  )

data2 <- data2 %>%
  mutate(INT_ARR_2022 = ifelse(is.na(INT_ARR_2022), 0, INT_ARR_2022)) %>%
  mutate(across(-c(`Country Name`, `Country Code`, WB_region), as.numeric))

data2 <- data2 %>%
  mutate(
    GDP = case_when(is.na(GDP) ~ NA_real_,GDP == 0 ~ 0,GDP > 0 ~ log(GDP),TRUE ~ log(1e-6)),
    perGDP = case_when(is.na(perGDP) ~ NA_real_,perGDP == 0 ~ 0,perGDP > 0 ~ log(perGDP),TRUE ~ log(1e-6)),
    PopDens = case_when(is.na(PopDens) ~ NA_real_,PopDens == 0 ~ 0,PopDens > 0 ~ log(PopDens),TRUE ~ log(1e-6))
  )

data2 <- data2 %>% rename(CountryCode = `Country Code`)

# 2b. Study-level VE (sheet "ve2"): convert to proportion and log(1-VE) for model
data1 <- read_excel(input_file, sheet = "ve2")

data1 <- data1 %>%
  mutate(
    ve_siv = ve_siv / 100,
    ve_siv_lower95ci = ve_siv_lower95ci / 100,
    ve_siv_upper95ci = ve_siv_upper95ci / 100,
    ve_siv_cidifference = ve_siv_upper95ci - ve_siv_lower95ci,
    logor_siv = log(1 - ve_siv)
  )

# 2c. Global country list and predictors (sheet "global"): all countries for prediction
df <- read_excel(input_file, sheet = "global")

df_new <- df %>%
  mutate(
    
    AvgSAT = apply(select(., AvgSAT_2022, AvgSAT_2023, AvgSAT_2024), 1, mean_simple),
    GDP = apply(select(., GDP_2022, GDP_2023, GDP_2024), 1, mean_exclude_zero),
    perGDP = apply(select(., perGDP_2022, perGDP_2023, perGDP_2024), 1, mean_exclude_zero),
    InfVax6m = apply(select(., InfVax6m_2022, InfVax6m_2023, InfVax6m_2024), 1, mean_exclude_zero),
    IVC65p   = apply(select(., IVC65p_2022, IVC65p_2023, IVC65p_2024), 1, mean_exclude_zero),
    PLGP     = apply(select(., PLGP_2022, PLGP_2023, PLGP_2024), 1, mean_exclude_zero),
    
    p14p     = apply(select(., p14p_2022, p14p_2023, p14p_2024), 1, mean_simple),
    p65p     = apply(select(., p65p_2022, p65p_2023, p65p_2024), 1, mean_simple),
    PopDens  = apply(select(., PopDens_2022, PopDens_2023), 1, mean_simple),
    re       = apply(select(., re_2022, re_2023, re_2024), 1, mean_simple)
  )

pre_data <- df_new %>%
  select(
    `Country Name`,
    `Country Code`,
    WB_region,
    BBS_pop_2022,
    AvgSAT,
    GDP, perGDP,
    InfVax6m,
    IVC65p,
    PLGP,
    p14p,
    p65p,
    PopDens,
    re,
    PM_2020,
    INT_ARR_2022,
    Income_region
  )

pre_data <- pre_data %>%
  mutate(
    `Country Name` = str_trim(`Country Name`),
    `Country Code` = str_trim(`Country Code`)
  )

pre_data <- pre_data %>%
  mutate(INT_ARR_2022 = ifelse(is.na(INT_ARR_2022), 0, INT_ARR_2022)) %>%
  mutate(across(-c(`Country Name`, `Country Code`, WB_region, Income_region), as.numeric))

pre_data <- pre_data %>%
  mutate(
    GDP = case_when(is.na(GDP) ~ NA_real_, GDP == 0 ~ 0, GDP > 0 ~ log(GDP), TRUE ~ log(1e-6)),
    perGDP = case_when(is.na(perGDP) ~ NA_real_, perGDP == 0 ~ 0, perGDP > 0 ~ log(perGDP), TRUE ~ log(1e-6)),
    PopDens = case_when(is.na(PopDens) ~ NA_real_, PopDens == 0 ~ 0, PopDens > 0 ~ log(PopDens), TRUE ~ log(1e-6))
  )
pre_data <- pre_data %>% rename(CountryCode = `Country Code`)
SD_pred <- apply(pre_data[, 4:15], 2, sd, na.rm = TRUE)

# Step 3: Build model inputs
# prep_pre_data: builds study outcomes (theta_hat, sigma2_hat_inv), study-by-country weight matrix (w),
#                and scaled predictors (x) + region IDs (super_region) for JAGS
prep_pre_data <- function(data1, pre_data, interactions = FALSE) {
  theta_hat <- data1$logor_siv
  se_ve <- (data1$ve_siv_upper95ci - data1$ve_siv_lower95ci) / 3.92
  denom <- (1.00 - data1$ve_siv)
  denom <- ifelse(denom <= 0 | !is.finite(denom), 1e-6, denom)
  var_log_rr <- (se_ve / denom)^2
  var_log_rr <- pmax(var_log_rr, 1e-10, na.rm = TRUE)
  sigma2_hat_inv <- 1.00 / var_log_rr
  sigma2_hat_inv <- pmin(sigma2_hat_inv, 1e6)
  sigma2_hat_inv[!is.finite(sigma2_hat_inv)] <- 1e4
  
  n <- length(theta_hat)
  w <- as.data.frame(data1[, 11:ncol(data1)])
  w[is.na(w)] <- 0
  
  all_countries <- pre_data$`CountryCode`
  
  w_full <- matrix(0, nrow = nrow(data1), ncol = length(all_countries))
  colnames(w_full) <- all_countries
  common <- intersect(colnames(w), all_countries)
  w_full[, match(common, all_countries)] <- as.matrix(w[, match(common, colnames(w))])

  w <- w_full
  n_c <- nrow(pre_data)
  
  pre_data <- pre_data %>%
    group_by(WB_region) %>%
    mutate(across(
      c(BBS_pop_2022, AvgSAT, perGDP, InfVax6m,
        IVC65p,
        PLGP,
        p14p, p65p, PopDens, re, PM_2020, INT_ARR_2022),
      ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)
    )) %>%
    ungroup()
  
  super_region <- rep(NA, times = n_c)
  unique_regions <- unique(pre_data$WB_region)
  for (j in seq_along(unique_regions)) {
    super_region[pre_data$WB_region == unique_regions[j]] <- j
  }
  
  n_s <- max(super_region, na.rm = TRUE)
  
  x1 <- pre_data$BBS_pop_2022
  x2 <- pre_data$AvgSAT
  x3 <- pre_data$perGDP
  x4 <- pre_data$InfVax6m
  x5 <- pre_data$IVC65p
  x6 <- pre_data$PLGP
  x7 <- pre_data$p14p
  x8 <- pre_data$p65p
  x9 <- pre_data$PopDens
  x10 <- pre_data$re
  x11 <- pre_data$PM_2020
  x12 <- pre_data$INT_ARR_2022
  
  preds <- cbind(
    scale(x1), scale(x2), scale(x3), scale(x4), 
    scale(x5), scale(x6),
    scale(x7), scale(x8), scale(x9), scale(x10),
    scale(x11), scale(x12)
  )
  result <- list(
    "theta_hat"      = theta_hat,
    "sigma2_hat_inv" = sigma2_hat_inv,
    "n"              = n,
    "n_s"            = n_s,
    "n_c"            = nrow(preds),
    "w"              = w,
    "super_region"   = super_region,
    "x"              = preds,
    "p_x"            = ncol(preds)
  )
  
  return(result)
}

# call_jags_full_model: runs the hierarchical model in JAGS (study theta ~ country eta1; country eta1 ~ covariates + region)
call_jags_full_model <- function(n, n_c, n_s, theta_hat, sigma2_hat_inv, x, p_x, w, super_region, n.chains, n.iter, thin) {
  ensure_jags_loaded()
  set.seed(8984)
  
  model_string <- "
  model {
    for (i in 1:n) {
      theta_hat[i] ~ dnorm(theta[i], sigma2_hat_inv[i])
      theta[i] ~ dnorm(mu_theta[i], sigma2_theta_inv)
      mu_theta[i] <- w[i,] %*% eta1
    }

    for (j in 1:n_c) {
      eta1[j] ~ dnorm(mu_eta1[j], sigma2_eta1_inv)
      mu_eta1[j] <- mu + x[j,] %*% beta + phi[super_region[j]]
    }

    for (k in 1:n_s) {
      phi[k] ~ dnorm(0.0, sigma2_phi_inv)
    }

    mu ~ dnorm(0.0, 0.0001)
    for (j in 1:p_x) {
      beta[j] ~ dnorm(0.0, sigma2_beta_inv)
    }
    sigma2_theta_inv ~ dgamma(0.01, 0.01)
    sigma2_beta_inv ~ dgamma(0.01, 0.01)
    sigma2_phi_inv ~ dgamma(0.01, 0.01)
    sigma2_eta1_inv ~ dgamma(0.01, 0.01)
  }
  "
  
  model_jags <- jags.model(
    textConnection(model_string),
    data = list(
      n = n,
      n_c = n_c,
      n_s = n_s,
      theta_hat = theta_hat,
      sigma2_hat_inv = sigma2_hat_inv,
      x = x,
      p_x = p_x,
      w = w,
      super_region = super_region
    ),
    n.chains = n.chains
  )
  
  update(model_jags, n.iter = n.iter)
  
  posterior_samples <- coda.samples(
    model_jags,
    variable.names = c(
      "theta", "eta1", "sigma2_theta_inv", "mu", 
      "beta", "sigma2_beta_inv", "phi", "sigma2_phi_inv", 
      "sigma2_eta1_inv"
    ),
    thin = thin,
    n.iter = n.iter
  )
  
  return(posterior_samples)
}

# eta_summary: from posterior eta1 (log scale), compute country-level VE median and 95% CrI
eta_summary <- function(data, posterior_samples){
  chains <- length(posterior_samples)
  final <- posterior_samples[[1]]
  for(j in 2:chains){
    final <- rbind(final, posterior_samples[[j]])
  }
  
  eta1 <- final[, (substring(colnames(final), 1, 4) == "eta1")]
  
  effectiveness_ve <- 1.00 - exp(eta1)
  
  output <- data.frame(
    CountryCode = data$CountryCode,
    VE_Effectiveness = matrixStats::colMedians(effectiveness_ve),
    LCI_VE_Effectiveness = matrixStats::colQuantiles(effectiveness_ve, probs = 0.025),
    UCI_VE_Effectiveness = matrixStats::colQuantiles(effectiveness_ve, probs = 0.975),
    Variance_Effectiveness = matrixStats::colVars(effectiveness_ve)
  )
  
  return(output)
}

# Step 4: Run main model (fit JAGS, then extract country-level VE)
model_input <- prep_pre_data(data1, pre_data)
posterior_samples <- call_jags_full_model(
  n = model_input$n,
  n_c = model_input$n_c,
  n_s = model_input$n_s,
  theta_hat = model_input$theta_hat,
  sigma2_hat_inv = model_input$sigma2_hat_inv,
  x = model_input$x,
  p_x = model_input$p_x,
  w = model_input$w,
  super_region = model_input$super_region,
  n.chains = 3,
  n.iter = 20000,
  thin = 10
)

VE_results_WB <- eta_summary(pre_data, posterior_samples)

# Step 5: Validation helpers (alternative data prep and JAGS model string; used by LOCO/K-fold/LOSO)
prep_data_CV <- function(data1, data2, interactions = FALSE) {
  theta_hat <- data1$logor_siv
  se_ve <- (data1$ve_siv_upper95ci - data1$ve_siv_lower95ci) / 3.92
  var_log_or <- (se_ve / (1.00 - data1$ve_siv))^2
  sigma2_hat_inv <- 1.00 / var_log_or
  n <- length(theta_hat)
  w <- as.matrix(data1[, 11:ncol(data1)])
  w[is.na(w)] <- 0
  n_c <- ncol(w)
  super_region <- as.integer(factor(data2$WB_region, levels = unique(data2$WB_region)))
  n_s <- length(unique(super_region))
  use_covariates <- TRUE  
  if (use_covariates) {
    predictors <- data2[, c("BBS_pop_2022","AvgSAT","perGDP","InfVax6m","IVC65p",
                            "PLGP","p14p","p65p","PopDens","re","PM_2020","INT_ARR_2022")]
    x <- scale(as.matrix(predictors))
    p_x <- ncol(x)
  } else {
    x <- NULL
    p_x <- NULL
  }
  
  list(
    theta_hat = theta_hat,
    sigma2_hat_inv = sigma2_hat_inv,
    n = n,
    n_c = n_c,
    n_s = n_s,
    w = w,
    super_region = super_region,
    x = x,
    p_x = p_x
  )
}

# Same JAGS model definition as in call_jags_full_model (used by LOCO/K-fold/LOSO)
build_model_string_full_model <- function() {
  return("
  model {
    for (i in 1:n) {
      theta_hat[i] ~ dnorm(theta[i], sigma2_hat_inv[i])
      theta[i] ~ dnorm(mu_theta[i], sigma2_theta_inv)
      mu_theta[i] <- w[i,] %*% eta1
    }

    for (j in 1:n_c) {
      eta1[j] ~ dnorm(mu_eta1[j], sigma2_eta1_inv)
      mu_eta1[j] <- mu + x[j,] %*% beta + phi[super_region[j]]
    }

    for (k in 1:n_s) {
      phi[k] ~ dnorm(0.0, sigma2_phi_inv)
    }

    mu ~ dnorm(0.0, 0.0001)
    for (j in 1:p_x) {
      beta[j] ~ dnorm(0.0, sigma2_beta_inv)
    }
    sigma2_theta_inv ~ dgamma(0.01, 0.01)
    sigma2_beta_inv ~ dgamma(0.01, 0.01)
    sigma2_phi_inv ~ dgamma(0.01, 0.01)
    sigma2_eta1_inv ~ dgamma(0.01, 0.01)
  }
  ")
}


run_leave_one_country_out_validation <- function(model_input, model_string, results_path) {
  library(foreach)
  library(doParallel)
  library(rjags)
  library(coda)
  
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  n_c <- ncol(model_input$w)
  cv_results <- foreach(cv = 1:n_c, .packages = c("rjags", "coda")) %dopar% {
    theta_hat_cv <- model_input$theta_hat
    
    studies_in_country <- which(model_input$w[, cv] > 0)
    
    if (length(studies_in_country) > 0) {
      theta_hat_cv[studies_in_country] <- NA
    }
    
    jags_data <- list(
      n = model_input$n,
      n_c = model_input$n_c,
      n_s = model_input$n_s,
      theta_hat = theta_hat_cv,
      sigma2_hat_inv = model_input$sigma2_hat_inv,
      x = model_input$x,
      p_x = model_input$p_x,
      w = model_input$w,
      super_region = model_input$super_region
    )
    
    model_jags <- jags.model(textConnection(model_string),
                             data = jags_data,
                             n.chains = 1)
    update(model_jags, n.iter = 5000)
    posterior_samples <- coda.samples(model_jags,
                                      variable.names = c("theta"),
                                      n.iter = 10000,
                                      thin = 10)
    return(posterior_samples)
  }
  
  stopCluster(cl)
  saveRDS(cv_results, results_path)
  return(cv_results)
}


run_kfold_validation <- function(model_input, model_string, results_path, k_folds = 5) {
  library(rjags)
  library(coda)
  
  set.seed(123)
  folds <- sample(rep(1:k_folds, length.out = model_input$n))
  cv_results <- list()
  
  for (fold in 1:k_folds) {
    test_indices <- which(folds == fold)
    theta_hat_train <- model_input$theta_hat
    theta_hat_train[test_indices] <- NA  
    
    jags_data <- list(
      n = model_input$n,
      n_c = model_input$n_c,
      n_s = model_input$n_s,
      theta_hat = theta_hat_train,
      sigma2_hat_inv = model_input$sigma2_hat_inv,
      w = model_input$w,
      super_region = model_input$super_region
    )
    
    if (!is.null(model_input$x) && model_input$p_x > 0) {
      jags_data$x <- model_input$x
      jags_data$p_x <- model_input$p_x
    } else {
      jags_data$x <- matrix(0, nrow = model_input$n_c, ncol = 0)
      jags_data$p_x <- 0
    }
    
    model_jags <- jags.model(textConnection(model_string),
                             data = jags_data,
                             n.chains = 1)
    update(model_jags, n.iter = 5000)
    posterior_samples <- coda.samples(model_jags,
                                      variable.names = c("theta"),
                                      n.iter = 10000,
                                      thin = 10)
    cv_results[[fold]] <- posterior_samples
  }
  
  saveRDS(cv_results, results_path)
  return(cv_results  )
}

# Step 6: Prepare for cross-validation (ensure w is matrix, get model string)
w_matrix <- as.matrix(model_input$w)
w_matrix[is.na(w_matrix)] <- 0
model_input$w <- w_matrix
model_string <- build_model_string_full_model()

# Step 7: Leave-One-Country-Out (LOCO) — for each country with data, set its studies to NA, refit, predict
run_loco_safe <- function(model_input,
                          model_string,
                          results_path,
                          parallel = FALSE,
                          n_cores = 2,
                          burnin = 5000,
                          n_iter = 10000,
                          thin = 10) {
  
  library(foreach)
  library(rjags)
  library(coda)
  
  n_c <- ncol(model_input$w)
  
  if (parallel) {
    library(doParallel)
    n_cores <- min(n_cores, parallel::detectCores() - 1)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    `%op%` <- `%dopar%`
    message("Running LOCO in PARALLEL mode with ", n_cores, " cores")
  } else {
    `%op%` <- `%do%`
    message("Running LOCO in SERIAL mode")
  }

  countries_with_data <- which(colSums(model_input$w > 0) > 0)
  
  cv_results <- foreach(
    cv = countries_with_data,
    .packages = c("rjags", "coda")
  ) %op% {
    
    cat("LOCO cv =", cv, "/", n_c, "\n")
    
    theta_hat_cv <- model_input$theta_hat
    studies_in_country <- which(model_input$w[, cv] > 0)
    
    if (length(studies_in_country) > 0) {
      theta_hat_cv[studies_in_country] <- NA
    }
    
    jags_data <- list(
      n = model_input$n,
      n_c = model_input$n_c,
      n_s = model_input$n_s,
      theta_hat = theta_hat_cv,
      sigma2_hat_inv = model_input$sigma2_hat_inv,
      x = model_input$x,
      p_x = model_input$p_x,
      w = model_input$w,
      super_region = model_input$super_region
    )
    
    t1 <- Sys.time()
    model_jags <- jags.model(
      textConnection(model_string),
      data = jags_data,
      n.chains = 1,
      quiet = TRUE
    )
    cat("  jags.model time:",
        round(difftime(Sys.time(), t1, units = "secs"), 1), "sec\n")
    
    t2 <- Sys.time()
    update(model_jags, n.iter = burnin)
    cat("  burn-in time:",
        round(difftime(Sys.time(), t2, units = "secs"), 1), "sec\n")
    
    t3 <- Sys.time()
    samples <- coda.samples(
      model_jags,
      variable.names = c("theta"),
      n.iter = n_iter,
      thin = thin
    )
    cat("  sampling time:",
        round(difftime(Sys.time(), t3, units = "secs"), 1), "sec\n")
    
    return(samples)
  }
  
  if (parallel) {
    stopCluster(cl)
  }
  
  saveRDS(list(cv_results = cv_results, countries_with_data = countries_with_data), results_path)
  message("LOCO finished, results saved to: ", results_path)
  
  return(list(cv_results = cv_results, countries_with_data = countries_with_data))
}

loco_results_obj <- run_loco_safe(
  model_input  = model_input,
  model_string = model_string,
  results_path = "results_loco_test.rds",
  parallel = FALSE  
)

loco_results <- loco_results_obj$cv_results
countries_with_data <- loco_results_obj$countries_with_data

# post_process_loco: from LOCO posterior theta, get study-level predicted VE and plot vs observed
post_process_loco <- function(cv_results, model_input, output_fig_path, plot_title, countries_with_data = NULL) {
  
  n <- model_input$n
  VE_post_median <- rep(NA_real_, n)
  
  if (is.null(countries_with_data))
    countries_with_data <- seq_along(cv_results)
  
  for (i in seq_along(cv_results)) {
    
    cv <- countries_with_data[i]
    posterior <- as.matrix(cv_results[[i]][[1]])
    
    studies <- which(model_input$w[, cv] > 0)
    
    cols <- paste0("theta[", studies, "]")
    cols_found <- intersect(cols, colnames(posterior))
    
    med <- apply(posterior[, cols_found, drop = FALSE], 2, median)
    
    idx <- as.integer(gsub("^theta\\[(\\d+)\\]$", "\\1", names(med)))
    
    VE_post_median[idx] <- 1 - exp(med)
  }
  
  VE_obs <- 1 - exp(model_input$theta_hat)
  
  corr <- cor(VE_obs, VE_post_median, use = "complete.obs")
  rmse <- sqrt(mean((VE_obs - VE_post_median)^2, na.rm = TRUE))
  
  df <- data.frame(VE_obs = VE_obs, VE_est = VE_post_median)
  
  library(ggplot2)
  
  p <- ggplot(df, aes(VE_obs, VE_est)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = plot_title,
      x = "Observed VE",
      y = "Predicted VE (Posterior Median)"
    ) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme_bw(base_size = 14)
  
  ggsave(output_fig_path, p, width = 6, height = 5)
  
  list(
    VE_post_median = VE_post_median,
    correlation = corr,
    rmse = rmse,
    plot = p
  )
}

loco_summary <- post_process_loco(
  loco_results_obj$cv_results,
  model_input,
  "VE_LOCO_WB.png",
  "LOCO CV Performance",
  countries_with_data = loco_results_obj$countries_with_data
)

print(loco_summary)

# Step 8: 5-fold cross-validation — split studies into 5 folds, leave one fold out, refit, predict test fold
library(rjags)
library(coda)

run_kfold_validation_with_indices <- function(model_input, model_string, results_path,
                                              k_folds = 5,
                                              n.chains = 1,
                                              n.iter_burnin = 1000,
                                              n.iter_sample = 2000,
                                              thin = 5,
                                              seed = 123) {
  set.seed(seed)
  n <- model_input$n
  folds <- sample(rep(1:k_folds, length.out = n))
  cv_results <- vector("list", k_folds)
  
  for (fold in 1:k_folds) {
    cat(sprintf("Running fold %d / %d\n", fold, k_folds))
    test_indices <- which(folds == fold)
    theta_hat_train <- model_input$theta_hat
    theta_hat_train[test_indices] <- NA
    
    jags_data <- list(
      n = model_input$n,
      n_c = model_input$n_c,
      n_s = model_input$n_s,
      theta_hat = theta_hat_train,
      sigma2_hat_inv = model_input$sigma2_hat_inv,
      w = model_input$w,
      super_region = model_input$super_region
    )
    if (!is.null(model_input$x) && model_input$p_x > 0) {
      jags_data$x <- model_input$x
      jags_data$p_x <- model_input$p_x
    } else {
      jags_data$x <- matrix(0, nrow = model_input$n_c, ncol = 0)
      jags_data$p_x <- 0
    }
    
    model_jags <- jags.model(textConnection(model_string), data = jags_data, n.chains = n.chains)
    update(model_jags, n.iter = n.iter_burnin)
    posterior_samples <- coda.samples(model_jags, variable.names = c("theta"), n.iter = n.iter_sample, thin = thin)
    
    posterior_mat <- as.matrix(posterior_samples)
    cv_results[[fold]] <- list(posterior = posterior_mat, test_indices = test_indices, fold = fold)
  }
  
  saveRDS(cv_results, results_path)
  return(cv_results)
}

kfold_results <- run_kfold_validation_with_indices(
  model_input = model_input,
  model_string = model_string,
  results_path = "results_5fold_with_indices.rds",
  k_folds = 5,
  n.chains = 1,
  n.iter_burnin = 1000,
  n.iter_sample = 2000,
  thin = 5
)

# post_process_kfold: from 5-fold posterior theta, get predicted VE per study and plot vs observed
post_process_kfold <- function(cv_results_or_path, theta_hat, output_fig_path, plot_title) {
  library(ggplot2)
  
  if (is.character(cv_results_or_path) && length(cv_results_or_path) == 1) {
    cv_results <- readRDS(cv_results_or_path)
  } else {
    cv_results <- cv_results_or_path
  }
  
  n <- length(theta_hat)
  VE_post_median <- rep(NA, n)
  
  for (fold_idx in seq_along(cv_results)) {
    fold_obj <- cv_results[[fold_idx]]
    posterior_mat <- as.matrix(fold_obj$posterior)
    test_idx <- fold_obj$test_indices
    
    param_names <- paste0("theta[", test_idx, "]")
    present <- intersect(param_names, colnames(posterior_mat))
    if (length(present) == 0) next
    
    medians <- apply(posterior_mat[, present, drop = FALSE], 2, median)
    
    for (nm in names(medians)) {
      idx <- as.integer(gsub("^theta\\[(\\d+)\\]$", "\\1", nm))
      VE_post_median[idx] <- medians[nm]
    }
  }
  
  VE_obs <- 1 - exp(theta_hat)
  VE_est <- 1 - exp(VE_post_median)
  
  valid_idx <- !is.na(VE_est) & !is.na(VE_obs)
  corr <- if (sum(valid_idx) >= 2) cor(VE_obs[valid_idx], VE_est[valid_idx], use = "complete.obs") else NA
  rmse <- sqrt(mean((VE_obs[valid_idx] - VE_est[valid_idx])^2, na.rm = TRUE))
  
  df <- data.frame(VE_obs = VE_obs[valid_idx], VE_est = VE_est[valid_idx])
  
  p <- ggplot(df, aes(x = VE_obs, y = VE_est)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    labs(title = plot_title,x = "Observed VE",y = "Predicted VE") +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme_bw(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", linewidth = 1),
      axis.line = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(output_fig_path, p, width = 6, height = 5)
  
  
  return(list(VE_post_median = VE_post_median, correlation = corr, rmse = rmse, plot = p))
}

kfold_summary <- post_process_kfold(
  kfold_results,
  model_input$theta_hat,
  "VE_5fold_WB.png",
  "5-Fold CV: Observed vs Predicted VE"
)
print(kfold_summary)

# Step 9: Leave-One-Study-Out (LOSO) — for each study, set it to NA, refit, predict that study
run_loso_validation <- function(model_input, model_string, results_path, n.chains=1, n.iter_burnin=2000, n.iter_sample=5000, thin=5) {
  library(rjags); library(coda)
  n <- model_input$n
  loso_results <- vector("list", n)
  for (i in 1:n) {
    cat("LOSO: leaving out study", i, "\n")
    theta_hat_cv <- model_input$theta_hat
    theta_hat_cv[i] <- NA
    jags_data <- list(
      n = model_input$n,
      n_c = model_input$n_c,
      n_s = model_input$n_s,
      theta_hat = theta_hat_cv,
      sigma2_hat_inv = model_input$sigma2_hat_inv,
      w = model_input$w,
      super_region = model_input$super_region
    )
    if (!is.null(model_input$x) && model_input$p_x > 0) {
      jags_data$x <- model_input$x; jags_data$p_x <- model_input$p_x
    } else { jags_data$x <- matrix(0, nrow=model_input$n_c, ncol=0); jags_data$p_x <- 0 }
    model_jags <- jags.model(textConnection(model_string), data=jags_data, n.chains=n.chains)
    update(model_jags, n.iter = n.iter_burnin)
    posterior_samples <- coda.samples(model_jags, variable.names = c("theta"), n.iter = n.iter_sample, thin = thin)
    loso_results[[i]] <- list(posterior = as.matrix(posterior_samples), left_out = i)
  }
  saveRDS(loso_results, results_path)
  return(loso_results)
}

loso_results <- run_loso_validation(
  model_input = model_input,
  model_string = model_string,
  results_path = "results_loso_newpriors.rds",
  n.chains = 1,
  n.iter_burnin = 1000,
  n.iter_sample = 2000,
  thin = 5
)

# Step 10: LOSO post-processing — extract posterior median log-OR for left-out study, convert to VE, compute correlation and RMSE
n <- length(loso_results)
log_or_post_median <- rep(NA, max(n, 1))

for (i in seq_along(loso_results)) {
  
  fold_obj <- loso_results[[i]]
  
  posterior_mat <- as.matrix(fold_obj$posterior)
  left_idx <- fold_obj$left_out

  col_idx <- grep(paste0("^theta\\[", left_idx, "\\]$"), colnames(posterior_mat))
  
  if (length(col_idx) == 1) {
    log_or_post_median[left_idx] <- median(posterior_mat[, col_idx], na.rm = TRUE)
  }
  
}

sum(is.na(log_or_post_median))

theta_hat <- model_input$theta_hat
VE_obs <- 1 - exp(theta_hat)
VE_est <- 1 - exp(log_or_post_median)

valid <- which(!is.na(VE_est) & !is.na(VE_obs))
pearson_r <- if(length(valid) >= 2) cor(VE_obs[valid], VE_est[valid], use="complete.obs") else NA
spearman_rho <- if(length(valid) >= 2) cor(VE_obs[valid], VE_est[valid], method="spearman", use="complete.obs") else NA
rmse <- sqrt(mean((VE_obs[valid] - VE_est[valid])^2, na.rm=TRUE))

cat("LOSO (corrected) Pearson r =", round(pearson_r, 4), "\n")
cat("LOSO (corrected) Spearman rho =", round(spearman_rho, 4), "\n")
cat("LOSO (corrected) RMSE =", round(rmse, 4), "\n")

# LOSO scatter plot: observed vs predicted VE with study labels for largest residuals
library(ggplot2)
idx_all <- seq_along(VE_obs)
df <- data.frame(
  VE_obs = VE_obs[valid],
  VE_est = VE_est[valid],
  idx    = idx_all[valid]
)
df$resid <- df$VE_obs - df$VE_est

p <- ggplot(df, aes(x = VE_obs, y = VE_est)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", color = "gray40"
  ) +
  geom_text(
    data = head(df[order(abs(df$resid), decreasing = TRUE), ], 6),
    aes(label = idx),
    vjust = -0.8,
    size = 3
  ) +
  labs(
    x = "Observed VE",
    y = "Predicted VE",
    title = "LOSO: Observed vs Predicted VE"
  ) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )


ggsave("LOSO_scatter_WB.png", p, width = 6, height = 5)


