# variables pre-setting
{
  cenPropAll <- c() # censoring proportion of the whole source RCT sample
  cenPropTrt <- c() # censoring proportion of the treatment arm in source RCT sample
  cenPropCtrl <- c() # censoring proportion of the control arm in source RCT sample
  
  # HR estimates of each replicate, when applying TADA with censoring adjustment
  HR_cen <- numeric()
  # sd(censoringHR_Boot)
  HR_cen_SE <- numeric()
  
  # HR estimates when applying TADA without censoring adjustment (only MoM adjustment ) 
  HR_noncen <- numeric()
  # sd(noncensoringHR_Boot) 
  HR_noncen_SE <- numeric()
  
  trueHR <- numeric() # simulated true HR generated with target IPD
  trueHR_SE <- numeric()
  
  # mean(censoringHR)
  HR_cen_average <- c()
  # mean
  HR_cen_SE_average <- c()
  
  # mean
  HR_noncen_average <- c()
  # mean
  HR_noncen_SE_average <- c()
  
  # mean
  trueHR_average <- c()
  # mean
  trueHR_SE_average <- c()
  
  # Bootstrapping
  CI_lower_cen <- c()
  # Bootstrapping
  CI_upper_cen <- c()
  
  # Bootstrapping 
  CI_lower_noncen <- c()
  # Bootstrapping 
  CI_upper_noncen <- c()
  
  CI_lower_cen_average <- c()
  CI_upper_cen_average <- c()
  
  CI_lower_noncen_average <- c()
  CI_upper_noncen_average <- c()
  
  CI_lower_true <- c()
  CI_upper_true <- c()
  
  CI_lower_true_average <- c()
  CI_upper_true_average <- c()
  
  bias_HR_cen <- c()
  coverage_HR_cen <- c()
  
  bias_HR_noncen <- c()
  coverage_HR_noncen <- c()
  
  Percentile_Lower_Bias_HR_censor <- c()
  Percentile_Upper_Bias_HR_censor <- c()
  
  Percentile_Lower_Bias_HR_noncensor <- c()
  Percentile_Upper_Bias_HR_noncensor <- c()
  
  
  TADA_Results <- data.frame()
}

# ----- Simulated data generation ----- #
generateTestDataTADA_TTE <- function(nStudyPopulation = 50000, 
                                     nStudy = 200,
                                     
                                     nTargetPopulation = 200000,  
                                     nTarget = 1000,
                                     
                                     alpha_event = 1.5,
                                     lambda_event = 0.1,
                                     
                                     beta_treatment = -0.4,
                                     beta_X1_trt = 0.5,
                                     beta_X2_trt = -0.3,
                                     beta_X3_trt = 0.2,
                                     beta_X1trt_all = c(1.2, 0.58, -0.9, -0.6), # 10 scenarios setting
                                     beta_X2trt_all = c(0.35, 0.2, -0.7, -0.1),
                                     beta_X3trt_all = c(1.3, 0.45, -0.3, -0.85),
                                     
                                     alpha_censor = 1.5,
                                     lambda_censor = 0.001,
                                     
                                     beta_trt_cen = 0.25,
                                     beta_X1_cen = -0.2,
                                     beta_X2_cen = 0.4,
                                     beta_X3_cen = -0.1,
                                     
                                     base_risk_cen = base_risk_cen) {
  
  studyData_All <- list()
  eventTime_All <- list()
  targetPopulationData_All <- list()
  aggregateTargetData_All <- list()
  
  # Step 0: Generate a common study population pool (X1, X2, X3)
  studyPopulationData <- data.frame(
    ID = 1:nStudyPopulation,
    X1 = rbinom(nStudyPopulation, 1, 0.45),
    X2 = rbinom(nStudyPopulation, 1, 0.65),
    X3 = rnorm(nStudyPopulation, 0, 1)
  )
  
  # --- Generate 10 fully independent studyData and targetData scenarios --- #
  for (i in 1:10) {
    # Determine which beta_trt interaction changes in this scenario
    if (i <= 4) {
      beta_X1trt <- beta_X1trt_all[i]
      beta_X2trt <- beta_X2trt_all[1]
      beta_X3trt <- beta_X3trt_all[1]
    } else if (i <= 7) {
      beta_X1trt <- beta_X1trt_all[1]
      beta_X2trt <- beta_X2trt_all[i - 3]
      beta_X3trt <- beta_X3trt_all[1]
    } else {
      beta_X1trt <- beta_X1trt_all[1]
      beta_X2trt <- beta_X2trt_all[1]
      beta_X3trt <- beta_X3trt_all[i - 6]
    }
    
    # ---------------- Study Data ------------------
    sampled_idx <- sample(1:nStudyPopulation, nStudy) # sample from population: index
    studyData <- studyPopulationData[sampled_idx, ]   # sample following index
    studyData$treatment <- rbinom(nStudy, 1, 0.5)     # treatment indicator: RCT setting
    
    # linear predictor 
    linpred_event <- with(studyData,
                          
                          beta_X1_trt * X1 + beta_X2_trt * X2 + beta_X3_trt * X3 +
                            beta_treatment * treatment +
                            beta_X1trt * (treatment * X1) +
                            beta_X2trt * (treatment * X2) +
                            beta_X3trt * (treatment * X3)
    )
    
    linpred_censor <- with(studyData,
                           
                           base_risk_cen + beta_trt_cen * treatment +
                             beta_X1_cen * X1 + beta_X2_cen * X2 + beta_X3_cen * X3
    )
    
    # event/censor time generation
    eventTime <- (-log(runif(nStudy)) / (lambda_event * exp(linpred_event)))^(1 / alpha_event) # T
    censoringTime <- (-log(runif(nStudy)) / (lambda_censor * exp(linpred_censor)))^(1 / alpha_censor) # C
    obsTime <- pmin(eventTime, censoringTime) # U = min(T, C)
    
    eventObserved <- as.integer(eventTime <= censoringTime) # event indicator
    censor <- 1 - eventObserved # censor indicator
    
    studyData_All[[i]] <- data.frame(
      studyData,
      eventTime = eventTime,
      censoringTime = censoringTime,
      obsTime = obsTime,
      eventObserved = eventObserved,
      censor = censor
    )
    
    # ---------------- Target Data ------------------
    targetPopulationData <- data.frame(
      ID = 1:nTargetPopulation,
      X1 = rbinom(nTargetPopulation, 1, 0.35),
      X2 = rbinom(nTargetPopulation, 1, 0.55),
      X3 = rnorm(nTargetPopulation, 0.21, 1.5)
    )
    targetPopulationData$treatment <- rbinom(nTargetPopulation, 1, 0.5)
    
    linpred_target <- with(targetPopulationData,
                           beta_X1_trt * X1 + beta_X2_trt * X2 + beta_X3_trt * X3 +
                             beta_treatment * treatment +
                             beta_X1trt * (treatment * X1) +
                             beta_X2trt * (treatment * X2) +
                             beta_X3trt * (treatment * X3)
    )
    
    eventTime_target <- (-log(runif(nTargetPopulation)) / (lambda_event * exp(linpred_target)))^(1 / alpha_event)
    
    targetPopulationData_All[[i]] <- data.frame(targetPopulationData, eventTime = eventTime_target)
    
    # Target AgD Summary (aggregated)
    targetSample_idx <- sample(1:nTargetPopulation, nTarget)
    targetSample <- targetPopulationData[targetSample_idx, ]
    aggregateTargetData_All[[i]] <- data.frame(
      N = nTarget,
      X1_prop = mean(targetSample$X1),
      X2_prop = mean(targetSample$X2),
      X3_mean = mean(targetSample$X3),
      X3_sd = sd(targetSample$X3)
    )
  }
  
  return(list(
    studyData = studyData_All,
    aggregateTargetData = aggregateTargetData_All,
    IPDTargetPopulationData = targetPopulationData_All
  ))
}  # end of function




# -------- function to calculate uncensored probability
calculateUncensoredProb <- function(time, linear_predictor, baseCumHaz) {
  idx <- which.min(abs(baseCumHaz$time - time))
  surv_prob <- exp(- baseCumHaz$hazard[idx] * exp(linear_predictor))
  return(surv_prob)
}


# Sensitivity Analysis: DGM for abnormal values in target population only
generateTestDataTADA_TTE_ABV <- function(
    nStudyPopulation   = 50000,
    nStudy             = 200,
    
    nTargetPopulation  = 200000,
    nTarget            = 1000,
    
    alpha_event        = 1.5,
    lambda_event       = 0.1,
    
    beta_treatment     = -0.4,
    beta_X1_trt        = 0.5,
    beta_X2_trt        = -0.3,
    beta_X3_trt        = 0.2,
    beta_X1trt_all     = c(1.2, 0.58, -0.9, -0.6),
    beta_X2trt_all     = c(0.35, 0.2, -0.7, -0.1),
    beta_X3trt_all     = c(1.3, 0.45, -0.3, -0.85),
    
    alpha_censor       = 1.5,
    lambda_censor      = 0.001,
    
    beta_trt_cen       = 0.25,
    beta_X1_cen        = -0.2,
    beta_X2_cen        = 0.4,
    beta_X3_cen        = -0.1,
    
    base_risk_cen,
    # abnormal value DGM
    ABV_Type               = "None",   # "None","X2_flip","X3_extreme","both"
    X2_Flips_Prop          = NA,       # e.g. 0.01 or 0.1
    X3_Extremes_Prop       = NA,       # e.g. 0.01 or 0.1
    X3_Extremes_Multiplier = NA        # e.g. 2 or 10
) {
  
  frac_X2 <- if (ABV_Type %in% c("X2_flip","both")) X2_Flips_Prop else 0
  frac_X3 <- if (ABV_Type %in% c("X3_extreme","both")) X3_Extremes_Prop else 0
  mult_X3 <- if (ABV_Type %in% c("X3_extreme","both")) X3_Extremes_Multiplier else 0
  
  studyData_All                <- list()
  aggregateTargetData_obs      <- list()  # abv AgD
  aggregateTargetData_true     <- list()  # clean AgD
  IPDTargetPopulationData_obs  <- list()  # abv IPD
  IPDTargetPopulationData_true <- list()  # clean IPD
  
  studyPopulationData <- data.frame(
    ID = 1:nStudyPopulation,
    X1 = rbinom(nStudyPopulation, 1, 0.45),
    X2 = rbinom(nStudyPopulation, 1, 0.65),
    X3 = rnorm(nStudyPopulation, 0, 1)
  )
  
  for (i in 1:10) {
    
    if (i <= 4) {
      beta_X1trt <- beta_X1trt_all[i]
      beta_X2trt <- beta_X2trt_all[1]
      beta_X3trt <- beta_X3trt_all[1]
    } else if (i <= 7) {
      beta_X1trt <- beta_X1trt_all[1]
      beta_X2trt <- beta_X2trt_all[i - 3]
      beta_X3trt <- beta_X3trt_all[1]
    } else {
      beta_X1trt <- beta_X1trt_all[1]
      beta_X2trt <- beta_X2trt_all[1]
      beta_X3trt <- beta_X3trt_all[i - 6]
    }
    
    sampled_idx <- sample(1:nStudyPopulation, nStudy)
    studyData   <- studyPopulationData[sampled_idx, ]
    studyData$treatment <- rbinom(nStudy, 1, 0.5)
    
    linpred_event <- with(studyData,
                          beta_X1_trt * X1 +
                            beta_X2_trt * X2 +
                            beta_X3_trt * X3 +
                            beta_treatment * treatment +
                            beta_X1trt * (treatment * X1) +
                            beta_X2trt * (treatment * X2) +
                            beta_X3trt * (treatment * X3)
    )
    linpred_censor <- with(studyData,
                           base_risk_cen +
                             beta_trt_cen * treatment +
                             beta_X1_cen * X1 +
                             beta_X2_cen * X2 +
                             beta_X3_cen * X3
    )
    
    eventTime     <- (-log(runif(nStudy)) / (lambda_event * exp(linpred_event)))^(1/alpha_event)
    censoringTime <- (-log(runif(nStudy)) / (lambda_censor * exp(linpred_censor)))^(1/alpha_censor)
    obsTime       <- pmin(eventTime, censoringTime)
    eventObserved <- as.integer(eventTime <= censoringTime)
    censor        <- 1 - eventObserved
    
    studyData_All[[i]] <- data.frame(
      studyData,
      eventTime      = eventTime,
      censoringTime  = censoringTime,
      obsTime        = obsTime,
      eventObserved  = eventObserved,
      censor         = censor
    )
    
    # -- clean Target IPD -- #
    tp_true <- data.frame(
      ID        = 1:nTargetPopulation,
      X1        = rbinom(nTargetPopulation, 1, 0.35),
      X2        = rbinom(nTargetPopulation, 1, 0.55),
      X3        = rnorm(nTargetPopulation, 0.21, 1.5)
    )
    tp_true$treatment <- rbinom(nTargetPopulation, 1, 0.5)
    
    linpred_t <- with(tp_true,
                      beta_X1_trt * X1 +
                        beta_X2_trt * X2 +
                        beta_X3_trt * X3 +
                        beta_treatment * treatment +
                        beta_X1trt * (treatment * X1) +
                        beta_X2trt * (treatment * X2) +
                        beta_X3trt * (treatment * X3)
    )
    eventTime_target <- (-log(runif(nTargetPopulation)) / (lambda_event * exp(linpred_t)))^(1/alpha_event)
    tp_true$eventTime <- eventTime_target
    IPDTargetPopulationData_true[[i]] <- tp_true
    
    # abnormal values in target IPD #
    # for various abv pattern:
    tp_obs <- tp_true
    
    if (scen$ABV_Type %in% c("X2_flip", "both")) {
      n2 <- floor(scen$X2_Flips_Prop * nTargetPopulation)
      if (n2 > 0) {
        idx2 <- sample(nTargetPopulation, n2)
        tp_obs$X2[idx2] <- 1 - tp_obs$X2[idx2]
      }
    }
    
    if (scen$ABV_Type %in% c("X3_extreme", "both")) {
      n3 <- floor(scen$X3_Extremes_Prop * nTargetPopulation)
      if (n3 > 0) {
        idx3 <- sample(nTargetPopulation, n3)
        sd3  <- sd(tp_obs$X3)
        tp_obs$X3[idx3] <- tp_obs$X3[idx3] + scen$X3_Extremes_Multiplier * sd3
      }
    }
    
    IPDTargetPopulationData_obs[[i]] <- tp_obs
    
    # clean target AgD #
    samp_true <- tp_true[sample(1:nTargetPopulation, nTarget), ]
    aggregateTargetData_true[[i]] <- data.frame(
      N       = nTarget,
      X1_prop = mean(samp_true$X1),
      X2_prop = mean(samp_true$X2),
      X3_mean = mean(samp_true$X3),
      X3_sd   = sd(samp_true$X3)
    )
    
    # abv target AgD #
    samp_obs <- tp_obs[sample(1:nTargetPopulation, nTarget), ]
    aggregateTargetData_obs[[i]] <- data.frame(
      N       = nTarget,
      X1_prop = mean(samp_obs$X1),
      X2_prop = mean(samp_obs$X2),
      X3_mean = mean(samp_obs$X3),
      X3_sd   = sd(samp_obs$X3)
    )
  }
  
  return(list(
    studyData                     = studyData_All,
    
    aggregateTargetData_abv       = aggregateTargetData_obs,
    aggregateTargetData_true      = aggregateTargetData_true,
    
    IPDTargetPopulationData_abv   = IPDTargetPopulationData_obs,
    IPDTargetPopulationData_true  = IPDTargetPopulationData_true
  ))
}

