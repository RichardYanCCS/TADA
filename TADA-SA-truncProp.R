rm(list = ls())
options(mc.cores = parallel::detectCores()) 

library(survival)        
library(tidyr)
library(dplyr)           
library(TransportHealth)
library(parallel)
library(doParallel)
library(foreach)

source("./assistFunctions.R")

rr <- 1            # scenario
n_study <- 200     # study sample size

cenProp <- 0.2     # study sample censoring proportion 
if (cenProp == 0.2) {
  base_risk_cen <- 2.5
} else if (cenProp == 0.3) {
  base_risk_cen <- 3.3
} else if (cenProp == 0.4) {
  base_risk_cen <- 3.7
} else if (cenProp == 0.5) {
  base_risk_cen <- 4.3
} else {
  stop("Censoring Proportion must be one of 0.2, 0.3, 0.4, or 0.5.")
}

nRep <- 500
BSRep <- 200


# sensitivity analysis: final weights truncation percentile 
cutq_final <- c(0.80, 0.90, 0.95, 0.99, NA)

summary_list <- list()

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

start_time <- Sys.time()

for (cutq in cutq_final) {
  
  print(paste0("[ Sensitivity Analysis - Truncation ] The final weights truncation setting is: ", cutq))
  
  trunc_obj <- if (is.na(cutq)) NULL else TransportHealth::trunc("quantile", cutq)
  
  HR_cen          <- numeric(nRep)
  HR_cen_SE       <- numeric(nRep)
  CI_lower_cen    <- numeric(nRep)
  CI_upper_cen    <- numeric(nRep)
  
  HR_noncen       <- numeric(nRep)
  HR_noncen_SE    <- numeric(nRep)
  CI_lower_noncen <- numeric(nRep)
  CI_upper_noncen <- numeric(nRep)
  
  trueHR          <- numeric(nRep)
  trueHR_SE       <- numeric(nRep)
  CI_lower_true   <- numeric(nRep)
  CI_upper_true   <- numeric(nRep)
  
  for (rep in seq_len(nRep)) {
    
    print(paste0("---- ", rep, "-th replicate starts ----"))

    testData <- generateTestDataTADA_TTE(nStudy = n_study, base_risk_cen = base_risk_cen)
    
    studyData    <- testData$studyData[[rr]][, -1]                # study sample IPD
    rownames(studyData) <- seq_len(nrow(studyData))
    aggregateData <- testData$aggregateTargetData[[rr]]           # target sample AgD 
    targetIPD <- testData$IPDTargetPopulationData[[rr]][,-1]      # target population IPD
    
    # --------- HR estimated by MoM participation weights only ------------- #
    
    # HR point estimate
    max_MoM_attempts <- 50       
    MoM_attempt      <- 0       
    
    repeat {
      MoM_attempt <- MoM_attempt + 1
      if (MoM_attempt > max_MoM_attempts) {
        stop(sprintf(
          "scenario = %d, replicate = %d: MoM only weights estimation has attempted %d times but failed all. Stop here.",
          rr, rep, max_MoM_attempts
        ))
      }
      
      results_MOM_only <- transportTADA(
        msmFormula           = Surv(obsTime, eventObserved) ~ treatment,
        
        propensityWeights    = rep(1, nrow(studyData)), # i.e., final weights = 1 * MoM weights
        participationWeights = NULL,                    # MoM weights estimated by transportTADA function
        
        matchingCovariates   = c("X1","X2","X3"),
        exOpt                = list(
          propensity    = NULL,
          participation = NULL,
          final         = trunc_obj        # since final weights, apply truncation rule
        ),
        family               = "coxph",
        studyData            = studyData,
        aggregateTargetData  = aggregateData
      )
      
      if (any(!is.finite(results_MOM_only$finalWeights))) {
        message(sprintf(
          "scenario = %d, replicate = %d: Attempt %d: non-finite MoM weights detected (Inf/-Inf), retrying...",
          rr, rep, MoM_attempt
        ))
        next
      }
      
      if (any(is.na(results_MOM_only$finalWeights))) {
        message(sprintf(
          "scenario = %d, replicate = %d: Attempt %d: NA value found in MoM weights, retrying...",
          rr, rep, MoM_attempt
        ))
        next
      }
      
      if (any(results_MOM_only$finalWeights <= 0)) {
        message(sprintf(
          "scenario = %d, replicate = %d: Attempt %d: non-positive MoM weights detected (<= 0), retrying...",
          rr, rep, MoM_attempt
        ))
        next
      }
      
      HR_noncen[rep] <- exp(results_MOM_only$msm$coefficients)
      break
    }
    
    
    # Bootstrap: MoM weights == final weights (parallel)
    HR_noncen_Boot <- foreach(
      b = 1:BSRep,
      .combine        = c,
      .packages       = c("survival", "dplyr", "tidyr", "TransportHealth"),
      .errorhandling = "pass"
    ) %dopar% {
      tryCatch({
        max_retry <- 50
        retry     <- 0
        hr_b      <- NA
        
        repeat {
          retry <- retry + 1
          if (retry > max_retry) {
            hr_b <- NA
            break
          }
          
          studyData_noncen_Boot <- studyData[sample.int(n = nrow(studyData), size = nrow(studyData), replace = TRUE), ]
          
          bootstrap_part_result <- transportTADA(
            msmFormula           = Surv(obsTime, eventObserved) ~ treatment,
            
            propensityWeights    = rep(1, nrow(studyData_noncen_Boot)),
            participationWeights = NULL,
            
            matchingCovariates   = c("X1", "X2", "X3"),
            exOpt                = list(
              propensity    = NULL,
              participation = NULL,
              final         = trunc_obj
            ),
            family               = "coxph",
            studyData            = studyData_noncen_Boot,
            aggregateTargetData  = aggregateData
          )
          
          fw <- bootstrap_part_result$finalWeights
          
          if (any(!is.finite(fw)) || any(is.na(fw)) || any(fw <= 0)) {
            message(sprintf(
              "scenario = %d, replicate = %d: The %d -th attempt detects invalid finalWeights (Inf/NA/<=0), retry...",
              rr, rep, retry
            ))
            next
          }
          
          hr_b <- exp(bootstrap_part_result$msm$coefficients)
          break
        }
        hr_b
        
      }, error = function(e) {NA} ) # end tryCatch
    } # end dopar
    
    
    HR_noncen_SE[rep]   <- sd(HR_noncen_Boot, na.rm = TRUE)
    CI_lower_noncen[rep] <- quantile(HR_noncen_Boot, 0.025, na.rm = TRUE)
    CI_upper_noncen[rep] <- quantile(HR_noncen_Boot, 0.975, na.rm = TRUE)
    
    # --- Final Weights == MoM Weights ends here ------------- #
    
    
    
    
    
    # --- time-varying IPCW ---------------------------------- #
    
    max_rep_attempts <- 50  
    rep_attempt      <- 0  
    
    repeat {
      rep_attempt <- rep_attempt + 1
      
      if (rep_attempt > max_rep_attempts) {
        stop(sprintf(
          "scenario = %d, replicate = %d: more than %d attempts due to invalid final weights, stop here.",
          rr, rep, max_rep_attempts
        ))
      }
      
      # time‐varying censoring weights
      all_times <- sort(unique(studyData$obsTime))
      unique_event_times <- all_times[ all_times < max(all_times) ]    # ensure interval length > 0 (avoid the last interval is (max_time, max_time) )
      
      
      # split study data into long format, break into intervals by event time
      studyData_long_event <- survSplit(
        studyData,
        cut   = unique_event_times,
        start = "Tstart",
        end   = "obsTime",
        event = "eventObserved",             # event
        id    = "ID"
      ) %>%
        filter(obsTime > Tstart) %>%      # ensure interval length > 0
        arrange(ID, obsTime)
      
      studyData_long_censor <- survSplit(
        studyData,
        cut   = unique_event_times,
        start = "Tstart",
        end   = "obsTime",
        event = "censor",                   # censor
        id    = "ID"
      ) %>%
        filter(obsTime > Tstart) %>%
        arrange(ID, obsTime)
      
      studyData_long_event$censor <- studyData_long_censor$censor
      studyData_long_event$Tstop <- studyData_long_censor$obsTime        # Tstop = obsTime
      
      studyData_long_event <- studyData_long_event %>% 
        filter(Tstop > Tstart) %>% 
        select(-c(eventTime, censoringTime, obsTime)) %>% 
        select(c(ID, X1, X2, X3, treatment, eventObserved, censor, Tstart, Tstop))
      
      # fit time-varying censoring model: risk of censoring ~ treatment + covariates
      censoring_cox_model <- coxph(
        Surv(Tstart, Tstop, censor) ~ treatment + X1 + X2 + X3,
        data    = studyData_long_event,
        timefix = FALSE
      )
      
      # Baseline cumulative hazards & calculate linear predictor
      baseline_hazard <- basehaz(censoring_cox_model, centered = FALSE)
      
      studyData_long_event <- studyData_long_event %>%
        mutate(
          linear_predictor = as.matrix(select(., treatment, X1, X2, X3)) %*% coef(censoring_cox_model)
        ) %>%
        left_join(baseline_hazard, by = c("Tstop" = "time")) %>% # studyData_long_event$Tstop == baseline_hazard$time
        arrange(ID, Tstop) %>%
        mutate(
          survival_function_censor = exp(- hazard * exp(linear_predictor)),
          time_varying_weights    = 1 / survival_function_censor
        ) 
      
      
      # convert time-fixed participation weights into long format and calculate the time-varying final weights
      # consider truncation
      part_wt_for_final <- transportTADA(
        msmFormula           = Surv(obsTime, eventObserved) ~ treatment,
        
        propensityWeights    = rep(1, nrow(studyData)),
        participationWeights = NULL,
        
        matchingCovariates   = c("X1","X2","X3"),
        exOpt                = list(propensity    = NULL,
                                    participation = NULL,
                                    final         = NULL), # no truncation here, truncate the final weights
        family               = "coxph",
        studyData            = studyData,
        aggregateTargetData  = aggregateData
      )
      participation_weights <- part_wt_for_final$participationWeights
      
      studyData_long_event <- studyData_long_event %>%
        left_join(
          data.frame(ID = seq_len(nrow(studyData)), part_wt = participation_weights),
          by = "ID"
        ) 
      
      if (!is.na(cutq)) {
        studyData_long_event <- studyData_long_event %>%
          mutate(raw_final_wt = time_varying_weights * part_wt) %>%
          { 
            thr <- quantile(.$raw_final_wt, probs = cutq, na.rm = TRUE)
            mutate(., final_weights = pmax(pmin(raw_final_wt, thr), 1e-6))
          }
      } else {
        studyData_long_event <- studyData_long_event %>%
          mutate(
            raw_final_wt   = time_varying_weights * part_wt,
            final_weights  = pmax(raw_final_wt, 1e-6)
          )
      }
      
      
      
      # ------ HR estimated by final weights = MoM weights * censoring weights ------------------------------------------------------------------------------- #
      
      # point estimate
      min_wt   <- min(studyData_long_event$final_weights, na.rm = TRUE)
      count_le0 <- sum(studyData_long_event$final_weights <= 0, na.rm = TRUE)
      
      if (count_le0 > 0) {
        message(sprintf(
          "scenario = %d, replicate = %d: a total of %d final weights ≤ 0 (min = %g), retry this replicate...",
          rr, rep, count_le0, min_wt
        ))
        next
      }
      
      weighted_cox_model_result <- survival::coxph(
        Surv(Tstart, Tstop, eventObserved) ~ treatment,
        data   = studyData_long_event,
        weight = studyData_long_event$final_weights
      )
      
      HR_cen[rep] <- exp(weighted_cox_model_result$coefficients)
      break
    }
    
    
    # Bootstrap for MoM weights * censoring weights
    HR_cen_Boot <- foreach(
      b = 1:BSRep,
      .combine  = c,
      .packages = c("survival","dplyr","tidyr","TransportHealth"),
      .errorhandling = "pass"
    ) %dopar% {
      set.seed(20250425 + rr * 1e5 + rep * 100 + b)
      
      tryCatch({
        max_retry <- 50
        retry     <- 0
        hr_b      <- NA
        
        repeat {
          retry <- retry + 1
          if (retry > max_retry) {
            hr_b <- NA  
            break
          }
          
          studyData_cen_Boot <- studyData[sample.int(n = nrow(studyData), size = nrow(studyData), replace = TRUE), ]
          rownames(studyData_cen_Boot) <- seq_len(nrow(studyData_cen_Boot))
          
          bootstrap_part_wt_for_final <- transportTADA(
            msmFormula           = Surv(obsTime, eventObserved) ~ treatment,
            
            propensityWeights    = rep(1, nrow(studyData_cen_Boot)),
            participationWeights = NULL,
            
            matchingCovariates   = c("X1","X2","X3"),
            exOpt                = list(propensity    = NULL,
                                        participation = NULL,
                                        final         = NULL),
            family               = "coxph",
            studyData            = studyData_cen_Boot,
            aggregateTargetData  = aggregateData
          )
          part_wt_for_final_Boot <- bootstrap_part_wt_for_final$participationWeights
          
          all_times_boot <- sort(unique(studyData_cen_Boot$obsTime))
          unique_event_times_Boot <- all_times_boot[ all_times_boot < max(all_times_boot) ]
          
          
          studyData_long_event_boot <- survSplit(
            studyData_cen_Boot,
            cut   = unique_event_times_Boot,
            start = "Tstart", end = "obsTime",
            event = "eventObserved", id = "ID"
          ) %>%
            filter(obsTime > Tstart) %>%
            arrange(ID, obsTime)
          
          studyData_long_censor_boot <- survSplit(
            studyData_cen_Boot,
            cut   = unique_event_times_Boot,
            start = "Tstart", end   = "obsTime",
            event = "censor", id = "ID"
          ) %>%
            filter(obsTime > Tstart) %>%
            arrange(ID, obsTime)
          
          studyData_long_event_boot$censor <- studyData_long_censor_boot$censor
          studyData_long_event_boot$Tstop  <- studyData_long_censor_boot$obsTime
          
          if (sum(studyData_long_event_boot$censor) == 0 ||
              sum(studyData_long_event_boot$eventObserved) == 0) {
            next
          }
          
          studyData_long_event_boot <- studyData_long_event_boot %>% 
            filter(Tstop > Tstart) %>% 
            select(-c(eventTime, censoringTime, obsTime)) %>% 
            select(c(ID, X1, X2, X3, treatment, eventObserved, censor, Tstart, Tstop))
          
          # fit time-varying censoring model
          censoring_cox_model_bootstrap <- coxph(
            Surv(Tstart, Tstop, censor) ~ treatment + X1 + X2 + X3,
            data    = studyData_long_event_boot,
            timefix = FALSE
          )
          beta_boot <- coef(censoring_cox_model_bootstrap)
          
          if (any(is.na(beta_boot)) || 
              any(!is.finite(beta_boot))
          ) {
            next
          }
          baseline_hazard_bootstrap <- basehaz(censoring_cox_model_bootstrap, centered = FALSE)
          
          studyData_long_event_boot <- studyData_long_event_boot %>%
            mutate(linear_predictor = as.matrix(select(., treatment, X1, X2, X3)) %*% beta_boot)
          if (any(is.na(studyData_long_event_boot$linear_predictor)) || 
              any(!is.finite(studyData_long_event_boot$linear_predictor))) {
            next
          }
          
          studyData_long_event_boot <- studyData_long_event_boot %>%
            left_join(baseline_hazard_bootstrap, by = c("Tstop" = "time")) %>%
            arrange(ID, Tstop) 

          if (any(is.na(studyData_long_event_boot$hazard)) ||
              any(!is.finite(studyData_long_event_boot$hazard)) ) {
            next
          }
          
          studyData_long_event_boot <- studyData_long_event_boot %>%
            mutate(censoring_survival   = pmax(exp(-hazard * exp(linear_predictor)), 1e-6))
          if (any(is.na(studyData_long_event_boot$baseline_survival)) ||
              any(is.na(studyData_long_event_boot$censoring_survival)) ||
              any(studyData_long_event_boot$censoring_survival <= 0) ||
              any(!is.finite(studyData_long_event_boot$baseline_survival)) ||
              any(!is.finite(studyData_long_event_boot$censoring_survival)) ) {
            next
          }
          
          studyData_long_event_boot <- studyData_long_event_boot %>%
            mutate(time_varying_weights_bootstrap = pmin(1 / censoring_survival, 1e6))
          if (any(is.na(studyData_long_event_boot$time_varying_weights_bootstrap)) ||
              any(studyData_long_event_boot$time_varying_weights_bootstrap <= 0) ||
              any(!is.finite(studyData_long_event_boot$time_varying_weights_bootstrap)) ) {
            next
          }
          
          studyData_long_event_boot <- studyData_long_event_boot %>%
            left_join(
              data.frame(ID = seq_len(nrow(studyData_cen_Boot)),
                         part_wt = part_wt_for_final_Boot),
              by = "ID"
            ) %>%
            mutate(raw_final_wt_bootstrap = time_varying_weights_bootstrap * part_wt)
          if (any(is.na(studyData_long_event_boot$raw_final_wt_bootstrap)) ||
              any(studyData_long_event_boot$raw_final_wt_bootstrap <= 0) ||
              any(!is.finite(studyData_long_event_boot$raw_final_wt_bootstrap))){
            next
          }
          
          if (!is.na(cutq)) {
            studyData_long_event_boot <- studyData_long_event_boot %>%
              { 
                thr_boot <- quantile(.$raw_final_wt_bootstrap, probs = cutq, na.rm = TRUE)
                mutate(., final_weights_bootstrap = pmax(pmin(raw_final_wt_bootstrap, thr_boot), 1e-6))
              }
          } else {
            studyData_long_event_boot <- studyData_long_event_boot %>%
              mutate(final_weights_bootstrap = pmax(raw_final_wt_bootstrap, 1e-6))
          }
          if (any(is.na(studyData_long_event_boot$final_weights_bootstrap)) ||
              any(!is.finite(studyData_long_event_boot$final_weights_bootstrap)) ||
              any(studyData_long_event_boot$final_weights_bootstrap <= 0)) {
            next
          }
          
          wc <- studyData_long_event_boot$final_weights_bootstrap
          
          fit <- coxph(Surv(Tstart, Tstop, eventObserved) ~ treatment,
                       data = studyData_long_event_boot,
                       weight = wc)
          
          res <- exp(coef(fit))
          break
          
        } # end repeat
        
        res
        
      }  # end foreach
      , error = function(e) {
        return(NA)  
      })
    }
    
    HR_cen_SE[rep]   <- sd(HR_cen_Boot, na.rm = TRUE)
    CI_lower_cen[rep] <- quantile(HR_cen_Boot, 0.025, na.rm = TRUE)
    CI_upper_cen[rep] <- quantile(HR_cen_Boot, 0.975, na.rm = TRUE)
    
    
    
    # ------  true model  --------------------------------------------------------------------------------- #
    
    res_true <- coxph(Surv(eventTime) ~ treatment, data = targetIPD, robust = T)
    
    logHR <- coef(res_true)["treatment"]  
    logSE <- sqrt(res_true$var[1,1]) 
    
    trueHR[rep] <- exp(logHR)
    trueHR_SE[rep] <- trueHR[rep] * logSE 
    
    CI_lower_true[rep] <- exp(logHR - 1.96*logSE)
    CI_upper_true[rep] <- exp(logHR + 1.96*logSE)
    
    # ----------------------------------------------------------------------------------------------------- #
    

  } # end rep loop
  

  
  
  # Average effect of all replicates in each scenario
  {
    # MoM weights * IPCW
    HR_cen_average <- mean(HR_cen, na.rm = T)
    HR_cen_SE_average <- mean(HR_cen_SE, na.rm = T)
    
    bias_HR_cen <- mean(HR_cen - trueHR, na.rm = T)
    coverage_HR_cen <-  sum(trueHR >= CI_lower_cen & trueHR <= CI_upper_cen, na.rm = T) / nRep
    
    CI_lower_cen_average <- mean(CI_lower_cen) # HR CI
    CI_upper_cen_average <- mean(CI_upper_cen)
    
    Percentile_Lower_Bias_HR_censor <- quantile(HR_cen - trueHR, 0.025) # bias CI
    Percentile_Upper_Bias_HR_censor <- quantile(HR_cen - trueHR, 0.975)
    
    # MoM weights only
    HR_noncen_average <- mean(HR_noncen, na.rm = T)
    HR_noncen_SE_average <- mean(HR_noncen_SE, na.rm = T)
    
    bias_HR_noncen <- mean(HR_noncen - trueHR, na.rm = T)
    coverage_HR_noncen <-  sum(trueHR >= CI_lower_noncen & trueHR <= CI_upper_noncen, na.rm = T) / nRep
    
    CI_lower_noncen_average <- mean(CI_lower_noncen) # HR CI
    CI_upper_noncen_average <- mean(CI_upper_noncen)
    
    Percentile_Lower_Bias_HR_noncensor <- quantile(HR_noncen - trueHR, 0.025) # bias CI
    Percentile_Upper_Bias_HR_noncensor <- quantile(HR_noncen - trueHR, 0.975)
    
    # True
    trueHR_average <- mean(trueHR, na.rm = T)
    trueHR_SE_average <- mean(trueHR_SE, na.rm = T)
    CI_lower_true_average <- mean(CI_lower_true)
    CI_upper_true_average <- mean(CI_upper_true)
    
    # Source of bias
    MonteCarlo_SE_cen <- sd(HR_cen)
    Bootstrap_SE_cen <- mean(HR_cen_SE)
    
    MonteCarlo_SE_noncen <- sd(HR_noncen)
    Bootstrap_SE_noncen <- mean(HR_noncen_SE)
  }
  
  raw_data_list <- list(
    HR_cen            = HR_cen,
    HR_cen_SE         = HR_cen_SE,
    CI_lower_cen      = CI_lower_cen,
    CI_upper_cen      = CI_upper_cen,
    HR_noncen         = HR_noncen,
    HR_noncen_SE      = HR_noncen_SE,
    CI_lower_noncen   = CI_lower_noncen,
    CI_upper_noncen   = CI_upper_noncen,
    trueHR            = trueHR,
    trueHR_SE         = trueHR_SE,
    CI_lower_true     = CI_lower_true,
    CI_upper_true     = CI_upper_true
  )
  
  summary_list[[length(summary_list) + 1]] <- data.frame(
    Scenario           = rr,
    N_study            = n_study,
    Trunc_cutoff       = if (is.na(cutq)) "No Truncation" else paste0(cutq*100, "%"),
    HR_censor          = sprintf("%.3f (%.3f, %.3f)",
                                 HR_cen_average, CI_lower_cen_average, CI_upper_cen_average),
    HR_noncensor       = sprintf("%.3f (%.3f, %.3f)",
                                 HR_noncen_average, CI_lower_noncen_average, CI_upper_noncen_average),
    HR_true            = sprintf("%.3f (%.3f, %.3f)",
                                 trueHR_average, CI_lower_true_average, CI_upper_true_average),
    Bias_censor        = sprintf("%.3f (%.3f, %.3f)",
                                 bias_HR_cen,
                                 Percentile_Lower_Bias_HR_censor,
                                 Percentile_Upper_Bias_HR_censor),
    Bias_noncensor     = sprintf("%.3f (%.3f, %.3f)",
                                 bias_HR_noncen,
                                 Percentile_Lower_Bias_HR_noncensor,
                                 Percentile_Upper_Bias_HR_noncensor),
    Coverage_censor    = sprintf("%.2f%%", 100 * coverage_HR_cen),
    Coverage_noncensor = sprintf("%.2f%%", 100 * coverage_HR_noncen),
    
    MonteCarlo_SE_cen  = sprintf("%.2f%%", MonteCarlo_SE_cen),
    Bootstrap_SE_cen   = sprintf("%.2f%%", Bootstrap_SE_cen),
    
    MonteCarlo_SE_noncen = sprintf("%.2f%%", MonteCarlo_SE_noncen),
    Bootstrap_SE_noncen  = sprintf("%.2f%%", Bootstrap_SE_noncen),
    
    stringsAsFactors   = FALSE
  )

}

stopCluster(cl)


end_time <- Sys.time()

execution_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
print(paste("Simulation Total Time:", round(execution_time, 2), "mins."))

summary_table <- do.call(rbind, summary_list)


