rm(list=ls())

library(TransportHealth)
library(survival)
library(ggplot2)
library(survminer)
library(grid)
library(tidyr)
library(dplyr)
library(tidyverse)
library(scales)

set.seed(20250428)

IPD <- read.csv("./processed_source_IPD.csv")
AgD <- read.csv("./processed_target_AgD.csv")

print(paste("Overall Censoring Proportion in Study Data: ", round(mean(IPD$CENS == 1), 2)))

knitr::kable(
  data.frame(
    Source     = c("IPD", "AgD"),
    Mean_Age   = c(mean(IPD$AGE),              AgD$AGE_MEAN),
    Sex_Prop   = c(mean(IPD$SEX == "Male"),    AgD$SEX_PROP)
  ),
  col.names = c("", "AGE mean", "SEX male(%)"),
  digits    = 3,
  format    = "markdown"
)

# ----- Time Varying IPCW Estimation ----- #

all_times <- sort(unique(IPD$OS))
unique_event_times <- all_times[ all_times < max(all_times) ]


# split study data into long format, break into intervals by event time
IPD_long_event <- survSplit(
  IPD,
  cut   = unique_event_times,
  start = "Tstart",
  end   = "OS",
  event = "DEATH",            
  id    = "ID"
) %>%
  filter(OS > Tstart) %>%      
  arrange(ID, OS)

IPD_long_censor <- survSplit(
  IPD,
  cut   = unique_event_times,
  start = "Tstart",
  end   = "OS",
  event = "CENS",                   
  id    = "ID"
) %>%
  filter(OS > Tstart) %>%
  arrange(ID, OS)

IPD_long_event$CENS <- IPD_long_censor$CENS
IPD_long_event$Tstop <- IPD_long_censor$OS        # Tstop = OS

IPD_long_event <- IPD_long_event %>% 
  filter(Tstop > Tstart) %>% 
  select(-c(OS, EVENT_REASON)) %>% 
  select(c(ID, AGE, SEX, RACE, SMOKE, ECOGBL, CENS, Tstart, Tstop, DEATH))

# fit time-varying censoring model: risk of censoring ~ covariates
censoring_cox_model <- coxph(
  Surv(Tstart, Tstop, CENS) ~  AGE + SEX + RACE + SMOKE + ECOGBL,
  data    = IPD_long_event,
  timefix = FALSE
)

# cox PH assumption test: global p = 0.43 > 0.05 cannot reject null assumption (PH satisfied)
plot(cox.zph(censoring_cox_model)) # alongside 0 smoothly

# Baseline cumulative hazards & calculate linear predictor
baseline_hazard <- basehaz(censoring_cox_model, centered = FALSE)

IPD_long_event <- IPD_long_event %>%
  mutate(linear_predictor = predict(censoring_cox_model, newdata = ., type    = "lp")) %>%
  left_join(baseline_hazard, by = c("Tstop" = "time")) %>% 
  arrange(ID, Tstop) %>%
  mutate(
    censoring_survival   = exp(- hazard * exp(linear_predictor)),
    time_varying_weights = 1 / censoring_survival
  )

# convert time-fixed participation weights into long format and calculate the time-varying final weights
# consider truncation
IPD$SEX <- ifelse(IPD$SEX == "Male", 1, 0) # SEX 1 = Male 0 = Female, switch to 0/1

part_wt_for_final <- transportTADA(
  msmFormula           = Surv(OS, DEATH) ~ 1,
  
  propensityWeights    = rep(1, nrow(IPD)),
  participationWeights = NULL,
  
  matchingCovariates   = c("AGE", "SEX"),
  exOpt                = list(propensity    = NULL,
                              participation = NULL,
                              final         = NULL),
  family               = "coxph",
  studyData            = IPD,
  aggregateTargetData  = AgD
)
participation_weights <- part_wt_for_final$participationWeights

part_wt_df <- as.data.frame(participation_weights)
rownames(part_wt_df) <- seq_len(nrow(part_wt_df))
colnames(part_wt_df) <- c("MoM_wt")


IPD_long_event <- IPD_long_event %>%
  left_join( data.frame(ID = seq_len(nrow(IPD)), part_wt = participation_weights), by = "ID") %>%
  mutate(raw_final_wt  = time_varying_weights * part_wt)


cutq_TADA <- 0.95
cutq_MoM <- 0.95
cutq_cen <- 0.95

IPD_long_event <- IPD_long_event %>%
  mutate(final_weights = pmax(pmin(raw_final_wt, 
                                   quantile(raw_final_wt, probs = cutq_TADA, na.rm = TRUE)), 1e-6))


weighted_cox_model_result <- survival::coxph(
  Surv(Tstart, Tstop, DEATH) ~ 1,          
  data   = IPD_long_event,
  weight = IPD_long_event$final_weights,
  timefix = F
)

trunc_obj <- TransportHealth::trunc("quantile", cutq_MoM)

results_MOM_only <- transportTADA(
  msmFormula           = Surv(OS, DEATH) ~ 1,
  propensityWeights    = rep(1, nrow(IPD)), # i.e., final weights = 1 * MoM weights
  participationWeights = NULL,             
  matchingCovariates   = c("AGE","SEX"),
  exOpt                = list(propensity = NULL, participation = NULL, final = trunc_obj),
  family               = "coxph",
  studyData            = IPD,
  aggregateTargetData  = AgD
)


# Generate weighted KM curves after various weighting adjustment
km_fit_origin1 <- survfit(Surv(Tstart, Tstop, DEATH) ~ 1, data = IPD_long_event, 
                          weights = pmin(IPD_long_event$time_varying_weights,
                                         quantile(IPD_long_event$time_varying_weights, probs = cutq_cen, na.rm  = TRUE)),
                          conf.int = T)     # with IPCW adjustment only

km_fit_adj <- survfit(Surv(Tstart, Tstop, DEATH) ~ 1, data = IPD_long_event, 
                      weights = IPD_long_event$final_weights,
                      conf.int = T)         # TADA full adjustment


km_fit_origin2 <- survfit(Surv(OS, DEATH) ~ 1, data = IPD,
                          conf.int = T)     # without any adjustment

km_fit_wo_adj <- survfit(Surv(OS, DEATH) ~ 1, data = IPD, 
                         weights = results_MOM_only$finalWeights,
                         conf.int = T)      # MoM adjustment only

# Convert survfit objects to data frames using surv_summary
orig_df_1   <- surv_summary(km_fit_origin1, data = IPD_long_event)    # Original KM curves data with censoring weights
orig_df_2   <- surv_summary(km_fit_origin2, data = IPD)               # Original KM curves 
adj_df    <- surv_summary(km_fit_adj, data = IPD_long_event)          # Transported and Censoring adjusted KM curves data
wo_adj_df <- surv_summary(km_fit_wo_adj, data = IPD)                  # Transported but Censoring unadjusted KM curves data



# overall survival rates at given time points
time_points <- c(6, 12, 18, 24, 30, 36)  # months

surv_rates_origin1 <- summary(km_fit_origin1, times = time_points)$surv
surv_rates_origin2 <- summary(km_fit_origin2, times = time_points)$surv
surv_rates_adj     <- summary(km_fit_adj, times = time_points)$surv
surv_rates_wo_adj  <- summary(km_fit_wo_adj, times = time_points)$surv

knitr::kable(
  data.frame(
    Month     = time_points,
    Original  = round(surv_rates_origin2, 3),
    MoM_only  = round(surv_rates_wo_adj,  3),
    IPCW_only = round(surv_rates_origin1, 3),
    TADA_full = round(surv_rates_adj,     3)
  ),
  col.names = c("Month", "Original", "Covariate-Adj", "Censor-Adj", "TADA (Full)"),
  align     = c("c","c","c","c","c"),
  format    = "markdown",
  digits    = 3
)


# Median OS Rate and CI's
extract_medianOS_CI <- function(km_fit) {
  median_table <- summary(km_fit)$table
  data.frame(
    Median_OS = median_table["median"],
    Lower_95CI = median_table["0.95LCL"],
    Upper_95CI = median_table["0.95UCL"]
  )
}

MedianOS_origin1_CI <- extract_medianOS_CI(km_fit_origin1)
MedianOS_origin2_CI <- extract_medianOS_CI(km_fit_origin2)
MedianOS_adj_CI     <- extract_medianOS_CI(km_fit_adj)
MedianOS_wo_adj_CI  <- extract_medianOS_CI(km_fit_wo_adj)

median_CI_results <- data.frame(
  Model = c("Original", "Covariates imbalance adjusted only", "Censoring adjusted only", "TADA fully adjusted"),
  
  Median_OS = c(round(MedianOS_origin2_CI$Median_OS,3),  # without any adjustment
                round(MedianOS_wo_adj_CI$Median_OS,3),   # MoM adjustment only
                
                round(MedianOS_origin1_CI$Median_OS,3),  # with IPCW adjustment only
                round(MedianOS_adj_CI$Median_OS, 3)),    # full TADA adjustment
  
  Lower_95CI = c(round(MedianOS_origin2_CI$Lower_95CI, 3),
                 round(MedianOS_wo_adj_CI$Lower_95CI, 3),
                 
                 round(MedianOS_origin1_CI$Lower_95CI, 3),
                 round(MedianOS_adj_CI$Lower_95CI, 3)),
  
  Upper_95CI = c(round(MedianOS_origin2_CI$Upper_95CI, 3),
                 round(MedianOS_wo_adj_CI$Upper_95CI, 3),
                 
                 round(MedianOS_origin1_CI$Upper_95CI, 3),
                 round(MedianOS_adj_CI$Upper_95CI, 3))
)

median_CI_results

# IPCW only checking
median_os_cen <- round(MedianOS_origin1_CI$Median_OS, 3)
ggplot(IPD, aes(x = OS, fill = factor(CENS))) +
  geom_histogram(
    bins     = 160,
    position = "identity",
    alpha    = 0.6
  ) +
  geom_vline(
    xintercept = median_os_cen,
    linetype   = "dashed",
    color      = "red",
    size       = 1.5
  ) +
  annotate(
    "text",
    x     = median_os_cen,
    y     = 11,                    
    label = paste0("Median OS = ", median_os_cen),
    color = "red",
    size  = 7,
    hjust = -0.1                   
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "0" = "#4E5138",  
      "1" = "#06436F"   
    ),
    labels = c("Event", "Censor")
  ) +
  scale_x_continuous(
    name   = "Months",
    breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36),
    limits = c(0, 36)
  ) +
  scale_y_continuous(
    name   = "Count",
    breaks = seq(0, 12, 3),
    limits = c(0, 12)
  ) +
  theme_minimal() +
  theme(
    text            = element_text(family = "Helvetica", size = 24),
    legend.position = "top",
    axis.title.x    = element_text(size = 24),
    axis.title.y    = element_text(size = 24),
    axis.text       = element_text(size = 24),
    legend.title    = element_text(size = 24),
    legend.text     = element_text(size = 24)
  )


# KM Curves
#   df_orig_2   : KM summary for the unadjusted ("Original") curve
#   df_wo_adj   : KM summary for the MoM-only ("Covariate-Adj") curve
#   df_orig_1   : KM summary for the IPCW-only ("Censoring-Adj") curve
#   df_adj      : KM summary for the full TADA-adjusted ("TADA-Full") curve
df_orig_2 <- orig_df_2
df_wo_adj <- wo_adj_df
df_orig_1 <- orig_df_1
df_adj <- adj_df

df_orig_2$Group <- "Original Unadjusted" #: No any adjustments"
df_wo_adj$Group <- "Covariates Adjusted" #: Only MoM and No Censoring Adjusted"
df_orig_1$Group <- "Censoring Adjusted"  #: With Censoring Weights"
df_adj$Group    <- "TADA"                #: Censoring and MoM Adjusted"

df_combined_1 <- rbind(df_orig_2, df_wo_adj)  
df_combined_2 <- rbind(df_orig_1, df_adj)     

# -----------------------

median_OS_1 <- summary(km_fit_origin2)$table["median"]
median_OS_2 <- summary(km_fit_wo_adj)$table["median"]
median_OS_3 <- summary(km_fit_origin1)$table["median"]
median_OS_4 <- summary(km_fit_adj)$table["median"]


levels_all <- c("Original",  
                "Transported (w/o Censoring Adjusted)", 
                "Original (with Censoring Adjusted)", 
                "TADA (Transported with Censoring Adjusted)")

df_combined_cov  <- bind_rows(
  df_orig_2 %>% mutate(Group = "Original"),        
  df_wo_adj %>% mutate(Group = "Transported (w/o Censoring Adjusted)")    
) %>% mutate(Group = factor(Group, levels = levels_all))

df_combined_cen  <- bind_rows(
  df_orig_2 %>% mutate(Group = "Original"),
  df_orig_1 %>% mutate(Group = "Original (with Censoring Adjusted)")    
) %>% mutate(Group = factor(Group, levels = levels_all))

df_combined_full <- bind_rows(
  df_orig_2 %>% mutate(Group = "Original"),
  df_adj    %>% mutate(Group = "TADA (Transported with Censoring Adjusted)")        
) %>% mutate(Group = factor(Group, levels = levels_all))

med_orig <- as.numeric(median_OS_1)  # Original unadjusted median
med_cov  <- as.numeric(median_OS_2)  # MoM-only median
med_cen  <- as.numeric(median_OS_3)  # IPCW-only median
med_full <- as.numeric(median_OS_4)  # Full TADA-adjusted median

curve_colors <- c(
  `Original` = "grey50",
  `Transported (w/o Censoring Adjusted)`= "orange",
  `Original (with Censoring Adjusted)`= "blue",
  `TADA (Transported with Censoring Adjusted)`    = "purple"
)

p_cov <- ggplot(df_combined_cov, aes(x = time, y = surv, 
                                     color = Group, group = Group)) +
  geom_step(size = 1) +
  geom_vline(xintercept = med_orig, color = "grey50", linetype = "dashed", size = 1) +
  geom_vline(xintercept = med_cov,  color = "red", linetype = "dashed", size = 1) +
  scale_color_manual(values = curve_colors, name = NULL) +
  scale_x_continuous(breaks = c(0,3,6,9,12,15,18,21,24,27,30,33,36),
                     limits = c(0, 36)) +  
  labs(
    x     = "Months",
    y     = "Survival Probability"
  ) +
  annotate(
    "text",
    x     = median_OS_2,
    y     = 0.3,                    
    label = paste0("Median OS = ", round(median_OS_2,3)),
    color = "red",
    size  = 7, 
    hjust = 1.1                   
  ) +
  annotate(
    "text",
    x     = median_OS_1,
    y     = 0.7,                    
    label = paste0("Median OS = ", round(median_OS_1,3)),
    color = "grey50", fontface = "bold",
    size  = 7,
    hjust = -0.1                   
  ) +
  theme_minimal() +
  theme(
    text            = element_text(family = "Helvetica", size = 24),
    legend.position = "top",
    plot.title      = element_text(size = 24, hjust = 0.5),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 24),
    legend.text     = element_text(size = 24)
  )

p_full <- ggplot(df_combined_full, aes(x = time, y = surv, 
                                       color = Group, group = Group)) +
  geom_step(size = 1) +
  geom_vline(xintercept = med_orig, color = "grey50", linetype = "dashed", size = 1) +
  geom_vline(xintercept = med_full, color = "red", linetype = "dashed", size = 1) +
  scale_color_manual(values = curve_colors, name = NULL) +
  scale_x_continuous(breaks = c(0,3,6,9,12,15,18,21,24,27,30,33,36),
                     limits = c(0, 36)) + 
  labs(
    x     = "Months",
    y     = "Survival Probability"
  ) + 
  annotate(
    "text",
    x     = median_OS_4,
    y     = 0.7,                    
    label = paste0("Median OS = ", round(median_OS_4,3)),
    color = "red",
    size  = 7, 
    hjust = -0.1                   
  ) + 
  annotate(
    "text",
    x     = median_OS_1,
    y     = 0.3,                    
    label = paste0("Median OS = ", round(median_OS_1,3)),
    color = "grey50", fontface = "bold",
    size  = 7,
    hjust = 1.1                   
  ) +
  theme_minimal() +
  theme(
    text            = element_text(family = "Helvetica", size = 24),
    legend.position = "top",
    plot.title      = element_text(size = 24, hjust = 0.5),
    axis.title      = element_text(size = 24),
    axis.text       = element_text(size = 24),
    legend.text     = element_text(size = 24)
  )

gridExtra::grid.arrange(
  p_cov,
  p_full,
  ncol = 1
)


knitr::kable(
  data.frame(
    Month     = time_points,
    Original  = round(surv_rates_origin2, 3),
    MoM_only  = round(surv_rates_wo_adj,  3),
    IPCW_only = round(surv_rates_origin1, 3),
    TADA_full = round(surv_rates_adj,     3)
  ),
  col.names = c("Month", "Original", "Covariate-Adj", "Censor-Adj", "TADA"),
  align     = c("c","c","c","c","c"),
  format    = "markdown",
  digits    = 3
)



# Archive 

# quantile_raw <- quantile(IPD_long_event$raw_final_wt,
#                          probs = c(0.90, 0.95, 0.99),
#                          na.rm = TRUE)
# 
# quantile_cen <- quantile(IPD_long_event$time_varying_weights,
#                          probs = c(0.90, 0.95, 0.99),
#                          na.rm = TRUE)
# 
# quantile_MoM <- quantile(part_wt_df$MoM_wt,
#                          probs = c(0.90, 0.95, 0.99),
#                          na.rm = TRUE)
# 
# 
# trunc_judge_TADA <- ggplot(IPD_long_event, aes(x = raw_final_wt)) +
#   geom_histogram(
#     aes(y = after_stat(count / sum(count))), bins  = 80, fill  = "steelblue", color = "white") +
#   geom_vline(xintercept = quantile_raw, linetype = "dashed", color = "red", size = 0.7) +
#   annotate( "text",
#     x = quantile_raw[1], y = 0.18, label = "90%",
#     angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_raw[2], y = 0.15, label = "95%",
#     angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_raw[3], y = 0.13, label = "99%", 
#     angle = 0, hjust = -0.3, vjust = 3.5, color = "red", size = 7) +
#   scale_x_continuous(limits = c(0, 4), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)) +
#   scale_y_continuous(limits = c(0, 0.25), labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Raw Final Weights",  y = "Percentage (%)") +
#   theme_minimal() +
#   theme(
#     # axis.title.x = element_text(hjust = 0.22),
#     axis.title = element_text(family = "Helvetica", size = 20),
#     axis.text  = element_text(family = "Helvetica", size = 20)
#   )
# 
# trunc_judge_MoM <- ggplot(part_wt_df, aes(x = MoM_wt)) +
#   geom_histogram(
#     aes(y = after_stat(count / sum(count))), bins  = 80, fill  = "steelblue", color = "white") +
#   geom_vline(xintercept = quantile_MoM, linetype = "dashed", color = "red", size = 0.7) +
#   annotate( "text",
#             x = quantile_MoM[1], y = 0.18, label = "90%",
#             angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_MoM[2], y = 0.15, label = "95%",
#     angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_MoM[3], y = 0.13, label = "99%", 
#     angle = 0, hjust = -0.3, vjust = 3.5, color = "red", size = 7) +
#   scale_x_continuous(limits = c(0, 4), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5,4)) +
#   scale_y_continuous(limits = c(0, 0.25), labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Raw MoM Weights",  y = "Percentage (%)") +
#   theme_minimal() +
#   theme(
#     # axis.title.x = element_text(hjust = 0.22),
#     axis.title = element_text(family = "Helvetica", size = 20),
#     axis.text  = element_text(family = "Helvetica", size = 20)
#   )
# 
# trunc_judge_cen <- ggplot(IPD_long_event, aes(x = time_varying_weights)) +
#   geom_histogram(
#     aes(y = after_stat(count / sum(count))), bins  = 80, fill  = "steelblue", color = "white") +
#   geom_vline(xintercept = quantile_cen, linetype = "dashed", color = "red", size = 0.7) +
#   annotate( "text", x = quantile_cen[1], y = 0.18, label = "90%",
#             angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_cen[2], y = 0.15, label = "95%",
#     angle = 0, hjust = -0.2, vjust = 3.5, color = "red", size  = 7) +
#   annotate(
#     "text", x = quantile_cen[3], y = 0.13, label = "99%", 
#     angle = 0, hjust = -0.3, vjust = 3.5, color = "red", size = 7) +
#   scale_x_continuous(limits = c(1, 2)) +
#   scale_y_continuous(limits = c(0, 0.25), labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Raw Censoring Weights",  y = "Percentage (%)") +
#   theme_minimal() +
#   theme(
#     # axis.title.x = element_text(hjust = 0.22),
#     axis.title = element_text(family = "Helvetica", size = 20),
#     axis.text  = element_text(family = "Helvetica", size = 20)
#   )
# 
# 
# gridExtra::grid.arrange(
#   trunc_judge_TADA, trunc_judge_MoM, trunc_judge_cen,
#   ncol = 1
# )


# Since the study cohort had a very low rate of censoring until approximately 10 months ago, 
# the IPCW weight distribution was concentrated at 1. 
# Therefore, application of IPCW alone did not significantly change the Kaplan-Meier median survival estimate (9.86 months, 95 % CI 8.90-11.14), 
# The effect of IPCW weighting on survival in the tail of the curve (>20 months) was more pronounced as censoring continued to occur.