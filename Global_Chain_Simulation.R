# library(dplyr)
# library(readr)
# library(lubridate)
# library(ggplot2)
# library(nloptr)
#
# setwd("~/UCLA PhD/SCALE Group/Research/Kidney Exchange/Global Chains Paper")
# # # To avoid having to rerun everything:
# load("data.RData")

# KIDPAN_DATA <- read.csv("KIDPAN_DATA_FULL.csv")
# #Donors <- read.csv("Donor_Data.csv")
#
# KIDPAN_DATA <- KIDPAN_DATA %>% filter(DAYSWAIT_CHRON_KI > 0, DON_TY %in% c("L", "F"))
#
#
# LIVING_DONORS <- read_delim("LIVING_DONOR_DATA.DAT",
#                             delim = "\t", escape_double = FALSE,
#                             col_names = FALSE, trim_ws = TRUE)
#
# colnames(LIVING_DONORS) <- c("REOP_BILIARY", "REOP_BILIARY_DT", "DON_DATE", "AGE_DON", "ETHCAT_DON", "REGION", "LUNG_RECOV",
#                              "KIDNEY_RECOV", "LIVER_RECOV", "DON_ORG2", "CITIZENSHIP", "LIV_DON_TY", "LIV_DON_TY_OSTXT",
#                              "EDUCATION", "KI_CREAT_PREOP", "BP_PREOP_SYST", "BP_PREOP_DIAST", "COD", "KI_PROC_TY", "MARITAL_STAT",
#                              "HEALTH_INS", "FUNC_STAT", "PHYSICAL_CAPACITY", "WORK_INCOME", "HIST_CANCER", "HIST_CANCER_OSTXT",
#                              "CANCER_FREE", "COD_OSTXT", "HIST_HYPER", "HYPER_DIET", "HYPER_DIUR", "HYPER_MEDS", "PREOP_URINE_RATIO",
#                              "DIABETES", "MACRO_FAT", "MICRO_FAT", "CONVERT_OPEN_KI", "NON_AUTO_BLOOD", "PRBC_UNITS", "PLATELETS_UNITS",
#                              "FFP_UNITS", "VASC_COMP_KI", "VASC_COMP_KI_INTER", "VASC_COMP_KI_INTER_OSTXT", "OTH_COMP_KI",
#                              "OTH_COMP_KI_INTER", "OTH_COMP_KI_INTER_OSTXT", "REOPERATION_KI", "REOP_BLEED_KI", "REOP_HERNIA_KI",
#                              "REOP_BOWEL_KI", "REOP_VASC_KI", "REOP_OTH_KI", "REOP_OTH_KI_OSTXT", "READMISSION_KI",
#                              "READMISSION_KI_REASON", "READMISSION_KI_OSTXT", "OTH_INTER_PROC_KI", "OTH_INTER_PROC_KI_OSTXT",
#                              "KI_CREAT_POSTOP", "BP_POSTOP_SYST", "BP_POSTOP_DIAST", "HYPERTENSION", "POSTOP_URINE_RATIO",
#                              "PREOP_BILI", "PREOP_SGOT_AST", "PREOP_ALK_PHOS", "PREOP_ALBUM", "PREOP_INR", "BIOPSY_LI",
#                              "LI_PROC_TY", "BILIARY_COMP", "BILIARY_COMP_GRADE", "VASC_COMP_LI", "VASC_COMP_LI_INTER",
#                              "VASC_COMP_LI_INTER_OSTXT", "OTH_COMP_LI", "OTH_COMP_LI_INTER", "OTH_COMP_LI_INTER_OSTXT",
#                              "REOPERATION_LI", "REOP_BLEED_LI", "REOP_HERNIA_LI", "REOP_BOWEL_LI", "REOP_VASC_LI", "REOP_OTH_LI",
#                              "REOP_OTH_LI_OSTXT", "REOP_LI_FAIL", "READMISSION_LI", "READMISSION_LI_REASON", "READMISSION_LI_OSTXT",
#                              "OTH_INTER_PROC_LI", "OTH_INTER_PROC_LI_OSTXT", "POSTOP_SGOT_AST", "POSTOP_ALK_PHOS", "POSTOP_ALBUM",
#                              "POSTOP_CREAT_LI", "POSTOP_INR", "POSTOP_SGPT_ALT", "POSTOP_BILI", "PREOP_FVC_BEFORE", "PREOP_FVC_AFTER",
#                              "PREOP_FEV1_BEFORE", "PREOP_FEV1_AFTER", "PREOP_FEF_BEFORE", "PREOP_FEF_AFTER", "PREOP_TLC_BEFORE",
#                              "PREOP_TLC_AFTER", "PREOP_LUNG_CAP", "PREOP_PAO2", "HIST_CIG", "PACK_YRS", "DUR_ABSTINENCE",
#                              "TOBACCO_USE", "LU_PROC_TY", "CONVERT_OPEN_LU", "INTRAOP_COMP", "INTRAOP_COMP_REASON",
#                              "SACRIFICE_LOBE", "ARRHYTHMIA", "ANESTHETIC_COMP", "INTRAOP_COMP_OSTXT", "LU_COMP", "LU_COMP_REASON",
#                              "LU_COMP_OSTXT", "THORAC_TUBES", "ARRHYTHMIA_POSTOP", "READMISSION_LU", "READMISSION_LU_REASON",
#                              "READMISSION_LU_OSTXT", "PREDON_HGT", "PREDON_WGT", "PREOP_SGPT_ALT", "PREOP_CREAT_LI",
#                              "PREOP_URINE_PROTEIN", "POSTOP_URINE_PROTEIN", "PX_STAT", "CITIZEN_COUNTRY", "DEATH_DT",
#                              "ORG_RECOVERY_DT", "INIT_DISCHARGE_DT", "REOP_BLEED_KI_DT", "REOP_HERNIA_KI_DT", "REOP_BOWEL_KI_DT",
#                              "REOP_VASC_KI_DT", "REOP_OTH_KI_DT", "READMISSION_KI_DT", "OTH_INTER_PROC_KI_DT", "POSTOP_TEST_DT",
#                              "REOP_BLEED_LI_DT", "REOP_HERNIA_LI_DT", "REOP_BOWEL_LI_DT", "REOP_VASC_LI_DT", "REOP_OTH_LI_DT",
#                              "REOP_LI_FAIL_DT", "READMISSION_LI_DT", "OTH_INTER_PROC_LI_DT", "READMISSION_LU_DT", "CMV_IGG",
#                              "CMV_IGM", "CMV_NUCLEIC", "EBV_IGG", "EBV_IGM", "HBV_CORE", "HBV_DNA", "HBV_SUR_ANTIGEN",
#                              "HCV_ANTIBODY", "HCV_RIBA", "HCV_RNA", "VIRUSES_TESTED", "CMV_TOTAL", "EBV_TOTAL", "HOME_STATE",
#                              "GENDER", "ABO", "YR_ENTRY_US", "WGT_KG", "DON_ORG", "STATUS_LDR", "VAL_DT_LDR", "DBW4", "DBW6",
#                              "DC1", "DC2", "DDP1", "DDP2", "DDPA1", "DDPA2", "DDR51", "DDR51_2", "DDR52", "DDR52_2", "DDR53",
#                              "DDR53_2", "DDQ1", "DDQ2", "DDQA1", "DDQA2", "RECOV_FACILITY_CODE", "DONOR_ID")
#
#
# LIVING_DONORS <- LIVING_DONORS %>% filter(KIDNEY_RECOV == 1)
#
# ## Record ETHCAT to list Race titles
# KIDPAN_DATA <- KIDPAN_DATA %>%
#   mutate(ETHCAT = case_when(
#     ETHCAT == 1 ~ "White",
#     ETHCAT == 2 ~ "Black",
#     ETHCAT == 4 ~ "Hispanic",
#     ETHCAT == 5 ~ "Asian",
#     ETHCAT == 6 ~ "Other",
#     ETHCAT == 7 ~ "Other",
#     ETHCAT == 9 ~ "Other",
#     ETHCAT == 998 ~ "Other",
#     TRUE ~ as.character(ETHCAT)
#   ))
#
# ## Record ETHCAT_DON to list Race titles
# KIDPAN_DATA <- KIDPAN_DATA %>%
#   mutate(ETHCAT_DON = case_when(
#     ETHCAT_DON == 1 ~ "White",
#     ETHCAT_DON == 2 ~ "Black",
#     ETHCAT_DON == 4 ~ "Hispanic",
#     ETHCAT_DON == 5 ~ "Asian",
#     ETHCAT_DON == 6 ~ "Other",
#     ETHCAT_DON == 7 ~ "Other",
#     ETHCAT_DON == 9 ~ "Other",
#     ETHCAT_DON == 998 ~ "Other",
#     TRUE ~ as.character(ETHCAT_DON)
#   ))
#
#
# LIVING_DONORS <- LIVING_DONORS %>%
#   mutate(ETHCAT_DON = case_when(
#     ETHCAT_DON == 1 ~ "White",
#     ETHCAT_DON == 2 ~ "Black",
#     ETHCAT_DON == 4 ~ "Hispanic",
#     ETHCAT_DON == 5 ~ "Asian",
#     ETHCAT_DON == 6 ~ "Other",
#     ETHCAT_DON == 7 ~ "Other",
#     ETHCAT_DON == 9 ~ "Other",
#     ETHCAT_DON == 998 ~ "Other",
#     TRUE ~ as.character(ETHCAT_DON)
#   ))
#
# # Convert INIT_DATE to Date format
# KIDPAN_DATA$INIT_DATE <- as.Date(KIDPAN_DATA$INIT_DATE, format = "%m/%d/%Y")
# LIVING_DONORS$DON_DATE <- as.Date(LIVING_DONORS$DON_DATE, format = "%m/%d/%Y")
#
# # Filter the dataset to include only dates later than first donor observation
# KIDPAN_DATA <- subset(KIDPAN_DATA, INIT_DATE > min(LIVING_DONORS$DON_DATE))
#
# ########################################################################################
# Patient_Data <- KIDPAN_DATA[, c("PT_CODE", "ETHCAT","DAYSWAIT_CHRON_KI", "INIT_DATE", "END_DATE", "DONOR_ID")]
# Donor_Data <- LIVING_DONORS[, c("ETHCAT_DON", "DONOR_ID", "DON_DATE")]
#
#
# Patient_Data <- Patient_Data %>% filter(ETHCAT %in% c("Asian", "Black", "Hispanic", "White"))
# Donor_Data <- Donor_Data %>% filter(ETHCAT_DON %in% c("Asian", "Black", "Hispanic", "White"), DONOR_ID %in% Patient_Data$DONOR_ID)
#
#
# Patient_Data$END_DATE <- as.Date(Patient_Data$END_DATE, format = "%m/%d/%Y")
# Patient_Data$INIT_DATE <- as.Date(Patient_Data$INIT_DATE, format = "%m/%d/%Y")
#
#
# Donor_Data <- Donor_Data %>% filter(DON_DATE > "2000-01-01")
# Patient_Data <- Patient_Data %>% filter(INIT_DATE > "2000-01-01")
# Required libraries
library(tidyverse)
library(lubridate)
library(ggplot2)
library(R6)

###########################################
# Parameter Tracker Class
###########################################

TemporalParameterTracker <- R6::R6Class("TemporalParameterTracker",
                                        public = list(
                                          # Storage for smoothed values for each group and metric
                                          smoothed_values = list(),
                                          
                                          # Single smoothing parameter
                                          alpha = NULL,  # smoothing factor
                                          
                                          # Current best estimates of zeta and xi
                                          current_zeta = 1.0,
                                          current_xi = 1.0,
                                          
                                          initialize = function(alpha = 0.03) {
                                            self$alpha <- alpha
                                          },
                                          
                                          init_smoothing_if_needed = function(group, metric, initial_value) {
                                            if (is.null(self$smoothed_values[[group]])) {
                                              self$smoothed_values[[group]] <- list()
                                            }
                                            
                                            if (is.null(self$smoothed_values[[group]][[metric]])) {
                                              self$smoothed_values[[group]][[metric]] <- initial_value
                                            }
                                          },
                                          
                                          # Apply simple exponential smoothing for a single observation Y_t
                                          exponential_update = function(group, metric, Y_t) {
                                            current_value <- self$smoothed_values[[group]][[metric]]
                                            
                                            # Simple exponential smoothing update
                                            S_t <- self$alpha * Y_t + (1 - self$alpha) * current_value
                                            
                                            # Handle numerical stability
                                            S_t <- if(is.finite(S_t)) S_t else current_value
                                            
                                            self$smoothed_values[[group]][[metric]] <- S_t
                                            return(S_t)
                                          },
                                          
                                          apply_smoothing = function(current_stats) {
                                            # For each group, apply simple exponential smoothing to T_current and Q_current
                                            for (i in seq_len(nrow(current_stats))) {
                                              group <- current_stats$ETHCAT[i]
                                              
                                              # Initialize if needed
                                              self$init_smoothing_if_needed(group, "T_current", current_stats$T_current[i])
                                              self$init_smoothing_if_needed(group, "Q_current", current_stats$Q_current[i])
                                              
                                              # Update smoothed values
                                              T_smoothed <- self$exponential_update(group, "T_current", current_stats$T_current[i])
                                              Q_smoothed <- self$exponential_update(group, "Q_current", current_stats$Q_current[i])
                                              
                                              # Replace original with smoothed
                                              current_stats$T_current[i] <- T_smoothed
                                              current_stats$Q_current[i] <- Q_smoothed
                                            }
                                            return(current_stats)
                                          },
                                          
                                          update_and_smooth = function(group_stats, current_arrival_rates) {
                                            # Apply smoothing
                                            smoothed_stats <- self$apply_smoothing(group_stats)
                                            
                                            # Compute parameters using smoothed statistics
                                            params <- self$compute_parameters(smoothed_stats, current_arrival_rates)
                                            
                                            self$current_zeta <- params$zeta
                                            self$current_xi <- params$xi
                                            
                                            return(list(
                                              zeta = self$current_zeta,
                                              xi = self$current_xi
                                            ))
                                          },
                                          
                                          compute_parameters = function(smoothed_stats, current_arrival_rates) {
                                            Q_min <- min(smoothed_stats$Q_current)
                                            W_max <- max(smoothed_stats$T_current)
                                            W_min <- min(smoothed_stats$T_current)
                                            delta_wait <- W_max - W_min
                                            delta_Q <- diff(range(smoothed_stats$Q_current))
                                            
                                            # Modified zeta computation to more strongly prioritize wait time differences
                                            zeta_raw <- if (delta_Q > 0 && delta_wait > 0) {
                                              # Base scaling on relative wait time difference
                                              wait_ratio <- W_max / W_min
                                              base_zeta <- (wait_ratio * Q_min) / delta_Q
                                              
                                              # Add progressive scaling for larger wait time differences
                                              scaling_factor <- log1p(delta_wait / W_min)
                                              base_zeta * scaling_factor
                                            } else {
                                              1.0
                                            }
                                            
                                            # Enforce minimum zeta to ensure sufficient priority differentiation
                                            zeta <- zeta_raw
                                            
                                            # Compute B* with modified formula to increase sensitivity to wait time differences
                                            B_star <- delta_wait * log1p(zeta * Q_min) - 
                                              (W_min * zeta * delta_Q) / (1 + zeta * Q_min)
                                            
                                            # Modified xi computation for stronger preference toward high-wait groups
                                            lambda_max <- max(current_arrival_rates$rates)
                                            gamma <- 1
                                            c_effective <- mean(smoothed_stats$c)
                                            
                                            xi_base <- if (abs(B_star) < 1e-9 || (c_effective - 1) <= 0) {
                                              2.0  # Higher default value
                                            } else {
                                              # Increase sensitivity to wait time differences
                                              base_xi <- (2 / abs(B_star)) * log1p(lambda_max * max(smoothed_stats$Q_current) /
                                                                                     (gamma * Q_min * (c_effective - 1)))
                                              # Scale xi based on relative wait time difference
                                              scaling <- log1p(wait_ratio)
                                              base_xi * scaling
                                            }
                                            
                                            # Ensure xi is large enough to create meaningful probability differences
                                            xi <- xi_base
                                            
                                            return(list(zeta = zeta, xi = xi))
                                          }
                                        )
)

###########################################
# Helper Functions
###########################################

select_patient_from_group <- function(waiting_patients, selected_group) {
  group_patients <- waiting_patients %>%
    filter(ETHCAT == selected_group)
  
  if(nrow(group_patients) == 0) return(NULL)
  
  # Find the patient with the highest wait_time
  max_wait <- max(group_patients$wait_time, na.rm = TRUE)
  eligible_patients <- group_patients %>%
    filter(wait_time == max_wait)
  
  # If multiple patients have the same max wait, pick one at random
  if(nrow(eligible_patients) > 1) {
    selected_patient <- eligible_patients[sample(nrow(eligible_patients), 1), ]
  } else {
    selected_patient <- eligible_patients
  }
  
  return(selected_patient)
}


compute_priority <- function(T_current, mu_current, domestic_rate, Q_current,
                             verbose = TRUE, zeta) {
  mu_current <- as.numeric(mu_current)
  T_current <- as.numeric(T_current)
  Q_current <- as.numeric(Q_current)
  
  theta <- T_current * log(1 + zeta * Q_current)
  
  theta[is.nan(theta) | is.infinite(theta)] <- max(theta[!is.nan(theta) & !is.infinite(theta)],
                                                   na.rm = TRUE)
  
  if (verbose) {
    cat("\nPriority Score Computation:\n")
    cat("Current Waitlist Avg:", T_current, "\n")
    cat("Queue Size:", Q_current, "\n")
    cat("Patient arrival rate:", mu_current, "\n")
    cat("Domestic match rate:", domestic_rate, "\n")
    cat("Combined Priority Score:", round(theta, 4), "\n")
  }
  
  return(theta)
}

compute_prob <- function(theta, xi) {
  theta_centered <- theta - mean(theta)
  
  exp_theta <- exp(theta_centered * xi)
  
  exp_theta[is.nan(exp_theta) | is.infinite(exp_theta)] <- max(exp_theta[!is.nan(exp_theta) & !is.infinite(exp_theta)], na.rm = TRUE)
  exp_theta <- pmax(exp_theta, 1e-10)
  
  probs <- exp_theta / sum(exp_theta)
  
  epsilon <- 1e-6
  probs <- (probs + epsilon) / sum(probs + epsilon)
  
  if(any(is.na(probs)) || sum(probs) < 0.99 || sum(probs) > 1.01) {
    probs <- rep(1/length(theta), length(theta))
  }
  
  cat("Softmax Probability:", round(probs, 4), "\n")
  return(probs)
}

###########################################
# Data Processing Functions
###########################################

process_initial_data <- function(Patient_Data, Donor_Data, foreign_donor_rate = 0.5) {
  # Process dates
  Patient_Data$END_DATE <- as.Date(Patient_Data$END_DATE, format = "%m/%d/%Y")
  Patient_Data$INIT_DATE <- as.Date(Patient_Data$INIT_DATE, format = "%m/%d/%Y")
  Donor_Data$DON_DATE <- as.Date(Donor_Data$DON_DATE, format = "%m/%d/%Y")
  
  # Filter data
  Patient_Data <- Patient_Data %>%
    filter(ETHCAT %in% c("Asian", "Black", "Hispanic", "White"),
           INIT_DATE > "1990-01-01")
  
  Donor_Data <- Donor_Data %>%
    filter(ETHCAT_DON %in% c("Asian", "Black", "Hispanic", "White"),
           DONOR_ID %in% Patient_Data$DONOR_ID,
           DON_DATE > "1990-01-01")
  
  # Calculate time window and generate foreign donors
  total_days <- as.numeric(max(Donor_Data$DON_DATE) - min(Donor_Data$DON_DATE))
  set.seed(123)
  num_events <- rpois(1, foreign_donor_rate * total_days)
  arrival_times <- sort(runif(num_events, 0, total_days))
  foreign_dates <- min(Donor_Data$DON_DATE) + arrival_times
  
  foreign_donors <- data.frame()
  for(arrival_date in foreign_dates) {
    template_donor <- Donor_Data[sample(nrow(Donor_Data), 1), ]
    new_donor <- template_donor
    new_donor$DON_DATE <- arrival_date
    new_donor$DONOR_ID <- paste0("F", sample(100000:999999, 1))
    new_donor$is_foreign <- TRUE
    foreign_donors <- rbind(foreign_donors, new_donor)
  }
  
  Donor_Data$is_foreign <- FALSE
  Donor_Data <- rbind(Donor_Data, foreign_donors) %>% arrange(DON_DATE)
  
  cat("Total time period:", total_days, "days\n")
  cat("Expected foreign donors:", foreign_donor_rate * total_days, "\n")
  cat("Actual foreign donors:", sum(Donor_Data$is_foreign), "\n")
  cat("Average foreign donors per day:", sum(Donor_Data$is_foreign) / total_days, "\n")
  cat("Original donors:", nrow(Donor_Data) - sum(Donor_Data$is_foreign), "\n")
  
  return(list(Patient_Data = Patient_Data, Donor_Data = Donor_Data))
}

compute_monthly_arrival_rates <- function(Patient_Data, current_date) {
  month_start <- lubridate::floor_date(current_date, "month")
  month_end <- lubridate::ceiling_date(current_date, "month") - lubridate::days(1)
  
  monthly_arrivals <- Patient_Data %>%
    filter(INIT_DATE >= month_start,
           INIT_DATE <= month_end) %>%
    group_by(ETHCAT) %>%
    summarize(
      arrivals = n()
    ) %>%
    mutate(
      # Convert to daily rate based on the number of days in the month
      days_in_month = as.numeric(month_end - month_start + 1),
      mu_current = arrivals / days_in_month
    )
  
  return(list(
    rates = monthly_arrivals$mu_current,
    groups = monthly_arrivals$ETHCAT
  ))
}


compute_monthly_domestic_rates <- function(Patient_Data, Donor_Data, current_date) {
  month_start <- lubridate::floor_date(current_date, "month")
  month_end <- lubridate::ceiling_date(current_date, "month") - lubridate::days(1)
  
  monthly_matches <- Patient_Data %>%
    filter(END_DATE >= month_start,
           END_DATE <= month_end,
           status == "matched",
           !(DONOR_ID %in% Donor_Data$DONOR_ID[Donor_Data$is_foreign])) %>%
    group_by(ETHCAT) %>%
    summarize(
      matches = n()
    ) %>%
    mutate(
      days_in_month = as.numeric(month_end - month_start + 1),
      domestic_rate = matches / days_in_month
    )
  
  all_groups <- data.frame(
    ETHCAT = unique(Patient_Data$ETHCAT)
  ) %>%
    left_join(monthly_matches, by = "ETHCAT") %>%
    mutate(
      domestic_rate = ifelse(is.na(domestic_rate), 0, domestic_rate)
    )
  
  return(list(
    rates = all_groups$domestic_rate,
    groups = all_groups$ETHCAT
  ))
}


###########################################
# Main Simulation Function
###########################################

run_simulation <- function(Patient_Data, Donor_Data, start_date = 1500) {
  parameter_tracker <- TemporalParameterTracker$new(alpha = 0.028)
  
  # Initialize patient data
  Patient_Data$wait_time <- 0
  Patient_Data$status <- "waiting"
  
  simulation_days <- as.numeric(max(Donor_Data$DON_DATE, na.rm = TRUE) -
                                  min(Patient_Data$INIT_DATE, na.rm = TRUE))
  current_date <- min(Patient_Data$INIT_DATE, na.rm = TRUE)
  
  tracking_data <- data.frame(
    day = integer(),
    ETHCAT = character(),
    avg_wait_time = numeric(),
    queue_length = integer(),
    arrival_rate = numeric(),
    domestic_rate = numeric()
  )
  
  current_arrival_rates <- compute_monthly_arrival_rates(Patient_Data, current_date)
  current_domestic_rates <- compute_monthly_domestic_rates(Patient_Data, Donor_Data, current_date)
  
  
  
  
  for(day in 1:simulation_days) {
    if(day %% 100 == 0) cat("\nProcessing day:", day, "\n")
    current_date <- current_date + days(1)
    
    # Update rates monthly
    if(day == 1 || month(current_date) != month(current_date - days(1))) {
      current_arrival_rates <- compute_monthly_arrival_rates(Patient_Data, current_date)
      current_domestic_rates <- compute_monthly_domestic_rates(Patient_Data, Donor_Data, current_date)
    }
    
    
    # Update wait times
    Patient_Data$wait_time[Patient_Data$status == "waiting" &
                             Patient_Data$INIT_DATE <= current_date] <-
      Patient_Data$wait_time[Patient_Data$status == "waiting" &
                               Patient_Data$INIT_DATE <= current_date] + 1
    
    # Process foreign donors
    if(day >= start_date) {
      foreign_donors_today <- Donor_Data %>%
        filter(is_foreign == TRUE, DON_DATE == current_date)
      
      for(i in 1:nrow(foreign_donors_today)) {
        waiting_patients <- Patient_Data %>%
          filter(status == "waiting",
                 INIT_DATE <= current_date,
                 END_DATE >= current_date)
        
        if(nrow(waiting_patients) > 0) {
          group_stats <- waiting_patients %>%
            group_by(ETHCAT) %>%
            summarize(
              T_current = mean(wait_time, na.rm = TRUE),
              Q_current = n(),
              T_max = max(wait_time, na.rm = TRUE),
              c = T_max / T_current
            )
          
          # Add current rates
          group_stats$mu_current <- current_arrival_rates$rates[match(group_stats$ETHCAT, current_arrival_rates$groups)]
          group_stats$domestic_rate <- current_domestic_rates$rates[match(group_stats$ETHCAT, current_domestic_rates$groups)]
          
          # Update parameters using Holt smoothing
          smoothed_params <- parameter_tracker$update_and_smooth(group_stats, current_arrival_rates)
          
          # Compute priorities and probabilities
          theta <- compute_priority(
            group_stats$T_current,
            group_stats$mu_current,
            group_stats$domestic_rate,
            group_stats$Q_current,
            zeta = smoothed_params$zeta
          )
          
          probs <- compute_prob(theta, smoothed_params$xi)
          
          if(day %% 100 == 0) {
            cat("\nDay", day, "matching statistics:\n")
            print(data.frame(
              ETHCAT = group_stats$ETHCAT,
              Wait_Time = group_stats$T_current,
              Queue_Length = group_stats$Q_current,
              Priority_Score = theta,
              Selection_Prob = probs
            ))
          }
          
          # Select group and patient
          selected_group <- sample(group_stats$ETHCAT, 1, prob = probs)
          matched_patient <- select_patient_from_group(waiting_patients, selected_group)
          
          if(!is.null(matched_patient)) {
            Patient_Data$status[Patient_Data$PT_CODE == matched_patient$PT_CODE] <- "matched"
          }
        }
      }
    }
    
    # Process regular matches and removals
    Patient_Data <- Patient_Data %>%
      mutate(
        status = case_when(
          status == "waiting" & END_DATE == current_date &
            !(DONOR_ID %in% Donor_Data$DONOR_ID[Donor_Data$is_foreign]) ~ "matched",
          status == "waiting" & END_DATE == current_date ~ "removed",
          TRUE ~ status
        )
      )
    
    # Record statistics weekly
    if(day %% 7 == 0) {
      current_waiting_patients <- Patient_Data %>%
        filter(status == "waiting",
               INIT_DATE <= current_date,
               END_DATE >= current_date)
      
      current_stats <- current_waiting_patients %>%
        group_by(ETHCAT) %>%
        summarize(
          avg_wait_time = mean(wait_time, na.rm = TRUE),
          queue_length = n()
        )
      
      current_stats$arrival_rate <- current_arrival_rates$rates[match(current_stats$ETHCAT,
                                                                      current_arrival_rates$groups)]
      current_stats$domestic_rate <- current_domestic_rates$rates[match(current_stats$ETHCAT,
                                                                        current_domestic_rates$groups)]
      
      tracking_data <- rbind(tracking_data,
                             data.frame(
                               day = rep(day, nrow(current_stats)),
                               current_stats
                             ))
    }
  }
  
  final_stats <- Patient_Data %>%
    group_by(ETHCAT) %>%
    summarize(
      avg_wait_time = mean(wait_time[status %in% c("matched", "removed")]),
      matched_rate = sum(status == "matched") / n(),
      removed_rate = sum(status == "removed") / n()
    )
  
  return(list(
    tracking_data = tracking_data,
    final_stats = final_stats,
    final_patient_data = Patient_Data
  ))
}

###########################################
# Visualization Functions
###########################################

create_visualizations <- function(results, start_date) {
  wait_time_plot <- results$tracking_data %>%
    filter(day <= 8000) %>%
    ggplot(aes(x = day, y = avg_wait_time, color = ETHCAT)) +
    geom_line() +
    geom_vline(xintercept = start_date, linetype = "dashed", color = "red") +
    labs(title = "Average Wait Times by Ethnic Category Over Time",
         x = "Simulation Day",
         y = "Average Wait Time (days)",
         color = "Ethnic Category") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")
  
  convergence_plot <- results$tracking_data %>%
    group_by(day) %>%
    summarize(
      max_diff = max(avg_wait_time) - min(avg_wait_time)
    ) %>%
    filter(day <= 8000) %>%
    ggplot(aes(x = day)) +
    geom_line(aes(y = max_diff), color = "red") +
    labs(title = "Convergence Over Time",
         x = "Simulation Day",
         y = "Max Waiting Time Difference") +
    theme_minimal()
  
  return(list(
    wait_time_plot = wait_time_plot,
    convergence_plot = convergence_plot
  ))
}

###########################################
# Run Simulation
###########################################

# Make sure you have Patient_Data and Donor_Data loaded.
# Example:
# Patient_Data <- ... # load patient data
# Donor_Data <- ...   # load donor data

start_date <- 1500
set.seed(1)
processed_data <- process_initial_data(Patient_Data, Donor_Data)
results <- run_simulation(processed_data$Patient_Data, processed_data$Donor_Data, start_date)
plots <- create_visualizations(results, start_date)
print(plots$wait_time_plot)
#print(plots$convergence_plot)
