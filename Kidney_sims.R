#################### IMPORTING DATASET ####################


#Includes one record per kidney, pancreas, or kidney-pancreas waiting list registration and/or transplant.


#If patient received a living donor transplant without being placed on the waiting list, there will be transplant information, but no waiting list information.

#If patient was listed for transplant, but not transplanted, there will be waiting list information but no transplant data.

#Registration events identified by WL_ID_CODE. Transplant events identified by TRR_ID_CODE.

library(mosaic)
library(readr)

KIDPAN_DATA <- read_delim("UNOS data/Delimited Text File 202206/Kidney_ Pancreas_ Kidney-Pancreas/KIDPAN_DATA.DAT", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)


colnames(KIDPAN_DATA) <- c("WL_ORG", "COD_WL", "COD_OSTXT_WL", "NUM_PREV_TX", "CURRENT_PRA", "PEAK_PRA", "USE_WHICH_PRA", "CREAT_CLEAR", "GFR", "DONATION", "ON_DIALYSIS", "MAX_KDPI_LOCAL_ZERO_ABDR", "MAX_KDPI_LOCAL_NON_ZERO_ABDR", "MAX_KDPI_IMPORT_ZERO_ABDR", "MAX_KDPI_IMPORT_NON_ZERO_ABDR", "C_PEPTIDE", "C_PEPTIDEDATE", "A2A2B_ELIGIBILITY", "A1", "A2", "B1", "B2", "DR1", "DR2", "ANTIBODY_TESTED", "GENDER", "ABO", "WGT_KG_TCR", "HGT_CM_TCR", "BMI_TCR", "CITIZENSHIP", "CITIZEN_COUNTRY", "PERM_STATE", "EDUCATION", "FUNC_STAT_TCR", "DGN_TCR", "DGN_OSTXT_TCR", "DGN2_TCR", "DGN2_OSTXT_TCR", "DIAB", "DRUGTRT_COPD", "TOT_SERUM_ALBUM", "C_PEPTIDE_PA_TCR", "HBA1C_PA_TCR", "SSDMF_DEATH_DATE", "INIT_CURRENT_PRA", "INIT_PEAK_PRA", "INIT_STAT", "INIT_WGT_KG", "INIT_HGT_CM", "INIT_CPRA", "END_CPRA", "INIT_EPTS", "END_EPTS", "REM_CD", "DAYSWAIT_CHRON", "END_STAT", "INIT_AGE", "ACTIVATE_DATE", "CREAT_CLEAR_DATE", "DEATH_DATE", "DIALYSIS_DATE", "END_DATE", "GFR_DATE", "INIT_DATE", "WT_QUAL_DATE", "ETHNICITY", "ETHCAT", "PT_CODE", "INIT_BMI_CALC", "END_BMI_CALC", "DAYSWAIT_ALLOC", "COMPOSITE_DEATH_DATE", "WLHR", "WLHL", "WLIN", "WLKI", "WLKP", "WLLI", "WLLU", "WLPA", "WLPI", "WLVC", "REGION", "INACT_REASON_CD", "BW4", "BW6", "C1", "C2", "DR51", "DR51_2", "DR52", "DR52_2", "DR53", "DR53_2", "DQ1", "DQ2", "WL_ID_CODE", "PERIP_VASC", "EXH_PERIT_ACCESS", "AGE_DIAB", "EXH_VASC_ACCESS", "YR_ENTRY_US_TCR", "WORK_INCOME_TCR", "ACADEMIC_PRG_TCR", "ACADEMIC_LEVEL_TCR", "MALIG_TCR_KI", "PRI_PAYMENT_TCR_KI", "MALIG_TCR_PA", "PRI_PAYMENT_TCR_PA", "PREV_TX", "PREV_KI_TX", "PREV_PA_TX", "ACADEMIC_LEVEL_TRR", "ACADEMIC_PRG_TRR", "FUNC_STAT_TRR", "MALIG_TRR", "MALIG_OSTXT_TRR", "MALIG_TY_TRR", "PERM_STATE_TRR", "PRI_PAYMENT_TRR_KI", "PRI_PAYMENT_CTRY_TRR_KI", "WORK_INCOME_TRR", "TX_DATE", "ACUTE_REJ_EPI_KI", "CREAT_TRR", "FIRST_WK_DIAL", "ORG_REC_ON", "PREV_PREG", "REC_ON_ICE", "REC_ON_PUMP", "SERUM_CREAT", "L_FIN_FLOW_RATE_TX", "L_FIN_RESIST_TX", "R_FIN_FLOW_RATE_TX", "R_FIN_RESIST_TX", "PRE_TX_TXFUS", "FIN_RESIST_TX", "TXHRT", "TXINT", "TXKID", "TXLIV", "TXLNG", "TXPAN", "TXVCA", "PRI_PAYMENT_CTRY_TCR_KI", "PREV_MALIG_TY", "PREV_MALIG_TY_OSTXT", "RDA1", "RDA2", "RDB1", "RDB2", "RDDR1", "RDDR2", "DON_RETYP", "RESUM_MAINT_DIAL_DT", "DA1", "DA2", "DB1", "DB2", "DDR1", "DDR2", "RA1", "RA2", "RB1", "RB2", "RDR1", "RDR2", "AMIS", "BMIS", "DRMIS", "HLAMIS", "NPKID", "NPPAN", "END_CPRA_DETAIL", "HAPLO_TY_MATCH_DON", "AGE_DON", "DDAVP_DON", "CMV_OLD_LIV_DON", "CMV_DON", "CMV_TEST_DON", "EBV_TEST_DON", "HBV_TEST_DON", "HCV_TEST_DON", "CMV_NUCLEIC_DON", "CMV_IGG_DON", "CMV_IGM_DON", "EBV_DNA_DON", "EBV_IGG_DON", "EBV_IGM_DON", "HBV_CORE_DON", "HBV_SUR_ANTIGEN_DON", "ETHCAT_DON", "COD_CAD_DON", "DEATH_CIRCUM_DON", "DEATH_MECH_DON", "CITIZENSHIP_DON", "HEP_C_ANTI_DON", "HCV_RNA_DON", "ABO_DON", "DON_TY", "GENDER_DON", "HOME_STATE_DON", "WARM_ISCH_TM_DON", "HCV_RIBA_DON", "HCV_ANTIBODY_DON", "LIV_DON_TY", "CITIZEN_COUNTRY_DON", "COD_OSTXT_DON", "CONTROLLED_DON", "CORE_COOL_DON", "NON_HRT_DON", "ANTIHYPE_DON", "BLOOD_INF_DON", "BLOOD_INF_CONF_DON", "BUN_DON", "CREAT_DON", "DOBUT_DON_OLD", "DOPAMINE_DON_OLD", "HTLV1_OLD_DON", "HTLV2_OLD_DON", "OTH_DON_MED1_OSTXT_DON_OLD", "OTH_DON_MED2_OSTXT_DON_OLD", "OTH_DON_MED3_OSTXT_DON_OLD", "OTHER_INF_DON", "OTHER_INF_CONF_DON", "OTHER_INF_OSTXT_DON", "PRETREAT_MED_DON_OLD", "PT_DIURETICS_DON", "PT_STEROIDS_DON", "PT_T3_DON", "PT_T4_DON", "PT_OTH2_OSTXT_DON", "PT_OTH3_OSTXT_DON", "PT_OTH4_OSTXT_DON", "PT_OTH1_OSTXT_DON", "PULM_INF_DON", "PULM_INF_CONF_DON", "SGOT_DON", "SGPT_DON", "TBILI_DON", "URINE_INF_DON", "URINE_INF_CONF_DON", "VASODIL_DON", "VDRL_DON", "CLIN_INFECT_DON", "HYPERTENS_DUR_DON", "CANCER_FREE_INT_DON", "CANCER_OTH_OSTXT_DON", "CONTIN_ALCOHOL_OLD_DON", "CONTIN_CIG_DON", "CONTIN_IV_DRUG_OLD_DON", "CONTIN_COCAINE_DON", "CONTIN_OTH_DRUG_DON", "DIET_DON", "DIURETICS_DON", "EXTRACRANIAL_CANCER_DON", "HIST_ALCOHOL_OLD_DON", "CANCER_SITE_DON", "HIST_CIG_DON", "DIABDUR_DON", "HIST_COCAINE_DON", "HIST_HYPERTENS_DON", "HIST_IV_DRUG_OLD_DON", "INSULIN_DEP_DON", "INTRACRANIAL_CANCER_DON", "OTHER_HYPERTENS_MED_DON", "HIST_CANCER_DON", "HIST_INSULIN_DEP_DON", "INSULIN_DUR_DON", "HIST_DIABETES_DON", "HIST_OTH_DRUG_DON", "SKIN_CANCER_DON", "DIABETES_DON", "LIV_DON_TY_OSTXT", "HEPARIN_DON", "ARGININE_DON", "INSULIN_DON", "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "BMI_DON_CALC", "KDPI", "KDRI_MED", "KDRI_RAO", "HBV_NAT_DON", "HCV_NAT_DON", "HIV_NAT_DON", "END_STAT_KI", "CREAT6M", "CREAT1Y", "DIAL_DATE", "RETXDATE_KI", "FAILDATE_KI", "PUMP_KI", "ABO_MAT", "AGE", "DISTANCE", "RESUM_MAINT_DIAL", "DIAL_TRR", "DIAG_KI", "DIAG_OSTXT_KI", "COLD_ISCH_KI", "GRF_STAT_KI", "GRF_FAIL_CAUSE_OSTXT_KI", "GRF_FAIL_CAUSE_TY_KI", "DWFG_KI", "PRVTXDIF_KI", "GTIME_KI", "GSTATUS_KI", "COD_KI", "COD_OSTXT_KI", "COD2_KI", "COD2_OSTXT_KI", "COD3_KI", "COD3_OSTXT_KI", "DAYSWAIT_CHRON_KI", "TX_PROCEDUR_TY_KI", "TRTREJ1Y_KI", "TRTREJ6M_KI", "MULTIORG", "PRI_PAYMENT_TRR_PA", "PRI_PAYMENT_CTRY_TRR_PA", "ART_RECON", "ART_RECON_OSTXT", "DUCT_MGMT", "DUCT_MGMT_OSTXT", "GRF_PLACEM", "PRE_AVG_INSULIN_USED_TRR", "PRE_AVG_INSULIN_USED_OLD_TRR", "ACUTE_REJ_EPI_PA", "PA_PRESERV_TM", "VASC_MGMT", "VEN_EXT_GRF", "INSULIN_PA", "INSULIN_RESUMED_DATE_PA", "INSULIN_DOSAGE_PA", "INSULIN_DURATION_PA", "METHOD_BLOOD_SUGAR_CONTROL_PA", "BLOOD_SUGAR_MEDICATION_PA", "BLOOD_SUGAR_MED_RESUMED_DATE_PA", "BLOOD_SUGAR_DIET_PA", "C_PEPTIDE_PA_TRR", "HBA1C_PA_TRR", "INSULIN_DOSAGE_OLD_PA", "PRI_PAYMENT_CTRY_TCR_PA", "PK_DA1", "PK_DA2", "PK_DB1", "PK_DB2", "PK_DDR1", "PK_DDR2", "ENTERIC_DRAIN", "ENTERIC_DRAIN_DT", "END_STAT_PA", "FAILDATE_PA", "DIAG_PA", "DIAG_OSTXT_PA", "GRF_STAT_PA", "GRF_FAIL_CAUSE_OSTXT_PA", "GRF_FAIL_CAUSE_TY_PA", "OTH_GRF_FAIL_CAUSE_OSTXT_PA", "GRF_VASC_THROMB_PA", "INFECT_PA", "BLEED_PA", "ANAST_LK_PA", "REJ_ACUTE_PA", "REJ_HYPER_PA", "BIOP_ISLET_PA", "PANCREATIT_PA", "REJ_CHRONIC_PA", "PX_NON_COMPL_PA", "RETXDATE_PA", "PRVTXDIF_PA", "GTIME_PA", "GSTATUS_PA", "COD_PA", "COD_OSTXT_PA", "COD2_PA", "COD2_OSTXT_PA", "COD3_PA", "COD3_OSTXT_PA", "DAYSWAIT_CHRON_PA", "TX_PROCEDUR_TY_PA", "TRTREJ1Y_PA", "TRTREJ6M_PA", "ORGAN", "CMV_IGG", "CMV_IGM", "EBV_SEROSTATUS", "HBV_CORE", "HBV_SUR_ANTIGEN", "HCV_SEROSTATUS", "HIV_SEROSTATUS", "CMV_STATUS", "HBV_SURF_TOTAL", "HIV_NAT", "HCV_NAT", "HBV_NAT", "PREV_TX_ANY", "PREV_TX_ANY_N", "TX_TYPE", "MED_COND_TRR", "PX_STAT", "PX_STAT_DATE", "PREV_KI_DATE", "FUNC_STAT_TRF", "SHARE_TY", "PSTATUS", "PTIME", "LOS", "PAYBACK", "ECD_DONOR", "AGE_GROUP", "MALIG", "MALIG_TY_OSTXT", "MALIG_TY", "HGT_CM_CALC", "WGT_KG_CALC", "BMI_CALC", "STATUS_TCR", "STATUS_TRR", "STATUS_DDR", "VAL_DT_DDR", "STATUS_LDR", "VAL_DT_LDR", "VAL_DT_TCR", "VAL_DT_TRR", "LT_ONE_WEEK_DON", "REJ_BIOPSY", "REJCNF_KI", "REJTRT_KI", "REJCNF_PA", "REJTRT_PA", "TRR_ID_CODE", "ADMISSION_DATE", "DISCHARGE_DATE", "COMPL_ABSC", "COMPL_ANASLK", "COMPL_PANCREA", "OTH_COMPL_OSTXT", "SURG_INCIS", "OPER_TECH", "EDUCATION_DON", "KI_CREAT_PREOP", "KI_PROC_TY", "PRI_PAYMENT_DON", "PRI_PAYMENT_CTRY_DON", "MEDICARE_DON", "MEDICAID_DON", "OTH_GOVT_DON", "PRIV_INS_DON", "HMO_PPO_DON", "SELF_DON", "DONATION_DON", "FREE_DON", "RECOV_OUT_US", "RECOV_COUNTRY", "PROTEIN_URINE", "LIPASE", "AMYLASE", "INOTROP_AGENTS", "CARDARREST_NEURO", "RESUSCIT_DUR", "INOTROP_SUPPORT_DON", "TATTOOS", "LT_KI_BIOPSY", "LT_KI_GLOMERUL", "RT_KI_BIOPSY", "RT_KI_GLOMERUL", "REFERRAL_DATE", "RECOVERY_DATE", "ADMIT_DATE_DON", "DONOR_ID", "HBSAB_DON", "EBV_IGG_CAD_DON", "EBV_IGM_CAD_DON", "HBV_DNA_DON", "CDC_RISK_HIV_DON", "INO_PROCURE_AGENT_1", "INO_PROCURE_AGENT_2", "INO_PROCURE_AGENT_3", "INO_PROCURE_OSTXT_1", "INO_PROCURE_OSTXT_2", "INO_PROCURE_OSTXT_3", "DATA_TRANSPLANT", "DATA_WAITLIST", "CTR_CODE", "OPO_CTR_CODE", "INIT_OPO_CTR_CODE", "END_OPO_CTR_CODE", "LISTING_CTR_CODE")


# Deceased donors only
data <- KIDPAN_DATA[KIDPAN_DATA$DON_TY == "C",]

### Compute number of days on waitlist
data$WL_days <- difftime(as.Date(data$END_DATE, format = "%m/%d/%Y"), as.Date(data$INIT_DATE, format = "%m/%d/%Y"),units = "days")


attach(data)




#######################################################################################################################################
# AIM: Use existing microsimulation models to evaluate potential inequality of organ 
#       allocation policies in existing deceased kidney organ donation programs.
#######################################################################################################################################


#######################################################################################################################################
################################################ POTENTIAL VARIABLES OF NOTE: #########################################################
#######################################################################################################################################

# ABO_DON -> DONOR BLOOD TYPE

# CITIZENSHIP -> CANDIDATE CITIZENSHIP
# CITIZENSHIP_DON -> DONOR CITIZENSHIP
# CITIZEN_COUNTRY	-> CANDIDATE COUNTRY OF PERMANENT RESIDENCE @ REGISTRATION
# CITIZEN_COUNTRY_DON -> 	DECEASED DONOR-HOME COUNTRY


# DAYSWAIT_ALLOC -> TIME USED FOR ALLOCATION PRIORITY (DAYS)
# DAYSWAIT_CHRON_KI -> TOTAL DAYS ON KIDNEY WAITING LIST

# ETHNICITY -> RECIPIENT ETHNICITY

# GENDER -> RECIPIENT GENDER
# GENDER_DON -> DONOR GENDER

# GSTATUS_KI	-> GRAFT FAILED (1=YES)-KIDNEY
# GTIME_KI	-> GRAFT LIFESPAN-KIDNEY-Days From Transplant to Failure/Death/Last Follow-Up

# HGT_CM_CALC	-> CALCULATED RECIPIENT HEIGHT(cm)
# HGT_CM_DON_CALC -> 	CALCULATED DONOR HEIGHT (CM)

# HIST_ALCOHOL_OLD_DON ->	DECEASED DONOR-HISTORY OF ALCOHOL DEPENDENCY
# HIST_CANCER_DON	-> DECEASED DONOR-HISTORY OF CANCER (Y/N)
# HIST_CIG_DON -> 	DECEASED DONOR-HISTORY OF CIGARETTES IN PAST @ >20PACK YRS
# HIST_COCAINE_DON ->	DECEASED DONOR-HISTORY OF COCAINE USE IN PAST
# HIST_DIABETES_DON ->	DECEASED DONOR-HISTORY OF DIABETES, INCL. DURATION OF DISEASE
# HIST_HYPERTENS_DON ->	DECEASED DONOR-HISTORY OF HYPERTENSION (Y,N)
# HIST_INSULIN_DEP_DON ->	DECEASED DONOR-INSULIN DEPENDENT DIABETES (Y,N)
# HIST_IV_DRUG_OLD_DON ->	DECEASED DONOR-HISTORY OF IV DRUG USE IN PAST
# HIST_OTH_DRUG_DON	-> DECEASED DONOR-HISTORY OF OTHER DRUG USE IN PAST

# INIT_AGE ->	CANDIDATE AGE IN YEARS AT TIME OF LISTING


# LIV_DON_TY ->	LIVING DONOR RELATION TO RECIPIENT

# NON_HRT_DON	DECEASED ->  DONOR-NON-HEART BEATING DONOR

# INIT_DATE -> DATE PLACED ON WAITING LIST
# END_DATE ->	EARLIEST OF DATES OF REMOVAL FROM WAITING LIST, TRANSPLANT, DEATH, OR TIME COPY OF DATA CREATED
# END_EPTS ->	Offer/Removal/Current Calculated EPTS (since 5/27/2014)


# DON_TY ->	DONOR TYPE - DECEASED OR LIVING
# DONATION	-> WL Will Receive Donation Points due to a Previous Living Organ Donation
# DONATION_DON ->	LIVING DONOR-DONATION AS SECONDARY PAYMENT TYPE

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################




#####################################################################################################################################
############################## SIMULATING THE TOP TRADING CYCLES AND CHAINS (TTCC) MECHANISM ########################################
#####################################################################################################################################
############################## https://rdrr.io/cran/matchingMarkets/man/ttcc.html ################################################### 
#####################################################################################################################################
source("ttcc.r")
#####################################################################################################################################
################################################ Sample Simulation ##################################################################
prefs <- matrix(c( 9,10, 1,NA,NA,NA,NA,
                  11, 3, 5, 6, 2,NA,NA,
                  2, 4, 5, 6, 7, 8,13,
                   5, 9, 1, 8,10, 3,13,
                   3, 7,11, 4, 5,NA,NA,
                   3, 5, 8, 6,NA,NA,NA,
                   6, 1, 3, 9,10, 1,13,
                   6, 4,11, 2, 3, 8,NA,
                   3,11,13,NA,NA,NA,NA,
                  11, 1, 4, 5, 6, 7,13,
                   3, 6, 5,11,NA,NA,NA,
                  11, 3, 9, 8,10,12,NA),
              byrow = FALSE, ncol = 12) 

priority <- 1:12


#####################################################################################################################################
#####################################################################################################################################

# Run the ttcc() function and produce patient-donor pairings
pairing <- ttcc(prefs = prefs, priority = priority)$matching


# Create graphical model of patient-donor pairings
library(igraph)

# Create a data frame for individuals and objects
individuals <- data.frame(id = unique(pairing$ind), label = paste("Patient", unique(pairing$ind)))
objects <- data.frame(id = unique(pairing$obj), label = paste("Kidney", unique(pairing$obj)))

nodes <- rbind(individuals, objects)

# Create unique ID for each vertex
node_ids <- 1:nrow(nodes)

graph <- graph.empty(n = nrow(nodes), directed = T)


V(graph)$name <- nodes$label
V(graph)$color <- ifelse(V(graph)$name %in% individuals$label, "coral1", "lightblue")
V(graph)$size <- ifelse(V(graph)$name %in% individuals$label, 10, 15)

layout <- layout_with_fr(graph)

# Add edges connecting patients to kidneys
for (i in 1:nrow(pairing)) {
  ind_id <- pairing$ind[i]
  obj_id <- pairing$obj[i]
  ind_vertex <- paste("Patient", ind_id)
  obj_vertex <- paste("Kidney", obj_id)
  
  # Check if the vertex is a patient
  if (ind_vertex %in% individuals$label) {
    # Check if the target vertex is not a patient
    if (!(obj_vertex %in% individuals$label)) {
      graph <- add_edges(graph, c(ind_vertex, obj_vertex))
    }
  }
}

# Plot graph
plot(graph, layout = layout, vertex.label = V(graph)$name, vertex.color = V(graph)$color,
     vertex.size = V(graph)$size, edge.arrow.size = 1,  vertex.label.dist = 1.9)

#######################################################################################################################################
#######################################################################################################################################



################################## Data summaries ##################################
hist(as.numeric(GTIME_KI), main="Histogram of Graft Lifespan for Kidney Transpants", xlab="Days")



######################## DONOR-RECIPIENT ABO MATCH LEVEL ######################## 
# Identical	1
# Compatible	2
# Incompatible 3
tally(ABO_MAT)


######################## Kidney Donor Profile Index (KDPI) ########################
KDPI_decimal <-(sapply(KDPI, function(x) {
  if (is.na(x) || x == ".") {
    NA
  } else {
    as.numeric(sub("%", "", x)) / 100
  }
}))


hist(as.numeric(KDPI_decimal), main = "Histogram of Kidney Donor Profile Index (KDPI)", xlab="Kidney Donor Profile Index (KDPI)")


######################## Estimated Post-Transplant Survival (EPTS) Score ########################
EPTS_decimal <-(sapply(END_EPTS, function(x) {
  if (is.na(x) || x == ".") {
    NA
  } else {
    as.numeric(sub("%", "", x)) / 100
  }
}))

hist(as.numeric(EPTS_decimal), main="Histogram of Estimated Post-Transplant Survival (EPTS) Score", xlab="Estimated Post-Transplant Survival")
#epts_model <- lm(as.numeric(WL_days)~EPTS_decimal)
