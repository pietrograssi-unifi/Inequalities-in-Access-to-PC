# ==============================================================================
# SCRIPT: Social Determinants of Palliative Care Access in Europe
# DATA SOURCE: Survey of Health, Ageing and Retirement in Europe (SHARE)
# AUTHORS: Pietro Grassi, MD Arianna Bellini, PhD Chiara Seghieri, PhD Daniele Vignoli
# LAST UPDATED: January 2, 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. ENVIRONMENT SETUP & LIBRARIES
# ------------------------------------------------------------------------------
rm(list=ls())

# Check and install required packages (using 'pacman' for efficiency)
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  tidyverse,    # Data manipulation & plotting (ggplot2, dplyr)
  haven,        # Read Stata (.dta) files
  fs,           # File system operations
  gtsummary,    # Publication-ready descriptive tables
  gt,           # Advanced table formatting
  lme4,         # Multilevel Mixed-Effects Models (GLMM)
  sjPlot,       # Visualization of regression models (Odds Ratios)
  performance,  # Model performance indices (ICC, R2)
  mice,         # Multiple imputation
  broom.mixed,  #
  mitml         #
)

# SET WORKING DIRECTORY
# NOTE: Ensure this path contains: easySHARE_rel9-0-0.dta AND all wave-specific subfiles (xt, ph)
setwd("/Users/pietro/Desktop/Sant'Anna/Ricerca/SHARE/Dati") 

# ------------------------------------------------------------------------------
# 1. OUTCOME EXTRACTION: END-OF-LIFE MODULE (XT)
# ------------------------------------------------------------------------------
# Objective: Identify deceased respondents and their palliative care usage.
# The 'xt' module is administered to a proxy respondent after the participant's death.

# Locate all End-of-Life (XT) files across waves
files_xt <- dir_ls(regexp = "sharew.*_xt\\.dta$")

# Define target columns to extract
target_cols <- c("mergeid", "country", "xt022_", "xt012_", "xt005_y")

message(">>> Importing End-of-Life (XT) modules...")

df_eol <- map_dfr(files_xt, function(f) {
  d <- read_dta(f)
  # Robust selection: select only columns present in the specific wave
  d_selected <- d %>%
    select(any_of(target_cols)) %>%
    mutate(wave_death = as.numeric(str_extract(f, "(?<=sharew)\\d+")))
  return(d_selected)
}) %>%
  mutate(
    # --- DEPENDENT VARIABLE: PALLIATIVE CARE ACCESS ---
    # Question: "Did [Name] receive hospice or palliative care services?"
    # Coding in SHARE: 1 = Yes, 5 = No.
    # Recoding: 1 = Yes, 0 = No.
    palliative_access = if (exists("xt022_", where = .)) {
      case_when(xt022_ == 1 ~ 1, xt022_ == 5 ~ 0, TRUE ~ NA_real_)
    } else { NA_real_ },
    
    # Place of Death (for descriptive purposes)
    place_death = if (exists("xt012_", where = .)) {
      to_character(xt012_) 
    } else { NA_character_ }
  ) %>%
  # Exclusion criteria: Remove cases with missing outcome data
  filter(!is.na(palliative_access))

message(paste(">>> Total deceased respondents with valid outcome:", nrow(df_eol)))

# ------------------------------------------------------------------------------
# 2. CLINICAL PROFILE CONSTRUCTION (PH MODULE)
# ------------------------------------------------------------------------------
# Objective: Reconstruct the clinical trajectory (Cancer vs. Organ Failure vs. Frailty).
# We retrieve chronic conditions from the last interview while alive (Physical Health module).

files_ph <- dir_ls(regexp = "sharew.*_ph\\.dta$")

message(">>> Importing Physical Health (PH) modules for clinical profiling...")

df_health_raw <- map_dfr(files_ph, function(f) {
  read_dta(f) %>%
    select(
      mergeid, 
      matches("ph006"), # Chronic diseases (dummy variables)
      matches("ph084")  # Pain occurrence
    ) %>%
    mutate(wave_interview = as.numeric(str_extract(f, "(?<=sharew)\\d+")))
})

# Define expected columns to handle version changes between waves (Legacy vs New codes)
# Codes: 1=Heart, 2=Stroke, 4=Cancer, 11=Lung, 14/16=Dementia
expected_cols <- c(
  "ph006_1", "ph006d1",   # Heart Attack / Problems
  "ph006_2", "ph006d2",   # Stroke
  "ph006_4", "ph006d4",   # Cancer
  "ph006_11", "ph006d11", # Lung Disease
  "ph006_16", "ph006d16", "ph006d14", # Alzheimer/Dementia
  "ph084_", "ph084"       # Pain
)

# Impute missing columns with NA to ensure stability
for(col in expected_cols) {
  if(!col %in% names(df_health_raw)) df_health_raw[[col]] <- NA
}

df_health_clean <- df_health_raw %>%
  # Remove Stata labels for numeric processing
  mutate(across(where(is.labelled), ~as.numeric(zap_labels(.)))) %>% 
  mutate(
    # 1. CANCER TRAJECTORY
    has_cancer = case_when(
      (!is.na(ph006_4) & ph006_4 == 1) | (!is.na(ph006d4) & ph006d4 == 1) ~ 1, 
      TRUE ~ 0),
    
    # 2. DEMENTIA / NEURODEGENERATIVE TRAJECTORY
    has_alzheimer = case_when(
      (!is.na(ph006_16) & ph006_16 == 1) | 
        (!is.na(ph006d16) & ph006d16 == 1) | 
        (!is.na(ph006d14) & ph006d14 == 1) ~ 1, 
      TRUE ~ 0),
    
    # 3. ORGAN FAILURE TRAJECTORY (Heart + Stroke + Lung)
    has_heart_stroke_lung = case_when(
      (!is.na(ph006_1) & ph006_1 == 1) | (!is.na(ph006d1) & ph006d1 == 1) ~ 1, 
      (!is.na(ph006_2) & ph006_2 == 1) | (!is.na(ph006d2) & ph006d2 == 1) ~ 1, 
      (!is.na(ph006_11) & ph006_11 == 1) | (!is.na(ph006d11) & ph006d11 == 1) ~ 1, 
      TRUE ~ 0
    ),
    
    # 4. PAIN SYMPTOM
    has_pain = case_when(
      (!is.na(ph084_) & ph084_ == 1) | (!is.na(ph084) & ph084 == 1) ~ 1, 
      TRUE ~ 0)
  ) %>%
  # Keep only the last available health interview before death
  group_by(mergeid) %>%
  arrange(wave_interview) %>%
  slice_tail(n = 1) %>% 
  ungroup() %>%
  select(mergeid, has_cancer, has_alzheimer, has_heart_stroke_lung, has_pain)

# ------------------------------------------------------------------------------
# 3a. SOCIO-ECONOMIC & CONTEXTUAL DATA (easySHARE)
# ------------------------------------------------------------------------------
# Objective: Retrieve standardized SES, Geography (Urban/Rural), and Social Support variables.

message(">>> Importing easySHARE for Socio-Economic variables...")

df_easy <- read_dta("easySHARE_rel9-0-0.dta") %>%
  select(
    mergeid, wave, female, age, isced1997_r, thinc_m, 
    eurod, hc012_, hc002_mod, adla,
    iv009_mod,     # Area
    sp002_mod,     # Help received
    hhsize,        # People in the household
    partnerinhh    # Partner in the household
  ) %>%
  mutate(
    # Harmonize Hospitalization (1=Yes, 0=No)
    hosp = case_when(hc012_ == 1 ~ 1, hc012_ == 5 ~ 0, TRUE ~ NA_real_),
    doc_total = hc002_mod,
    thinc = thinc_m,
    
    # Clean Area Variable (Remove negative missing codes)
    area_clean = case_when(
      iv009_mod >= 1 & iv009_mod <= 5 ~ iv009_mod,
      TRUE ~ NA_real_
    ),
    
    # Social Support Dummy (1=Yes, 0=No)
    # sp002_mod: "Received help from others"
    social_support = case_when(
      sp002_mod == 1 ~ 1, 
      sp002_mod == 0 | sp002_mod == 5 ~ 0, 
      TRUE ~ 0
    )
  ) %>%
  group_by(mergeid) %>%
  arrange(wave) %>%
  slice_tail(n = 1) %>% 
  ungroup()

# ------------------------------------------------------------------------------
# 3b. MACRO-INDICATOR: PC COVERAGE (Sánchez-Cárdenas et al., 2021)
# ------------------------------------------------------------------------------

df_pc_coverage <- tribble(
  ~country, ~pc_coverage_score,
  "Netherlands",    51,
  "Germany",        43,
  "Switzerland",    40,
  "Belgium",        38,
  "Czech Republic", 36,
  "Denmark",        36,
  "Italy",          36,
  "Poland",         35,
  "Spain",          35,
  "Austria",        34,
  "France",         33,
  "Hungary",        25,
  "Portugal",       25,
  "Israel",         20,
  "Lithuania",      19,
  "Sweden",         19,
  "Greece",         17,
  "Luxembourg",     17,
  "Slovenia",       15,
  "Finland",        13,
  "Malta",          11,
  "Cyprus",         10,
  "Romania",        21,
  "Bulgaria",       8,
  "Croatia",        8,
  "Estonia",        8,
  "Latvia",         10,
  "Slovakia",       10
)

# ------------------------------------------------------------------------------
# 4. DATA MERGING & THEORETICAL OPERATIONALIZATION
# ------------------------------------------------------------------------------

message(">>> Merging datasets and creating analysis variables (Macro-Meso-Micro Framework)...")

df_final_merged <- df_eol %>%
  left_join(df_easy, by = "mergeid") %>%
  left_join(df_health_clean, by = "mergeid") %>%
  mutate(country_label_temp = as.character(as_factor(country))) %>%
  left_join(df_pc_coverage, by = c("country_label_temp" = "country")) %>%
  filter(!is.na(isced1997_r)) %>%
  filter(!is.na(country))

df_analysis <- df_final_merged %>%
  mutate(
    # --- LEVEL 1: MICRO (Individual Predisposing & Need Factors) ---
    age = as.numeric(age),
    female = as_factor(female),
    country_lab = as_factor(country),
    wave_death_factor = as.factor(wave_death), # Time Trend Control
    household_size = as.numeric(hhsize),
    has_partner = factor(if_else(partnerinhh == 1, "Yes", "No")),
    pc_coverage_score = as.numeric(pc_coverage_score),
    
    # Social Gradient (SES - Predisposing)
    education = factor(case_when(
      isced1997_r <= 2 ~ "Low",
      isced1997_r <= 4 ~ "Medium",
      isced1997_r >= 5 ~ "High"
    ), levels = c("Low", "Medium", "High")),
    
    # Clinical Trajectories (Need Factors - Hierarchy)
    disease_trajectory = factor(case_when(
      has_cancer == 1 ~ "1. Cancer",
      has_alzheimer == 1 ~ "2. Dementia/Neuro",
      has_heart_stroke_lung == 1 ~ "3. Organ Failure",
      TRUE ~ "4. Frailty/Other"
    ), levels = c("1. Cancer", "2. Dementia/Neuro", "3. Organ Failure", "4. Frailty/Other")),
    
    # Functional Status (Need Factors)
    disability_status = factor(if_else(adla > 0, "Dependent", "Independent")),
    hospital_user = factor(if_else(hosp == 1, "Yes", "No"), levels = c("Yes", "No")),
    
    # --- LEVEL 2: MESO (Structural & Enabling Factors) ---
    # Geographic Gradient (Distance Decay Proxy)
    living_area = factor(case_when(
      iv009_mod %in% c(1, 2) ~ "Urban/City",
      iv009_mod %in% c(3, 4) ~ "Town",
      iv009_mod == 5 ~ "Rural",
      TRUE ~ NA_character_
    ), levels = c("Urban/City", "Town", "Rural")),
    
    # Social Capital
    has_social_support = factor(if_else(social_support == 1, "Yes", "No")),
    
    # Income
    income_tertile = ntile(thinc, 3), 
    income_group = factor(case_when(
      income_tertile == 1 ~ "Low",
      income_tertile == 2 ~ "Medium",
      income_tertile == 3 ~ "High"
    ), levels = c("Low", "Medium", "High")),
    
    # --- LEVEL 3: MACRO (Systemic - Welfare Regimes) ---
    # Based on Esping-Andersen typology. 
    welfare_group = factor(case_when(
      country %in% c(11, 12, 13, 14, 15, 16) ~ "Continental",          # DE, AT, NL, FR, CH, BE
      country %in% c(17, 18) ~ "Nordic",                               # SE, DK
      country %in% c(23, 25, 28, 29, 31, 32, 33, 34, 35) ~ "Southern", # ES, IT, GR, PT, IL
      country %in% c(20, 30:60) ~ "Eastern",                           # PL, CZ, SI, EE, HR etc.
      TRUE ~ "Continental" # Fallback to avoid NA
    ), levels = c("Continental", "Nordic", "Southern", "Eastern"))
  )

# EXPORT PROCESSED DATA FOR REPRODUCIBILITY
write_csv(df_analysis, "share_palliative_analysis_data.csv")
message(">>> Processed dataset exported to 'share_palliative_analysis_data.csv'.")

message(">>> Generating Descriptive Table (Table 1)...")

table1 <- df_analysis %>%
  select(
    palliative_access,  
    pc_coverage_score,
    welfare_group,          
    living_area, has_social_support, 
    household_size, has_partner,
    age, female, education, disease_trajectory, disability_status, hospital_user 
  ) %>%
  mutate(across(where(is.factor), droplevels)) %>% 
  tbl_summary(
    by = welfare_group,     
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      palliative_access ~ "Access to Palliative Care",
      living_area ~ "Area of Residence",
      has_social_support ~ "Social Support Received",
      household_size ~ "Household Size",
      has_partner ~ "Living with Partner",
      age ~ "Age at Death",
      pc_coverage_score ~ "Integration Capacity Score (ICS)",
      female ~ "Gender",
      education ~ "Education Level",
      disease_trajectory ~ "Disease Trajectory",
      disability_status ~ "Functional Disability",
      hospital_user ~ "Hospitalized in Last Year"
    ),
    missing = "no" 
  ) %>%
  add_overall() %>% 
  add_p(test = all_categorical() ~ "chisq.test") %>%       
  bold_labels() %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Table 1. Descriptive Characteristics of the Deceased Population by Welfare Regime**")

print(table1)

table1 %>%
  as_gt() %>%
  gtsave(filename = "Table1.png", path = getwd())

# ------------------------------------------------------------------------------
# 5. MULTIPLE IMPUTATION STRATEGY (MICE)
# ------------------------------------------------------------------------------
message(">>> Preparing data for Multiple Imputation...")

df_for_imputation <- df_analysis %>%
  select(
    mergeid, palliative_access, 
    education, living_area, has_social_support, disability_status, hospital_user,
    household_size, has_partner, age, female, wave_death_factor, disease_trajectory,
    welfare_group, country_lab, pc_coverage_score
  )

message(paste(">>> Total sample size for analysis:", nrow(df_for_imputation)))

# --- MODEL 0: NULL MODEL (BASELINE) ---
m0_null <- glmer(
  palliative_access ~ 1 + (1 | country_lab), 
  data = df_for_imputation, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

m0_fixed <- glm(palliative_access ~ 1, 
                data = df_for_imputation, 
                family = binomial)

lrt_test <- anova(m0_null, m0_fixed)

message(">>> Likelihood Ratio Test for Country Effect:")
print(lrt_test)

icc0 <- performance::icc(m0_null)$ICC_adjusted
message(paste(">>> Baseline ICC (Model 0):", round(icc0, 4)))

# --- RUNNING MICE ---
message(">>> Running Multiple Imputation (m=5)... this may take a moment.")
pred_matrix <- quickpred(df_for_imputation, mincor = 0.1)
imp_data <- mice(df_for_imputation, predictorMatrix = pred_matrix, m = 5, maxit = 5, seed = 123, print = TRUE)

# ------------------------------------------------------------------------------
# 6. POOLED MULTILEVEL MODELS
# ------------------------------------------------------------------------------
message(">>> Running Pooled Models on Imputed Data...")

glmer_opts_robust <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# --- MODEL 1: MICRO ---
fit_m1_imp <- with(imp_data, 
                   glmer(palliative_access ~ 
                           education + disease_trajectory + disability_status + hospital_user + 
                           age + female + wave_death_factor + household_size + has_partner + 
                           (1 | country_lab), 
                         family = binomial,
                         control = glmer_opts_robust)
)

# --- MODEL 2: FULL ---
fit_m2_imp <- with(imp_data, 
                   glmer(palliative_access ~ 
                           welfare_group + living_area + has_social_support + household_size + has_partner +
                           education + disease_trajectory + disability_status + hospital_user + 
                           age + female + wave_death_factor + 
                           (1 | country_lab), 
                         family = binomial,
                         control = glmer_opts_robust)
)

# --- MODEL 3: COVERAGE ---
fit_m3_coverage <- with(imp_data, 
                        glmer(palliative_access ~ 
                                pc_coverage_score +
                                living_area + has_social_support + 
                                household_size + has_partner + 
                                education + disease_trajectory + disability_status + hospital_user + 
                                age + female + wave_death_factor + 
                                (1 | country_lab), 
                              family = binomial,
                              control = glmer_opts_robust)
)

# ------------------------------------------------------------------------------
# 7. EXPORT RESULTS
# ------------------------------------------------------------------------------

get_pooled_icc <- function(pooled_model) {
  vars <- sapply(pooled_model$analyses, function(x) {
    as.data.frame(VarCorr(x))[1, "vcov"]
  })
  mean_var <- mean(vars)
  icc <- mean_var / (mean_var + (pi^2/3))
  return(sprintf("%.2f", icc)) 
}

icc_m1 <- get_pooled_icc(fit_m1_imp)
icc_m2 <- get_pooled_icc(fit_m2_imp)
icc_m3 <- get_pooled_icc(fit_m3_coverage)

prepare_model_table <- function(model) {
  tbl_regression(model, exponentiate = TRUE) %>%
    bold_p() %>%
    modify_table_body(
      ~ .x %>%
        mutate(
          est_txt = sprintf("%.2f", estimate),
          low_txt = sprintf("%.2f", conf.low),
          high_txt = sprintf("%.2f", conf.high),
          final_str = paste0(est_txt, " (", low_txt, ", ", high_txt, ")"),
          final_str = ifelse(is.na(estimate), "", final_str)
        ) %>%
        mutate(estimate = final_str)
    ) %>%
    modify_column_hide(columns = c(conf.low, conf.high, p.value))
}

t1_ready <- prepare_model_table(fit_m1_imp)
t2_ready <- prepare_model_table(fit_m2_imp)
t3_ready <- prepare_model_table(fit_m3_coverage)

table_final_comparison <- tbl_merge(
  tbls = list(t1_ready, t2_ready, t3_ready),
  tab_spanner = c("**Model 1**", "**Model 2**", "**Model 3**")
) %>%
  modify_fmt_fun(
    starts_with("estimate") ~ function(x) x 
  ) %>%
  modify_table_body(~ .x %>% 
                      add_row(
                        label = "ICC (Country)", 
                        estimate_1 = icc_m1,
                        estimate_2 = icc_m2,
                        estimate_3 = icc_m3,
                        row_type = "label"
                      )
  ) %>%
  modify_header(
    estimate_1 ~ "**OR (95% CI)**",
    estimate_2 ~ "**OR (95% CI)**",
    estimate_3 ~ "**OR (95% CI)**"
  ) %>%
  bold_labels()

print(table_final_comparison)

table_final_comparison %>%
  as_gt() %>%
  gtsave("Table2.png")

# ------------------------------------------------------------------------------
# 8. VISUALIZATION
# ------------------------------------------------------------------------------
message(">>> Generating Plots...")

df_plot <- complete(imp_data, 1)

# ------------------------------------------------------------------------------
# A. FIGURE 1 & 2: MODEL 2 (WELFARE REGIMES)
# ------------------------------------------------------------------------------

m2_plot <- glmer(
  palliative_access ~ 
    welfare_group + living_area + has_social_support + 
    household_size + has_partner +   
    education + disease_trajectory + disability_status + hospital_user + 
    age + female + wave_death_factor + 
    (1 | country_lab), 
  data = df_plot, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

my_labels <- c(
  "has_partnerYes"                      = "Living with Partner: Yes (ref: No)",
  "household_size"                      = "Household Size (continuous)",
  "has_social_supportYes"               = "Social Support: Yes (ref: No)",
  "female1. female"                     = "Female (ref: Male)",
  "living_areaRural"                    = "Area: Rural (ref: Urban)",
  "living_areaTown"                     = "Area: Town (ref: Urban)",
  "welfare_groupSouthern"               = "Welfare: Southern (ref: Continental)",
  "welfare_groupNordic"                 = "Welfare: Nordic (ref: Continental)",
  "welfare_groupEastern"                = "Welfare: Eastern (ref: Continental)",
  "pc_coverage_score"                   = "Integration Capacity Score",
  "disease_trajectory2. Dementia/Neuro" = "Diag: Dementia (ref: Cancer)",
  "disease_trajectory3. Organ Failure"  = "Diag: Organ Failure (ref: Cancer)",
  "disease_trajectory4. Frailty/Other"  = "Diag: Frailty (ref: Cancer)",
  "educationHigh"                       = "Edu: High (ref: Low)",
  "educationMedium"                     = "Edu: Medium (ref: Low)",
  "disability_statusIndependent"        = "Disability: Independent (ref: Dependent)",
  "hospital_userNo"                     = "Hospital User: No (ref: Yes)",
  "age"                                 = "Age (Years)" 
)

# --- FIGURE 1: FOREST PLOT (MODEL 2) ---
p1_base <- plot_model(
  m2_plot, 
  type = "est", 
  sort.est = TRUE,  
  transform = "exp", 
  show.values = FALSE,  
  show.p = FALSE,       
  rm.terms = c("wave_death_factor4", "wave_death_factor5", "wave_death_factor6", 
               "wave_death_factor7", "wave_death_factor8", "wave_death_factor9",
               "wave_death_factor[4-9]", "Intercept"), 
  vline.color = "black",
  colors = c("firebrick", "grey50", "#0072B2"),
  title = "Model 2: Determinants of Access (Welfare Regimes)",
  axis.title = "Odds Ratio (95% CI)"
)

p1_final <- p1_base +
  scale_x_discrete(labels = my_labels) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11, color = "black"), 
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted")
  ) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  labs(y = "Odds Ratio (Log Scale)")

print(p1_final)
ggsave("Figure1.pdf", plot = p1_final, width = 10, height = 8)

# --- FIGURE 2: RANDOM EFFECTS (MODEL 2) ---
p2 <- plot_model(m2_plot, type = "re", 
                 title = "Country-Level Variation (Model 2)",
                 sort.est = TRUE, 
                 grid = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_sjplot()

print(p2)
ggsave("Figure2.pdf", plot = p2, width = 8, height = 6)

# ------------------------------------------------------------------------------
# B. FIGURE 3: MODEL 3
# ------------------------------------------------------------------------------

m3_plot <- glmer(
  palliative_access ~ 
    pc_coverage_score + 
    living_area + has_social_support + 
    household_size + has_partner +   
    education + disease_trajectory + disability_status + hospital_user + 
    age + female + wave_death_factor + 
    (1 | country_lab), 
  data = df_plot, 
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

p3_base <- plot_model(
  m3_plot, 
  type = "est", 
  sort.est = TRUE,  
  transform = "exp", 
  show.values = FALSE,
  show.p = FALSE,
  rm.terms = c("wave_death_factor4", "wave_death_factor5", "wave_death_factor6", 
               "wave_death_factor7", "wave_death_factor8", "wave_death_factor9",
               "wave_death_factor[4-9]", "Intercept"), 
  vline.color = "grey30",
  colors = c("firebrick", "grey50", "#0072B2"),
  title = "Model 3: Sensitivity Analysis (ICS)",
  axis.title = "Odds Ratio (95% CI)"
)

p3_final <- p3_base +
  scale_x_discrete(labels = my_labels) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 11, color = "black"), 
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black")

print(p3_final)
ggsave("Figure3.pdf", plot = p3_final, width = 10, height = 8)

# ------------------------------------------------------------------------------
# APPENDIX. INTERACTION HOSPITAL * DIAGNOSIS
# ------------------------------------------------------------------------------
fit_m_interact <- with(imp_data, 
                       glmer(palliative_access ~ 
                               hospital_user * disease_trajectory + 
                               welfare_group + living_area + has_social_support + 
                               household_size + has_partner +   
                               education + disability_status + 
                               age + female + wave_death_factor + 
                               (1 | country_lab), 
                             family = binomial,
                             control = glmer_opts_robust)
)

pool_interact <- pool(fit_m_interact)
summary(pool_interact) %>% 
  select(term, estimate, p.value) %>% 
  filter(str_detect(term, ":")) 
