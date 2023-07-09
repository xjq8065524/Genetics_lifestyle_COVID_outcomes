# Load basic packages -----------------------------------------------------
library( tidyverse)
library( lubridate)
library( adjustedCurves)
library( survival)
library( survminer)
library(patchwork)
library( ggpmisc)
devtools::install_github("https://github.com/RobinDenz1/adjustedCurves")

# Universal functions -----------------------------------------------------
temp_func <- function( raw_model){
  
  temp <- summary( raw_model )
  
  p <- temp$coefficients %>% as.data.frame()%>% janitor::clean_names() %>% pull( pr_z)
  
  Cox_results <- 
    temp$conf.int %>% as_tibble( rownames = "Vars") %>% janitor::clean_names() %>%   
    mutate( coef = temp$coefficients[,"coef"],
            se = temp$coefficients[,4],
            n = temp$n,
            nevent = temp$nevent,
            HR_CI = paste(format(round(exp_coef, 2), nsmall = 2, scientific = FALSE), " (",
                          format(round(lower_95, 2), nsmall = 2), " to ",
                          format(round(upper_95, 2), nsmall = 2), ")",sep = ""),
            p_value = p) 
  return( Cox_results)
  
}
missingness_func <- function( df){
  
  df %>% summarise( across( everything(), ~ sum( is.na(.))/ n()))
  
  
}
# Import all PRS scores ---------------------------------------------------
Standard_and_Enhanced_PRS <- readRDS("D:/DPhil/UK_Biobank_opioid_application/Add_Basket_4000278/Standard_and_Enhanced_PRS.Rds") %>% janitor::clean_names()

Scaled_Standard_and_Enhanced_PRS <- 
  Standard_and_Enhanced_PRS %>% 
  mutate( across( .cols = !f_eid, ~ scale(.)[ ,1])) %>% 
  janitor::clean_names()

# Import PRS dictionary ---------------------------------------------------
PRS_dictionary <- 
  readxl::read_excel("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Covid_PRS/PRS_dictionary.xlsx", sheet = "Sheet1") %>% 
  mutate( Trait_code = stringr::str_extract(string = Description, pattern = "(?<=\\().*(?=\\))"))

Enhaned_PRS_dictionary <- 
  readxl::read_excel("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Covid_PRS/PRS_dictionary.xlsx", sheet = "Sheet2") %>% 
  mutate( Trait_code = stringr::str_extract(string = Description, pattern = "(?<=\\().*(?=\\))"))
  
Original_source_trait_ICD_10 <- 
  readxl::read_excel("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Covid_PRS/Original_source_trait_ICD_10.xlsx") %>% 
  mutate( Trait = str_replace_all(  Trait, "[^[:alnum:]]", "_")) %>% 
  filter( Trait %in% c( "Atrial_fibrillation", "Coronary_artery_disease", "Ischaemic_stroke", "Venous_thromboembolic_disease"))

distinct_df <- 
  Original_source_trait_ICD_10 %>% 
  group_by( Trait, Field_id) %>% 
  filter( row_number() == 1) %>% 
  ungroup() %>% 
  mutate( Field_id = str_replace_all( Field_id,"[.]", "_")) %>% 
  filter( Trait %in% c( "Atrial_fibrillation", "Coronary_artery_disease", "Ischaemic_stroke", "Venous_thromboembolic_disease"))

Enhanced_distinct_df <- 
  tibble( Trait = c( "Atrial_fibrillation", "Coronary_artery_disease", "Ischaemic_stroke", "Venous_thromboembolic_disease"),
          Field_id = c( "f_26213_0_0", "f_26228_0_0", "f_26249_0_0", "f_26290_0_0"))

# Import baseline raw datasets ---------------------------------------------------------
# UKBB England participants
baseline_df <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID-19_MAIN_df/baseline_df.rds")

table( baseline_df_addtional$PA_status)

# Import Covid-19 testing data --------------------------------------------
#England Covid-19 test results
covid19_result_england <- 
  data.table::fread( input = "D:/DPhil/UK_Biobank_opioid_application/COVID_test_results/Version_20230405/covid19_result_england.txt") %>% 
  mutate( specdate = lubridate::dmy( specdate))

# Import HES --------------------------------------------------------------

# (attention: different HES version)
hesin <- read_delim("D:/DPhil/UK_Biobank_opioid_application/HES/Version_20230405/hesin.txt",
                    delim = "\t",
                    escape_double = FALSE,
                    col_types = cols(epistart = col_date(format = "%d/%m/%Y"),
                                     epiend = col_date(format = "%d/%m/%Y"),
                                     elecdate = col_date(format = "%d/%m/%Y"),
                                     admidate = col_date(format = "%d/%m/%Y"),
                                     disdate = col_date(format = "%d/%m/%Y"),
                                     admisorc = col_character(),
                                     disdest = col_character()),
                    guess_max = 5000,
                    trim_ws = TRUE)

hesin_diag <- read_delim("D:/DPhil/UK_Biobank_opioid_application/HES/Version_20230405/hesin_diag.txt",
                         delim = "\t",
                         escape_double = FALSE,
                         col_types = cols(diag_icd9 = col_character(),
                                          diag_icd9_nb = col_character(),
                                          diag_icd10_nb = col_character()),
                         trim_ws = TRUE)

# Import GP vaccination (diagnosis source) -----------------------------------------
vaccine_diagnosis_source <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Breakthrough_infection/Derived_dataset/vaccine_diagnosis_source.rds")

# Import GP medication and diagnosis data  ----------------------------------------------------------
tpp_gp_medication <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_prescription/tpp_gp_medication.rds")
emis_gp_medication <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_prescription/emis_gp_medication.rds")
combined_gp_medication <- bind_rows( tpp_gp_medication, emis_gp_medication)
all_medication_dmd_code <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_prescription/all_medication_dmd_code.rds")
medication_dictionary_dpa <- readxl::read_excel("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Vac_effect/Derived_datasets/Candidate_medication_list_dpa_v3.xlsx")

emis_gp_clinical_extraction <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_diagnosis/emis_gp_clinical_extraction.rds")
tpp_gp_clinical_extraction <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_diagnosis/tpp_gp_clinical_extraction.rds")
combined_gp_diagnosis <- bind_rows( emis_gp_clinical_extraction, tpp_gp_clinical_extraction)
all_clinical_emis_code <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_diagnosis/all_clinical_emis_code.rds")
all_clinical_tpp_code <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID_19_EMIS_TPP_diagnosis/all_clinical_tpp_code.rds")

# Import death registry ---------------------------------------------------
death <- read_delim("D:\\DPhil\\UK_Biobank_opioid_application\\Death\\Version_20230405\\death.txt", 
                    delim = "\t",
                    col_types = cols(date_of_death = col_date(format = "%d/%m/%Y")))

# Construction cohort -----------------------------------------------------
#coivd-19 comtemporary negative cohort
comtemporary_func <- function(){
  
  set.seed(1)
  
  start.time <- Sys.time()
  # fully_vaccinated cohort
  fully_vaccinated_cohort <- 
    vaccine_diagnosis_source %>% 
    filter( code == "Dose_2") %>% 
    arrange( event_dt) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    filter( event_dt >= as.Date( "2020/03/01"), event_dt <= as.Date( "2021/09/30")) %>% 
    ungroup()
  
  # covid-infection cohort
  covid_infection_cohort_england <- 
    covid19_result_england %>% 
    filter( result == 1) %>% 
    arrange( specdate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    filter( specdate <= as.Date( "2022-09-30")) %>% 
    mutate( index_date = specdate) %>% 
    left_join( select( fully_vaccinated_cohort, eid, event_dt), by = "eid") %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    filter( is.na( date_of_death) | date_of_death > index_date) %>%  # exclude  death
    mutate( BT_index = case_when( event_dt < index_date ~ "BT_infection", TRUE ~ "infection")) %>% 
    select( -date_of_death)
  
  distribution <- covid_infection_cohort_england$index_date
  
  # comtemporary control cohort
  comtemporary_control <- 
    England_participants_ID %>% #all participants from England
    select( eid) %>% 
    filter( !eid %in% covid_infection_cohort_england$eid) %>%   # exclude infected people
    left_join( death, by = "eid") %>% 
    rowwise() %>% 
    mutate( index_date = sample( distribution, 1)) %>%   # assign random index date
    ungroup() %>% 
    filter( is.na( date_of_death) | date_of_death > index_date) %>%  # exclude  death
    mutate( result = 0, BT_index = "comtemporary_control") %>% 
    select( eid, index_date, result, BT_index)
  
  combined_cohort <- 
    covid_infection_cohort_england %>% 
    bind_rows( comtemporary_control)
  
  return( combined_cohort)
  
  print( Sys.time() - start.time)
}
comtemporary_infection_control_cohort <- comtemporary_func()

# Generate historical medication covariates  -----------------------------------------------------
# link other covariates
covariate_infected_generate_func <- function( cohort){
  
  # link to GP medication
  GP_medication_long <- 
    cohort %>% 
    select( eid, index_date) %>% 
    left_join( combined_gp_medication, by = "eid") %>%
    filter( issue_date <= index_date, issue_date >= index_date - 365) %>%  # medication history within years before index date
    left_join( all_medication_dmd_code, by = "dmd_code") %>% 
    left_join( medication_dictionary_dpa, by = c("top_level" = "Top_level_medication_classes")) %>% 
    filter( include == 1)
  
  GP_medication_wide <- 
    GP_medication_long %>% 
    group_by( eid, Manual_category) %>% 
    filter( row_number() == 1) %>% 
    ungroup() %>% 
    mutate( medication = "1") %>% 
    pivot_wider( id_cols = eid,names_from = Manual_category, values_from = medication) %>% 
    janitor::clean_names()
  
  medication_name_list <<-  setdiff( names( GP_medication_wide), "eid") # to environment
    
  return( GP_medication_wide)
  
  
}
covariates_df_list <-  covariate_infected_generate_func( cohort = comtemporary_infection_control_cohort)

# Curate baseline charateriestics -----------------------------------------
link_covariates_func <- function( infection_status){
  
  linked_baseline <- 
    comtemporary_infection_control_cohort %>% 
    filter( BT_index %in% infection_status) %>% 
    select( eid, index_date, BT_index, origin) %>% 
    left_join( Standard_and_Enhanced_PRS, by = c("eid" = "f_eid")) %>% # genetic information
    mutate( across( .cols = !any_of( c("eid", "index_date", "BT_index", "origin")), ~ scale(.)[ ,1])) %>%  
    left_join( lancet_lifestyle_tidy_df , by = c("eid" = "f.eid")) %>% # lifestyle factor
    mutate( across( .cols = starts_with("Tidy_") , ~ case_when(is.na(.) ~ 0, TRUE ~ .))) %>% 
    mutate( Tidy_total_score_CAT = case_when( Tidy_total_score <= 2 ~ "3", 
                                              Tidy_total_score %in% c(3,4) ~ "2", 
                                              Tidy_total_score %in% c(5:9) ~ "1", 
                                              TRUE ~ NA_character_)) %>% 
    mutate( Tidy_total_score_BI = case_when( Tidy_total_score_CAT %in% c( "2", "3") ~ "Healthy", 
                                             Tidy_total_score_CAT %in% c( "1") ~ "Unhealthy", 
                                             TRUE ~ NA_character_)) %>% 
    mutate( Coutinuous_tidy_total_score = -as.numeric(Tidy_total_score)) %>% 
    mutate( across( .cols = starts_with("Tidy_") , ~ as.factor(.))) %>%  
    left_join( baseline_df, by = "eid") %>%  # baseline information
    left_join( baseline_education, by = "eid") %>%  # addtional baseline information
    mutate( recruitment_year = year( recruitment_date) - birth_year,
            age_year = year( index_date) - birth_year,
            IMD = case_when( is.na( IMD) ~ mean( IMD, na.rm = TRUE), TRUE ~ IMD),
            bmi = case_when( is.na( bmi) ~ mean( bmi, na.rm = TRUE), TRUE ~ bmi),
            ethnicity = replace_na( ethnicity, "NA"),
            age_binary = case_when( age_year >= 65 ~ "1", TRUE ~ "0"),
            sex_binary = case_when( sex == "Male" ~ "1", TRUE ~ "0"),
            bmi_binary = case_when( bmi >= 30 ~ "1", TRUE ~ "0") ,
            ethnicity_binary = case_when( ethnicity %in% c( "White") ~ "0", TRUE ~ "1"),
            ethnicity_adj = as.factor(ethnicity_binary),
            infection_type_binary = case_when( BT_index == "infection" ~ "1", TRUE ~ "0"),
            infection_severity_binary = case_when( origin == 1 ~ "1", TRUE ~ "0")) %>% 
    left_join( covariates_df_list, by = "eid") %>% # medication
    mutate( antithrombotics_drug_binary = case_when( antithrombotics == 1 | anticoagulants == 1 ~ "1",TRUE ~ "0")) %>% 
    mutate( AF_quantile= cut_number( f_26212_0_0, n = 5, na.rm = TRUE, labels = FALSE),
            CAD_quantile= cut_number( f_26227_0_0, n = 5, na.rm = TRUE, labels = FALSE),
            ISS_quantile= cut_number( f_26248_0_0, n = 5, na.rm = TRUE, labels = FALSE),
            VTE_quantile= cut_number( f_26289_0_0, n = 5, na.rm = TRUE, labels = FALSE),
            AF_groups = case_when( AF_quantile == "1" ~ "Low", AF_quantile %in% c( "2", "3", "4") ~ "Intermediate", AF_quantile =="5" ~ "High"),
            CAD_groups = case_when( CAD_quantile == "1" ~ "Low", CAD_quantile %in% c( "2", "3", "4") ~ "Intermediate", CAD_quantile =="5" ~ "High"),
            ISS_groups = case_when( ISS_quantile == "1" ~ "Low", ISS_quantile %in% c( "2", "3", "4") ~ "Intermediate", ISS_quantile =="5" ~ "High"),
            VTE_groups = case_when( VTE_quantile == "1" ~ "Low", VTE_quantile %in% c( "2", "3", "4") ~ "Intermediate", VTE_quantile =="5" ~ "High"))
  
  return( linked_baseline)
  
}

linked_infection_baseline <- link_covariates_func( infection_status = c( "BT_infection", "infection"))

covariates_list_infection <- c(  "age_year", "sex","ethnicity_binary", "bmi_binary", "IMD", "Education",
                                 "Tidy_smoking" ,"Tidy_alcohol_intake",  "Tidy_PA", "Tidy_TV_view", "Tidy_sleep_time",    
                                 "Tidy_fruit_vegetable", "Tidy_oil_fish", "Tidy_red_meat",  "Tidy_process_meat" ,  
                                 "Tidy_total_score_CAT")

baseline_table <- tableone::CreateTableOne( data = filter(linked_infection_baseline, !is.na(CAD_groups)), 
                                vars = covariates_list_infection, 
                                includeNA = TRUE,
                                strata = "AF_groups",
                                test = FALSE,
                                addOverall = TRUE) 

table_format <- 
  print(baseline_table, showAllLevels = TRUE, smd = TRUE) %>% 
  as.data.frame() %>% 
  add_row( level = NA, .before = 3) %>% 
  add_row( level = NA, .before = 6) %>% 
  add_row( level = NA, .before = 9) %>% 
  add_row( level = NA, .before = 13) %>% 
  add_row( level = NA, .before = 19)%>% 
  add_row( level = NA, .before = 20)%>% 
  add_row( level = NA, .before = 23)%>% 
  add_row( level = NA, .before = 26)%>% 
  add_row( level = NA, .before = 29)%>% 
  add_row( level = NA, .before = 32)%>% 
  add_row( level = NA, .before = 35)%>% 
  add_row( level = NA, .before = 38)%>% 
  add_row( level = NA, .before = 41)%>% 
  add_row( level = NA, .before = 44)%>% 
  add_row( level = NA, .before = 47)

  
write.table( table_format, "clipboard", sep="\t", row.names=TRUE)
# Curate taget outcomes  -----------------------------------------------
#target disease outcomes
target_outcomes_func <- function( df,  lookback_window, window){
  
  #historical outcomes index
  clinical_sequelae_history <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid") %>% 
    left_join( select( hesin_diag, eid, ins_index, diag_icd10), by = c("eid", "ins_index")) %>% 
    filter( index_date > admidate, admidate >= index_date - lookback_window ) %>% 
    mutate( first_three_charater_diag_icd10 = str_sub( diag_icd10, 1, 3))
  
  match_on_first_three_history <- 
    clinical_sequelae_history %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 3), by = c( "first_three_charater_diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  match_on_first_four_history <- 
    clinical_sequelae_history %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 4), by = c( "diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  clinical_target_outcomes_history <- 
    bind_rows( match_on_first_three_history, match_on_first_four_history) %>% 
    arrange( admidate) %>% 
    group_by( eid, Trait) %>% 
    filter( row_number() == 1) %>% 
    ungroup()
  
  clinical_sequelae_wide_history <- 
    clinical_target_outcomes_history %>% 
    pivot_wider( id_cols = c( eid), names_from = Trait, names_prefix = "history_", values_from = c( admidate))   
  
  #follow-up outcomes index
  clinical_sequelae_after <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid") %>% 
    left_join( select( hesin_diag, eid, ins_index, diag_icd10), by = c("eid", "ins_index")) %>% 
    filter( index_date < admidate ) %>% 
    mutate( first_three_charater_diag_icd10 = str_sub( diag_icd10, 1, 3))
  
  match_on_first_three_after <- 
    clinical_sequelae_after %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 3), by = c( "first_three_charater_diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  match_on_first_four_after <- 
    clinical_sequelae_after %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 4), by = c( "diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  clinical_target_outcomes_after <- 
    bind_rows( match_on_first_three_after, match_on_first_four_after) %>% 
    arrange( admidate) %>% 
    group_by( eid, Trait) %>% 
    filter( row_number() == 1) %>% 
    ungroup()
  
  clinical_sequelae_wide_after <- 
    clinical_target_outcomes_after %>% 
    pivot_wider( id_cols = c( eid), names_from = Trait, names_prefix = "prefix_", values_from = c( admidate))   
  
  names( clinical_sequelae_wide_after)
  
  total_cohort <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    left_join( clinical_sequelae_wide_history, by = "eid") %>% 
    left_join( clinical_sequelae_wide_after, by = "eid") %>% 
    mutate( across( .cols = starts_with("history_"), ~ case_when( !is.na(.x)  ~ 1, TRUE ~ 0), .names = "history_outcome_{.col}")) %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + window, as.Date( "2022-09-30")), 
                                            TRUE ~ pmin( date_of_death, index_date + window, as.Date( "2022-09-30"))),
            across( .cols = starts_with("prefix_"), ~ case_when( !is.na(.x) & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
            across( .cols = starts_with("prefix_"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}"),
            across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x )))  %>% 
    select( -index_date, -date_of_death)
    
  
  return( total_cohort)
  
}

comtemporary_target_outcomes <- target_outcomes_func( df = linked_infection_baseline, 
                                                      lookback_window = 365,
                                                      window = 90)

linked_infection_target_outcomes <- 
  linked_infection_baseline %>%  
  left_join( comtemporary_target_outcomes, by = "eid") %>% 
  filter( !is.na( f_26248_0_0))

pre_omicron_linked_infection_target_outcomes <- 
  linked_infection_target_outcomes %>% 
  filter( index_date <= as.Date("2021/12/01"))

after_omicron_linked_infection_target_outcomes <- 
  linked_infection_target_outcomes %>% 
  filter( index_date > as.Date("2021/12/01"))

#target disease outcomes (primary diagnoses only)
primary_target_outcomes_func <- function( df, lookback_window, window){
  
  #historical outcomes index
  clinical_sequelae_history <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid") %>% 
    left_join( select( hesin_diag, eid, ins_index, diag_icd10), by = c("eid", "ins_index")) %>% 
    filter( index_date > admidate, admidate >= index_date - lookback_window) %>% 
    mutate( first_three_charater_diag_icd10 = str_sub( diag_icd10, 1, 3))
  
  match_on_first_three_history <- 
    clinical_sequelae_history %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 3), by = c( "first_three_charater_diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  match_on_first_four_history <- 
    clinical_sequelae_history %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 4), by = c( "diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  clinical_target_outcomes_history <- 
    bind_rows( match_on_first_three_history, match_on_first_four_history) %>% 
    arrange( admidate) %>% 
    group_by( eid, Trait) %>% 
    filter( row_number() == 1) %>% 
    ungroup()
  
  clinical_sequelae_wide_history <- 
    clinical_target_outcomes_history %>% 
    pivot_wider( id_cols = c( eid), names_from = Trait, names_prefix = "history_", values_from = c( admidate))   
  
  
  #follow-up outcomes index
  clinical_sequelae_after <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( hesin, eid, ins_index, admidate), by = "eid") %>% 
    left_join( select( hesin_diag, eid, level, ins_index, diag_icd10), by = c("eid", "ins_index")) %>% # attention: level == 1
    filter( index_date < admidate, admidate <= index_date + window, level %in% c(1,2,3)) %>% 
    mutate( first_three_charater_diag_icd10 = str_sub( diag_icd10, 1, 3))
  
  match_on_first_three_after <- 
    clinical_sequelae_after %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 3), by = c( "first_three_charater_diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  match_on_first_four_after <- 
    clinical_sequelae_after %>% 
    left_join( filter( Original_source_trait_ICD_10, str_length( ICD_10) == 4), by = c( "diag_icd10" = "ICD_10")) %>% 
    filter( !is.na(Trait))
  
  clinical_target_outcomes_after <- 
    bind_rows( match_on_first_three_after, match_on_first_four_after) %>% 
    arrange( admidate) %>% 
    group_by( eid, Trait) %>% 
    filter( row_number() == 1) %>% 
    ungroup()
  
  clinical_sequelae_wide_after <- 
    clinical_target_outcomes_after %>% 
    pivot_wider( id_cols = c( eid), names_from = Trait, names_prefix = "prefix_", values_from = c( admidate))   
  
  names( clinical_sequelae_wide_after)
  
  total_cohort <- 
    df %>% 
    select( eid, index_date) %>% 
    left_join( select( death, eid, date_of_death), by = "eid") %>% 
    left_join( clinical_sequelae_wide_history, by = "eid") %>% 
    left_join( clinical_sequelae_wide_after, by = "eid") %>% 
    mutate( date_of_death = case_when( date_of_death <  index_date ~ as.Date( "2021-09-30"), 
                                       date_of_death ==  index_date ~ date_of_death + 1, 
                                       TRUE ~ date_of_death)) %>% 
    mutate( across( .cols = starts_with("history_"), ~ case_when( .x < index_date & .x >= index_date - lookback_window ~ 1, TRUE ~ 0), .names = "history_outcome_{.col}")) %>% 
    mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + window, as.Date( "2021-09-30")), 
                                            TRUE ~ pmin( date_of_death, index_date + window, as.Date( "2021-09-30"))),
            across( .cols = starts_with("prefix_"), ~ case_when( .x >= index_date & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
            across( .cols = starts_with("prefix_"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}"),
            across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x ))) %>% 
    select( -index_date, -date_of_death)
  
  
  return( total_cohort)
  
}

comtemporary_primary_target_outcomes <- primary_target_outcomes_func( df = linked_infection_baseline, 
                                                                      lookback_window = 365,
                                                                      window = 90)

linked_infection_primary_target_outcomes <- 
  linked_infection_baseline %>% 
  select( eid, index_date) %>% 
  left_join( Standard_and_Enhanced_PRS, by = c("eid" = "f_eid")) %>% 
  mutate( across( .cols = !any_of( c("eid", "index_date")), ~ scale(.)[ ,1])) %>%  
  left_join( comtemporary_primary_target_outcomes, by = "eid") %>% 
  left_join( select( linked_infection_baseline, eid, index_date, all_of(covariates_list_infection), Coutinuous_tidy_total_score,
                     f_22009_0_1, f_22009_0_2, f_22009_0_3, f_22009_0_4, f_22009_0_5,
                     f_22009_0_6, f_22009_0_7, f_22009_0_8, f_22009_0_9, f_22009_0_10), by = c( "eid", "index_date"))
# Run model (continuous PRS) ---------------------------------------------------------------
#target disease outcomes
Model_target_outcomes_func <- function( df = linked_infection_target_outcomes, 
                                        outcomes = "Atrial_fibrillation", 
                                        exposure = "f_26212_0_0",
                                        prevalent_index = c(0,1)
                                        ){

  FL_Outcomes <- outcomes
  PRS_exposure <- exposure
  
  complete_df <- 
    df %>%
    filter( .data[[paste( "history_outcome_history_", FL_Outcomes, sep = "")]] %in% prevalent_index) 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  output_HR <- survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() == 1)
  
  follow_up_days <- paste( "follow_up_days_prefix_", FL_Outcomes, sep = "")
  follow_up_events <- paste( "incident_outcome_prefix_", FL_Outcomes, sep = "")
  
  output_summary_statistics <- 
    complete_df %>% 
    summarise( sum_follow = sum( .data[[follow_up_days]]),
               cases = sum( .data[[follow_up_events]]),
               rate = cases/ sum_follow * 1000 * 365) 
  
  output <- 
    bind_cols( output_HR, output_summary_statistics) %>% 
    mutate( cumulative_rate = (cases/ n) * 100)
  
  return( output)
  
}

Continuous_HR_among_infection <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_target_outcomes_func,
           df = linked_infection_target_outcomes,
           prevalent_index = c(0,1),
           .id = "target_outcome")

pre_omicron_Continuous_HR_among_infection <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_target_outcomes_func,
           df = pre_omicron_linked_infection_target_outcomes,
           prevalent_index = c(0,1),
           .id = "target_outcome")

after_omicron_Continuous_HR_among_infection <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_target_outcomes_func,
           df = after_omicron_linked_infection_target_outcomes,
           prevalent_index = c(0,1),
           .id = "target_outcome")

write.table( after_omicron_Continuous_HR_among_infection, "clipboard", sep="\t", row.names=FALSE)

# advanced approach to calcualte coef and se in subgroups
Update_model_target_outcomes_subgroup_func <- function( df = linked_infection_target_outcomes, 
                                                        outcomes = "Atrial_fibrillation", 
                                                        exposure = "f_26212_0_0"){
  FL_Outcomes <- outcomes
  PRS_exposure <- exposure
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = "")) 
  
  temp_cox_func <- function( var = "age_binary"){
    
    interaction_formula <- 
      as.formula( paste( 
        paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
        paste( c( PRS_exposure, paste( PRS_exposure, var, sep = ":"),
                  "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                  "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
        sep = ""))
    
    running_df <- 
      df %>% 
      rename( "subgroup_var" = all_of(var))
    
    follow_up_days <- paste( "follow_up_days_prefix_", FL_Outcomes, sep = "")
    follow_up_events <- paste( "incident_outcome_prefix_", FL_Outcomes, sep = "")
    
    output_summary_statistics <- 
      running_df %>% 
      group_by( subgroup_var) %>% 
      summarise( sum_follow = sum( .data[[follow_up_days]]),
                 cases = sum( .data[[follow_up_events]]),
                 rate = (cases/ sum_follow) * 1000 * 365,
                 n = n()) %>% 
      ungroup()
    
    subgroup_model <- survival::coxph( formula = interaction_formula, data = df) 
    p <- subgroup_model %>% temp_func( ) %>% filter( row_number() == n())
    
    coef <- summary(subgroup_model)$coef
    
    # Coefficients for subgroup A (reference level)
    coef_A <- coef[, "coef"][c(PRS_exposure)]
    
    # Coefficients for subgroup B
    coef_B <- coef[, "coef"][c(PRS_exposure)] + coef[, "coef"][c(paste( PRS_exposure, ":", var, "1", sep = ""))]
    
    
    # Variance for subgroup A (reference level)
    vcov_matrix <- vcov(subgroup_model)
    var_A <- diag(vcov_matrix)[c(PRS_exposure)] %>% sqrt()
    
    # Variance for subgroup B
    var_B <- var_A + 
      2 * vcov_matrix[PRS_exposure, c(paste( PRS_exposure, ":", var, "1", sep = ""))] +
      diag(vcov_matrix)[paste( PRS_exposure, ":", var, "1", sep = "")] %>% sqrt()
    
    # Confidence intervals for subgroup A
    ci_A_lower <- coef_A - 1.96 * var_A
    ci_A_upper <- coef_A + 1.96 * var_A
    ci_A <- tibble(group = "0", point = coef_A, Lower = ci_A_lower, Upper = ci_A_upper)
    
    # Confidence intervals for subgroup B
    ci_B_lower <- coef_B - 1.96 * var_B
    ci_B_upper <- coef_B + 1.96 * var_B
    ci_B <- tibble(group = "1", point = coef_B, Lower = ci_B_lower, Upper = ci_B_upper)
    
    combined <- bind_rows( ci_A, ci_B)
    
    output <- 
      output_summary_statistics %>% 
      left_join( combined, by = c("subgroup_var" = "group")) %>% 
      mutate( across( c(point, Lower, Upper), ~ exp(.), .names = "Exp_{.col}")) %>% 
      mutate( cumulative_rate = (cases/ n) * 100) %>% 
      mutate( p_interaction_value = p$p_value)
    
    return( output)
    
  }
  
  varname_list <- 
    c( "age_binary", "sex_binary", "bmi_binary", "ethnicity_binary", "infection_type_binary", "infection_severity_binary", "antithrombotics_drug_binary") %>% 
    set_names()

  subgroup_HR <- 
    map_df( varname_list, temp_cox_func, .id = "group") %>% mutate( group_name = str_c( group, "_",subgroup_var)) %>% 
    rename( exp_coef = Exp_point, lower_95 = Exp_Lower, upper_95 = Exp_Upper)
  
  overall_HR <- 
    survival::coxph( formula = formula_input, data = df) %>% temp_func() %>% 
    filter( row_number() == 1) %>% 
    mutate( group_name = "overll") 
  
  
  output <- bind_rows( subgroup_HR, overall_HR)
  
  return( output)
  
  
}


Continuous_HR_among_infection_subgroup <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Update_model_target_outcomes_subgroup_func,
           df = linked_infection_target_outcomes,
           .id = "target_outcome")

subgroup_forest_plot <- function( df = Continuous_HR_among_infection_subgroup, 
                                  outcome = "Atrial_fibrillation"){
  
  library(forestplot)
  
  temp_df <- 
    df %>% 
    mutate( group_name = factor( group_name, 
                                 levels = c( "age_binary_1", "age_binary_0", 
                                             "sex_binary_1", "sex_binary_0", 
                                             "bmi_binary_1", "bmi_binary_0", 
                                             "ethnicity_binary_1", "ethnicity_binary_0",
                                             "antithrombotics_drug_binary_1", "antithrombotics_drug_binary_0",
                                             "infection_severity_binary_1", "infection_severity_binary_0",
                                             "infection_type_binary_1", "infection_type_binary_0",
                                             "overll"),
                                 labels = c( ">= 65 years", "< 65 years", 
                                             "Male", "Female", 
                                             ">= 30", "< 30", 
                                             "White", "Other ethnic",
                                             "User", "Non-user",
                                             "Inpatient", "Outpatient",
                                             "Not or partial", "Complete",
                                             "Overall")))  %>% 
    arrange( group_name) %>%
    mutate( characteristics = case_when( group_name == ">= 65 years" ~ "Age",
                                         group_name == "Male" ~ "Sex",
                                         group_name == ">= 30" ~ "BMI",
                                         group_name == ">= 65 years" ~ "Age",
                                         group_name == "White" ~ "Ethnicity",
                                         group_name == "User" ~ "Antithrombotics",
                                         group_name == "Inpatient" ~ "PCR test settings",
                                         group_name == "Not or partial" ~ "Vaccincation")) %>% 
    filter( target_outcome == outcome) %>% 
    mutate( HR_label = format(round(exp_coef, 2), nsmall = 2),
            n = as.character(n),
            nevent = as.character(nevent),
            rate = format(round(rate, 2), nsmall = 2),
            p_interaction_value = format(round(p_interaction_value, 2), nsmall = 2),
            p_interaction_value = case_when( p_interaction_value == "0.00" ~ "<0.01", TRUE ~ p_interaction_value)) %>% 
    rename( mean = exp_coef,
            lower = lower_95,
            upper = upper_95,
            study = group_name,
            number = n,
            case = nevent,
            incidence = rate,
            OR = HR_label) 
    
  
  Continuous_HR_among_infection_subgroup$group_name
    
  base_data <- 
    temp_df %>% 
    filter( study != "Overall") %>% 
    select( characteristics, study,  cases, incidence, OR, p_interaction_value, mean, lower, upper) %>% 
    mutate( p_interaction_value = case_when( is.na(characteristics) ~ NA_character_, TRUE ~ p_interaction_value),
            cases = as.character(cases))

  summary <- 
    temp_df %>% 
    filter( study == "Overall") %>%
    select( study, OR, mean, lower, upper) %>% 
    mutate( summary = TRUE)
  
  header <- tibble(characteristics = c( "", "Characteristics"),
                   study = c( "", "Subgroup"),
                   #number = c( "No. of ", "persons"),
                   cases = c( "No. of ", "events"),
                   incidence = c("Incidence", "rate"),
                   OR = c( "Hazard", "ratio"),
                   p_interaction_value = c( "P", "interaction"),
                   summary = TRUE)
  
  cochrane_output_df <- bind_rows(header,
                                  base_data,
                                  summary)
  
  cochrane_output_df %>% 
    forestplot(labeltext = c(characteristics, study, cases, incidence, OR, p_interaction_value), 
               is.summary = summary,
               
               clip = c(0.5, 3), 
               xticks = log(c( 0.5, 1, 2, 3)),
               xlog = TRUE, 
               boxsize = 0.25,
               
               graphwidth = unit(25, "mm"),
               graph.pos = 5,
               colgap = unit(4, "mm"),
               #row property
               lineheight = unit(6, "mm"),
               hrzl_lines = list("1" = gpar(lty=1, lwd=0.5, columns=c(1:7)),
                                 "3" = gpar(lty=1, lwd=0.5, columns=c(1:7)),
                                 # "4" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "5" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "6" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "7" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "8" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "9" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "10" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "11" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "12" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "13" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "14" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "15" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 # "16" = gpar(lty=1, lwd=0.4, columns=c(1:5), col = "grey"),
                                 "17" = gpar(lty=1, lwd=0.5, columns=c(1:4, 7))),
               
               col=fpColors(box=c("black"), line=c("black"), zero = "black"),
               
               txt_gp = fpTxtGp( label = gpar( fontfamily = "sans", fontsize = 8),
                                 ticks = gpar( fontfamily = "sans", fontsize = 10)),
               new_page = TRUE
               )
  
  
  
}

subgroup_AF <- subgroup_forest_plot( outcome = "Atrial_fibrillation")
subgroup_CAD <- subgroup_forest_plot( outcome = "Coronary_artery_disease")
subgroup_VTE <- subgroup_forest_plot( outcome = "Venous_thromboembolic_disease")
subgroup_ISS <- subgroup_forest_plot( outcome = "Ischaemic_stroke")

pdf(file = "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Manuscripts/COVID_19_PRS/Submission/NC/Second_round/Plots/Forest_ISS.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 3.5)

subgroup_ISS

dev.off()


# Run model (categorical PRS) ------------------------------------------------
# three categories
Model_three_categories <- function( df, outcomes, exposure){
  
  
  PRS_exposure <- exposure 
  FL_Outcomes <- outcomes 
  
  complete_df <- 
    df %>%
    mutate( across( .cols = starts_with("f_"), ~ as.factor( case_when( . > quantile( ., 0.8, na.rm = TRUE) ~ "3", 
                                                                       . > quantile( ., 0.2, na.rm = TRUE) & . <= quantile( ., 0.8, na.rm = TRUE) ~ "2",
                                                                       . <= quantile( ., 0.2, na.rm = TRUE) ~ "1", TRUE ~ NA_character_)))) 
  
  # run model
  adj_formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  output <- survival::coxph( formula = adj_formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <= 2)
  
  return( output)
  
  
}
Model_three_categories_statistics <- function( df,outcomes, exposure){
  
  PRS_exposure <- exposure 
  FL_Outcomes <- outcomes 
  
  complete_df <- 
    df %>%
    mutate( across( .cols = starts_with("f_"), ~ as.factor( case_when( . > quantile( ., 0.66, na.rm = TRUE) ~ "3", 
                                                                       . > quantile( ., 0.33, na.rm = TRUE) & . <= quantile( ., 0.66, na.rm = TRUE) ~ "2",
                                                                       . <= quantile( ., 0.33, na.rm = TRUE) ~ "1", TRUE ~ NA_character_)))) 
  
  output <- 
    complete_df %>% 
    group_by( !!sym(PRS_exposure)) %>% 
    summarise( n = n(), case = sum( .data[[paste( "incident_outcome_prefix_", FL_Outcomes, sep = "")]]))
  
  return(output)
}

# all cohort
HR_three_categories <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_three_categories,
           df = linked_infection_target_outcomes,
           .id = "target_outcome") 

HR_three_categories_summary_statistics <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
          distinct_df$Field_id %>% purrr::set_names(),
          Model_three_categories_statistics,
          df = linked_infection_target_outcomes,
          .id = "target_outcome")

write.table( HR_three_categories_summary_statistics, "clipboard", sep="\t", row.names=TRUE)

# pre omicron cohort
pre_omicron_HR_three_categories <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_three_categories,
           df = pre_omicron_linked_infection_target_outcomes,
           .id = "target_outcome")

# after omicron cohort
after_omicron_HR_three_categories <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_three_categories,
           df = after_omicron_linked_infection_target_outcomes,
           .id = "target_outcome") 

write.table( after_omicron_HR_three_categories, "clipboard", sep="\t", row.names=TRUE)

# Run model (lifestyle) ---------------------------------------------------------------
# estimate HR (continuous lifestyle score)
Model_continous_lifestyle_func <- function( df = linked_infection_target_outcomes, 
                                            outcomes = "Ischaemic_stroke", 
                                            exposure = "f_26248_0_0"){
  
  PRS_exposure <- exposure 
  FL_Outcomes <- outcomes 
  
  complete_df <- 
    df %>%
    as.data.frame() # extremely important step for using functions in the survminer package
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( "Coutinuous_tidy_total_score", "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  all_PRS <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() == 1)
  
  return( all_PRS)
}

Continous_lifestyle_HR <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_continous_lifestyle_func,
           df = linked_infection_target_outcomes,
           .id = "target_outcome")

pre_omicron_continous_lifestyle_HR <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_continous_lifestyle_func,
           df = pre_omicron_linked_infection_target_outcomes,
           .id = "target_outcome")

after_omicron_continous_lifestyle_HR <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_continous_lifestyle_func,
           df = after_omicron_linked_infection_target_outcomes,
           .id = "target_outcome")

write.table( after_omicron_continous_lifestyle_HR, "clipboard", sep="\t", row.names=TRUE)

# estimate HR (Unfavourable vs Intermediate vs favourable)
Model_lifestyle_category_func <- function( df = linked_infection_target_outcomes, 
                                           outcomes = "Ischaemic_stroke", 
                                           exposure = "Tidy_total_score_CAT"){
  
  FL_Outcomes <- outcomes 
  
  complete_df <- df 

  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( exposure, "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  all_PRS <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=2)
  
  return( all_PRS)
  

  
}

# all cohort
lifestyle_category_HR <- 
  map_df( c("Atrial_fibrillation", "Coronary_artery_disease", "Venous_thromboembolic_disease", "Ischaemic_stroke") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = linked_infection_target_outcomes,
          .id = "target_outcome") 

# pre omicron cohort
pre_omicron_lifestyle_category_HR <- 
  map_df( c("Atrial_fibrillation", "Coronary_artery_disease", "Venous_thromboembolic_disease", "Ischaemic_stroke") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = pre_omicron_linked_infection_target_outcomes,
          .id = "target_outcome")

# after omicron cohort
after_omicron_lifestyle_category_HR <- 
  map_df( c("Atrial_fibrillation", "Coronary_artery_disease", "Venous_thromboembolic_disease", "Ischaemic_stroke") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = after_omicron_linked_infection_target_outcomes,
          .id = "target_outcome")

write.table( after_omicron_lifestyle_category_HR, "clipboard", sep="\t", row.names=TRUE)
# Plot KM -----------------------------------------------------------------
#  KM curve for PRS
KM_PRS_plot_func <- function( df, outcomes, exposure, manual_label, y_limit = 0.05 ){
  
  FL_Outcomes <- outcomes 
  PRS_exposure <- exposure
  

  complete_df <- 
    df %>%
    mutate( across( .cols = starts_with("f_"), ~ as.factor( case_when( . > quantile( ., 0.66, na.rm = TRUE) ~ "3", 
                                                                       . > quantile( ., 0.33, na.rm = TRUE) & . <= quantile( ., 0.66, na.rm = TRUE) ~ "2",
                                                                       . <= quantile( ., 0.33, na.rm = TRUE) ~ "1", TRUE ~ NA_character_)))) 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure), collapse="+"),
      sep = ""))
  
  KM_model <- survminer::surv_fit( formula_input, data = complete_df)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    #ylim = c(0, y_limit),
                                    xlim = c(0, 91),
                                    #ncensor.plot = TRUE,
                                    fun = "event",
                                    xlab = "Days after the infection",
                                    ylab = paste("Cumulative incidence for", manual_label),
                                    #break.time.by = 15,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#00a087FF", "#DF8F44FF", "#DC0000FF"),
                                    legend.title = "",
                                    #legend.labs = c( "Low (reference)", paste("Intermediate: HR", Intermediate), paste("High: HR", High)),
                                    legend.labs = c( "Low genetic risk", paste("Intermediate genetic risk"), paste("High genetic risk")),
                                    # tables.y.text.col = TRUE,
                                    # tables.y.text = FALSE,
                                    # fontsize = unit( 2.5, "cm"),
                                    # tables.height = 0.2,
                                    size = 0.7,
                                    surv.linewidth = 0.2,
                                    # ggtheme = theme_classic()
  )
  
  p1  <-  
    KM_plot$plot + 
    scale_x_continuous( expand = c(0, 0), limits = c(0, 90), breaks = seq(0, 90, by = 30))+
    scale_y_continuous( expand = c(0, 0), limits  = c(0, y_limit))+
    theme_classic() +
    guides(color = guide_legend(direction = "vertical"), fill = FALSE)+
    theme( text = element_text( family = "sans", color = "black"), 
           #axis.title = element_blank(),
           panel.grid.minor.y = element_blank(),
           panel.grid.minor.x = element_blank(),
           legend.background = element_rect(fill='transparent'),
           legend.position = c(0.3, 0.8))
  
  # p2 = KM_plot$table + theme( plot.title=element_text(size=10),
  #                             text = element_text( family = "sans"),
  #                             axis.text = element_text( family = "sans", color = "black", size = 12),
  #                             axis.text.x=element_blank(),
  #                             axis.ticks.x=element_blank()) + labs( x = "")
  # 
  # p3 = KM_plot$cumevents + theme( plot.title=element_text(size=10),
  #                                 text = element_text( family = "sans", size = 6),
  #                                 axis.text = element_text( family = "sans", color = "black"),
  #                                 axis.text.x=element_blank(),
  #                                 axis.ticks.x=element_blank()) + labs( x = "")
  
  # combined <- cowplot::plot_grid( p1, p2, p3, align = "v", ncol =1,rel_heights = c(2.5,1,1))
  
  return( p1)
  
  
} 

KM_PRS_among_infection <- 
  pmap( list(distinct_df$Trait %>% purrr::set_names(),
             distinct_df$Field_id %>% purrr::set_names(),
             c("AF", "CAD", "ISS", "VTE"),
             c(0.03, 0.03, 0.005,0.005)),
        KM_PRS_plot_func,
        df = linked_infection_target_outcomes)

KM_PRS_among_infection$Coronary_artery_disease


pdf(file = "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Manuscripts/COVID_19_PRS/Submission/NC/Second_round/Plots/PRS_KM.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 11)

KM_PRS_among_infection$Atrial_fibrillation + 
  KM_PRS_among_infection$Coronary_artery_disease + 
  KM_PRS_among_infection$Venous_thromboembolic_disease + 
  KM_PRS_among_infection$Ischaemic_stroke + 
  plot_layout( ncol = 1) 


dev.off()


#  KM curve for lifestyle
KM_lifestyle_plot_func <- function( df, outcomes, exposure, manual_label, y_limit = 0.05 ){
  
  PRS_exposure <- exposure 
  FL_Outcomes <- outcomes  
  
  complete_df <- 
    df %>%
    mutate( across( .cols = starts_with("f_"), ~ as.factor( case_when( . > quantile( ., 0.8, na.rm = TRUE) ~ "3", 
                                                                       . > quantile( ., 0.2, na.rm = TRUE) & . <= quantile( ., 0.8, na.rm = TRUE) ~ "2",
                                                                       . <= quantile( ., 0.2, na.rm = TRUE) ~ "1", TRUE ~ NA_character_)))) 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure), collapse="+"),
      sep = ""))
  
  KM_model <- survminer::surv_fit( formula_input, data = complete_df)
  
  KM_plot <- survminer::ggsurvplot( KM_model,
                                    conf.int = TRUE,
                                    risk.table = TRUE,
                                    cumevents = TRUE,
                                    xlim = c(0, 91),
                                    fun = "event",
                                    xlab = "Days after the infection",
                                    ylab = paste("Cumulative incidence for", manual_label),
                                    #break.time.by = 15,
                                    censor.shape = 124,
                                    censor = FALSE,
                                    palette = c("#DC0000FF", "#DF8F44FF", "#00a087FF" ),
                                    legend.title = "",
                                    #legend.labs = c( "Low (reference)", paste("Intermediate: HR", Intermediate), paste("High: HR", High)),
                                    legend.labs = c( "Unfavorable lifestyle", paste("Moderate lifestyle"), paste("Favorable lifestyle")),
                                    # tables.y.text.col = TRUE,
                                    # tables.y.text = FALSE,
                                    # fontsize = unit( 2.5, "cm"),
                                    # tables.height = 0.2,
                                    size = 0.7
                                    # ggtheme = theme_classic()
  )
  
  p1  <-  
    KM_plot$plot + 
    scale_x_continuous( expand = c(0, 0), limits = c(0, 90), breaks = seq(0, 90, by = 30))+
    scale_y_continuous( expand = c(0, 0), limits  = c(0, y_limit))+
    theme_classic() +
    guides(color = guide_legend(direction = "vertical"), fill = FALSE)+
    theme( text = element_text( family = "sans", color = "black"), 
           #axis.title = element_blank(),
           panel.grid.minor.y = element_blank(),
           panel.grid.minor.x = element_blank(),
           legend.background = element_rect(fill='transparent'),
           legend.position = c(0.3, 0.8))
  
  
  return( p1)
  
  
} 

KM_lifestyle_among_infection <- 
  pmap( list(distinct_df$Trait %>% purrr::set_names(),
             exposure = "Tidy_total_score_CAT",
             c("AF", "CAD", "ISS", "VTE"),
             c(0.03, 0.03, 0.005,0.005)),
        KM_lifestyle_plot_func,
        df = linked_infection_target_outcomes)

KM_lifestyle_among_infection$Atrial_fibrillation

pdf(file = "D:/Document_sync/OneDrive - Nexus365/DPhil_training/Manuscripts/COVID_19_PRS/Submission/NC/Second_round/Plots/Lifestyle_KM.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 11)

KM_lifestyle_among_infection$Atrial_fibrillation + 
  KM_lifestyle_among_infection$Coronary_artery_disease + 
  KM_lifestyle_among_infection$Venous_thromboembolic_disease + 
  KM_lifestyle_among_infection$Ischaemic_stroke + 
  plot_annotation( ) +
  plot_layout( ncol = 1) 

dev.off()  

# Sensitivity analyses (Enhanced PRS) ----------------------------------------------------
Enhanced_PRS_outcomes_func <- function( df = linked_infection_target_outcomes, outcomes = "Atrial_fibrillation", exposure = "f_26213_0_0"){
  
  FL_Outcomes <- outcomes
  PRS_exposure <- exposure
  
  complete_df <- df 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "f_22009_0_1", "f_22009_0_2", "f_22009_0_3"), collapse="+"),
      sep = ""))
  
  output_HR <- survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() == 1)
  
  follow_up_days <- paste( "follow_up_days_prefix_", FL_Outcomes, sep = "")
  follow_up_events <- paste( "incident_outcome_prefix_", FL_Outcomes, sep = "")
  
  output_summary_statistics <- 
    df %>% 
    summarise( sum_follow = sum( .data[[follow_up_days]]),
               cases = sum( .data[[follow_up_events]]),
               rate = cases/ sum_follow * 1000 * 365) 
  
  output <- 
    bind_cols( output_HR, output_summary_statistics) %>% 
    mutate( cumulative_rate = (cases/ n) * 100)
  
  return( output)
  
}

Enhanced_continuous_HR_among_infection <- 
  map2_df( Enhanced_distinct_df$Trait %>% purrr::set_names(),
           Enhanced_distinct_df$Field_id %>% purrr::set_names(),
           Enhanced_PRS_outcomes_func,
           df = linked_infection_target_outcomes,
           .id = "target_outcome")

# Sensitivity analyses (incident outcome only)
# incident outcomes
Model_target_incident_func <- function( df, 
                                        outcomes, 
                                        exposure, 
                                        prevalent_index){
  
  FL_Outcomes <- outcomes
  PRS_exposure <- exposure
  
  complete_df <- 
    df %>%
    filter( .data[[paste( "history_outcome_history_", FL_Outcomes, sep = "")]] == prevalent_index) 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "ethnicity_binary", "bmi_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9"), collapse="+"),
      sep = ""))
  
  output_HR <- survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() == 1)
  
  return( output_HR)
  
}

Incident_continuous_HR_among_infection <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_target_incident_func,
           df = linked_infection_target_outcomes,
           prevalent_index = 0,
           .id = "target_outcome")

Model_continous_incident_lifestyle_func <- function( df = linked_infection_target_outcomes, 
                                            outcomes = "Ischaemic_stroke", 
                                            exposure = "f_26248_0_0",
                                            prevalent_index){
  
  PRS_exposure <- exposure 
  FL_Outcomes <- outcomes 
  
  complete_df <- 
    df %>%
    filter( .data[[paste( "history_outcome_history_", FL_Outcomes, sep = "")]] == prevalent_index) 
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_prefix_", FL_Outcomes, ",", "incident_outcome_prefix_", FL_Outcomes, ")", "~", sep = ""), 
      paste( c( "Coutinuous_tidy_total_score", "age_year", "sex", "ethnicity_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  all_PRS <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() == 1)
  
  return( all_PRS)
}


Incident_continous_lifestyle_HR <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_continous_incident_lifestyle_func,
           df = linked_infection_target_outcomes,
           prevalent_index = 0,
           .id = "target_outcome")


# Sensitivity analyses (primary diagnosis only)
Continuous_primary_HR_among_infection <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_target_outcomes_func,
           df = linked_infection_primary_target_outcomes,
           .id = "target_outcome")


Sensivity_continous_lifestyle_HR <- 
  map2_df( distinct_df$Trait %>% purrr::set_names(),
           distinct_df$Field_id %>% purrr::set_names(),
           Model_continous_lifestyle_func,
           df = linked_infection_primary_target_outcomes,
           .id = "target_outcome")

all_diabetes_cases <- 
  hesin_diag %>% 
  filter( str_detect( diag_icd10, "E119")) %>% 
  left_join( select( hesin, eid, ins_index, admidate), by = c(  "eid", "ins_index"))

table( all_diabetes_cases$diag_icd10)

all_diabetes_first_occurence_cases <- 
  all_diabetes_cases %>% 
  arrange( admidate) %>% 
  group_by( eid) %>% 
  filter( row_number() == 1)

baseline_df

temp <- 
  linked_infection_target_outcomes %>% 
  select( eid, index_date, f_26212_0_0, f_26227_0_0, f_26248_0_0, f_26289_0_0, Tidy_total_score) %>% 
  left_join( baseline_df, by = c("eid")) %>% #baseline charateristics
  left_join( baseline_education, by = c("eid")) %>% 
  left_join( select( death, eid, date_of_death), by = c( "eid" = "eid")) %>% 
  left_join( select( all_diabetes_first_occurence_cases, eid, admidate), by = c( "eid" = "eid")) %>% 
  mutate( Tidy_total_score = -as.numeric( Tidy_total_score)) %>% 
  mutate( recruitment_year = year( recruitment_date) - birth_year,
          age_year = year( index_date) - birth_year,
          IMD = case_when( is.na( IMD) ~ mean( IMD, na.rm = TRUE), TRUE ~ IMD),
          bmi = case_when( is.na( bmi) ~ mean( bmi, na.rm = TRUE), TRUE ~ bmi),
          ethnicity = replace_na( ethnicity, "NA"),
          age_binary = case_when( age_year >= 65 ~ "1", TRUE ~ "0"),
          sex_binary = case_when( sex == "Male" ~ "1", TRUE ~ "0"),
          bmi_binary = case_when( bmi >= 30 ~ "1", TRUE ~ "0") ,
          ethnicity_binary = case_when( ethnicity %in% c( "White") ~ "1", TRUE ~ "0")) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + 90, as.Date( "2022-09-30")), 
                                          TRUE ~ pmin( date_of_death, index_date + 90, as.Date( "2022-09-30"))),
          outcome_index = case_when( !is.na( admidate) & admidate >= index_date & admidate <= index_date + 90 ~ 1, TRUE ~ 0),
          follow_up_days = follow_up_end_date - index_date) 

NCO_func <- function( exposure){
  
  PRS_exposure <- exposure
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "survival::Surv( follow_up_days", ",", "outcome_index", ")", "~", sep = ""), 
      paste( c( PRS_exposure, 
                "age_year", "sex", "ethnicity_binary", "bmi_binary", "IMD", "Education", 
                "f_22009_0_1", "f_22009_0_2", "f_22009_0_3", "f_22009_0_4", "f_22009_0_5", "f_22009_0_6", "f_22009_0_7", "f_22009_0_8", "f_22009_0_9", "f_22009_0_10"), collapse="+"),
      sep = ""))
  
  output_HR <- survival::coxph( formula = formula_input, data = temp) %>% temp_func() %>% filter( row_number() == 1)
  
  
  
}

NCO_HR_among_infection <- map_df( distinct_df$Field_id %>% purrr::set_names(distinct_df$Trait), 
                                  NCO_func, 
                                  .id = "target_outcome")







