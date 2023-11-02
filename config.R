require("svMisc")
require("ggplot2")
require("dplyr")
require("tidyr")
require("reticulate")
library("stringr")
library("jsonlite")
library("RColorBrewer")
library("knitr")
library("doParallel")
library("foreach")
library("gtsummary")
library("labelled")

options(warn=-1)
# Only show code in HTML output, not in word documents for example
if (!knitr::is_html_output()) knitr::opts_chunk$set(echo=FALSE)
opts_knit$set(eval.after = "fig.cap")


# Function to change label colors of CONSORT plot
change_label_color <- function(g, n, txt_color, box_color){
    for(i in 1:n){
        g[["children"]][[1]][["children"]][[1]][["children"]][[1]][["children"]][["label"]][["children"]][[i]][["txt_gp"]][["col"]] = txt_color
        g[["children"]][[1]][["children"]][[1]][["children"]][[1]][["children"]][["label"]][["children"]][[i]][["box_gp"]][["fill"]] = box_color
    }
    return(g)
}
# For parallel computing
num_cores = detectCores() - 1
cluster = makeCluster(num_cores)

# Parameters
SC_MRI_MAX_BEFORE = 180 # Maximal number of days that spinal cord MRI may have taken place before treatment initiation
SC_MRI_MAX_AFTER = 180 # Maximal number of days that spinal cord MRI may have taken place after treatment initiation
EDSS_MAX_BEFORE = 90 # Maximal number of days that EDSS was determined before baseline SC MRI
EDSS_MAX_AFTER = 90 # Maximal number of days that EDSS was determined after baseline SC MRI
BRAIN_MRI_MAX_BEFORE = 90 # Maximal number of days that baseline brain MRI took place before baseline SC MRI
BRAIN_MRI_MAX_AFTER = 90 # Maximal number of days that baseline brain MRI took place after baseline SC MRI
MIN_INTERVAL = 180 # Minimal interval between baseline SC MRI and the next SC MRI to take into account
HDMT_CUTOFF = 0.95 # Proportion of follow-up that patient needed to have used a HDMT to be classified as part of the HDMT-group
LDMT_CUTOFF = 0.95 # Proportion of follow-up that patient needed to have used a LDMT to be classified as part of the LDMT-group
SENSITIVITY = TRUE # If current set of parameters are for a sensitivity analysis.

BIO_EFFECTIVENESS = list("ALE" = 5 * 365, "OCR" = 0.5*365, "RTX" = 0.5*365, "CLA" = 96 * 7, "NAT" = 60, "OFA" = 0.5 * 365, "MIT" = 0.5 * 365) # Biological effectiveness of DMTs in days

hdmt = c("NAT", "OCR", "MIT", "RTX", "CLA", "ALE", "OFA") # High efficacy-therapies
idmt = c("DMF", "SIP", "FIN", "PON", "OZA") # Intermediate efficacy-therapies
ldmt = c("GLA", "MTX", "AZA", "IVIG", "IFN", "TER") # Low efficacy-therapies
allowed = c("FAM", "PRED", "PLAS", "VIT") # Allowed add-on therapies

### Helper functions ###
# Function for adding pre-matching SMDs for included variables
unmatched_smd <- function(data, variable, by, ...){
    if(variable == 'GENDER'){
        c(tableone.all['`factor(GENDER)`M', 'Std..Mean.Diff.'], tableone.all['`factor(GENDER)`F', 'Std..Mean.Diff.'])
    }else if(variable == 'BRAIN_MRI_T2_LESIONCAT_LAB'){
        c(tableone.all['`factor(BRAIN_MRI_T2_LESIONCAT)`0', 'Std..Mean.Diff.'], tableone.all['`factor(BRAIN_MRI_T2_LESIONCAT)`1', 'Std..Mean.Diff.'], tableone.all['`factor(BRAIN_MRI_T2_LESIONCAT)`2', 'Std..Mean.Diff.'], tableone.all['`factor(BRAIN_MRI_T2_LESIONCAT)`3', 'Std..Mean.Diff.'])
    }else if(variable == 'COUNTRY'){
        tableone.all[grepl('COUNTRY', rownames(tableone.all)), 'Std..Mean.Diff.']
    }else{
        tableone.all[variable, 'Std..Mean.Diff.']
    }
}
# Function for adding post-matching SMDs for included variables
matched_smd <- function(data, variable, by, ...){
    if(variable == 'GENDER'){
        c(tableone.matched['`factor(GENDER)`M', 'Std..Mean.Diff.'], tableone.matched['`factor(GENDER)`F', 'Std..Mean.Diff.'])
    }else if(variable == 'BRAIN_MRI_T2_LESIONCAT_LAB'){
        c(tableone.matched['`factor(BRAIN_MRI_T2_LESIONCAT)`0', 'Std..Mean.Diff.'], tableone.matched['`factor(BRAIN_MRI_T2_LESIONCAT)`1', 'Std..Mean.Diff.'], tableone.matched['`factor(BRAIN_MRI_T2_LESIONCAT)`2', 'Std..Mean.Diff.'], tableone.matched['`factor(BRAIN_MRI_T2_LESIONCAT)`3', 'Std..Mean.Diff.'])
    }else if(variable == 'COUNTRY'){
        tableone.matched[grepl('COUNTRY', rownames(tableone.matched)), 'Std..Mean.Diff.']
    }else{
        tableone.matched[variable, 'Std..Mean.Diff.']
    }
}
# Generate main patient characteristics table
generate_table1 <- function(matched, matched_data, unmatched_data){
    # Refactor/recode variables for table
    matched_data$HIGH_GROUP = recode(as.factor(matched_data$HIGH_GROUP), '1'='Matched hDMT', '0'='Matched lDMT')
    matched_data = subset(matched_data, select = -c(subclass) )
    unmatched_data$HIGH_GROUP = recode(as.factor(unmatched_data$HIGH_GROUP), '1'='hDMT', '0'='lDMT')
    names(unmatched_data) <- names(matched_data)
    combined_data = rbind(data.frame(matched_data), data.frame(unmatched_data))
    combined_data$INCLUDED_FOLLOWUP = combined_data$INCLUDED_FOLLOWUP/365.25
    combined_data$DISEASE_DURATION = as.numeric(combined_data$DISEASE_DURATION)/365.25
    combined_data$COUNTRY = factor(combined_data$COUNTRY)
    combined_data$BRAIN_MRI_T2_LESIONCAT_LAB = recode(combined_data$BRAIN_MRI_T2_LESIONCAT ,'0'='0', '1'='1-2', '2'='3-8', '3'='9+')
    table_variables <- c('HIGH_GROUP', 'AGE_BASELINE', 'GENDER', 'EDSS', "BRAIN_MRI_T2_LESIONCAT_LAB", "SC_MRI_LESIONCOUNT", "RELAPSES_1Y", "RELAPSES_2Y", "DISEASE_DURATION", "INCLUDED_FOLLOWUP", "COUNTRY", "CURRENT_CENTER_CODE") # FIRST_TREATMENT
    smd_cat_variables <- c('GENDER', "BRAIN_MRI_T2_LESIONCAT_LAB", "COUNTRY")
    combined_data <- combined_data %>% select(all_of(table_variables))

    var_label(combined_data) <- list(AGE_BASELINE = "Age", EDSS = "EDSS", GENDER = "Gender", RELAPSES_2Y = "Relapses in past 2 years", RELAPSES_1Y = "Relapses in past year", SC_MRI_LESIONCOUNT = "Baseline spinal cord lesions", BRAIN_MRI_T2_LESIONCAT_LAB = "Baseline brain lesions", DISEASE_DURATION = "Years since first symptoms", COUNTRY = "Country", INCLUDED_FOLLOWUP = "Follow-up in years", CURRENT_CENTER_CODE = "CENTER") #FIRST_TREATMENT = "Initiated treatment"

    table1 <- combined_data %>% mutate(HIGH_GROUP = forcats::fct_rev(HIGH_GROUP)) %>% tbl_summary(by=HIGH_GROUP,
                                                                                                  statistic = list(RELAPSES_1Y ~ "{median} ({p25}, {p75})", RELAPSES_2Y ~ "{median} ({p25}, {p75})", EDSS ~ "{median} ({p25}, {p75})", SC_MRI_LESIONCOUNT ~ "{median} ({p25}, {p75})", DISEASE_DURATION ~ "{mean} ({sd})"),
                                                                                                  type = list(RELAPSES_1Y ~ 'continuous', RELAPSES_2Y ~ 'continuous', SC_MRI_LESIONCOUNT ~ 'continuous', EDSS ~ 'continuous', GENDER ~ 'categorical'), # FIRST_TREATMENT ~ 'categorical'
                                                                                                  digits = list(EDSS ~ 1, SC_MRI_LESIONCOUNT ~ 1, RELAPSES_1Y ~ 1, RELAPSES_2Y ~ 1, INCLUDED_FOLLOWUP ~ 1, DISEASE_DURATION ~ 1),
                                                                                                  missing = "ifany",
                                                                                                  missing_text = "(Missing)"
                                                                                                  )  %>% add_stat(fns = everything() ~ unmatched_smd, location = all_of(smd_cat_variables) ~ 'level') %>%
   add_stat(fns = everything() ~ matched_smd, location = all_of(smd_cat_variables) ~ 'level') %>%
    modify_header(
                  list(
                       add_stat_1 ~ "**SMD**",
                       add_stat_2 ~ "**SMD**",
                       stat_3 ~ "**hDMT**, N = {n}",
                       stat_4 ~ "**lDMT**, N = {n}"
                  )
                  ) %>%
    modify_footnote(c(add_stat_1, add_stat_2) ~ "Standardized mean difference") %>%
    modify_table_body(~.x %>% dplyr::relocate(add_stat_1, .after = stat_2)) %>%
    modify_spanning_header(c(stat_1, stat_2) ~ "Unmatched", c(stat_3, stat_4) ~ "Matched") %>%
    modify_fmt_fun(update = list(add_stat_1 ~ function(x) style_number(abs(x), digits=2), add_stat_2 ~ function(x) style_number(abs(x), digits=2)))

    return(table1)
}

# Function to refactor data in dictionary with treatment name as key and years used-up as value
get_treatment_dict <- function(treatment_data, follow_up_time){

  # In case of treatment dict calculations for intervals
  follow_up_time_vec = c()
  if(length(follow_up_time) > 1){
    follow_up_time_vec = follow_up_time
    follow_up_time = sum(follow_up_time)
  }

  group_treatments = c()
  for(i in 1:length(treatment_data)){
    treatments <- fromJSON(gsub("\'","\"", treatment_data[[i]]))

    if(length(follow_up_time_vec) > 1){
      max_time = follow_up_time_vec[[i]]
    }else{
      max_time = follow_up_time
    }

    if(length(names(treatments)) == 0){ next }
      for(j in 1:length(names(treatments))){
        key = names(treatments)[[j]]
        if(key %in% names(group_treatments)){
          group_treatments[key] = group_treatments[key] + min(treatments[key][[1]], max_time)
          #group_treatments[key] = group_treatments[key] + treatments[key][[1]]
        }else{
          group_treatments[key] = min(treatments[key][[1]], max_time)
          #group_treatments[key] = treatments[key][[1]]
        }
      }
  }
  group_treatments["None"] = follow_up_time - sum(group_treatments)
  return(group_treatments)
}

# Function to cluster therapies which make up less than 2% of total into 'Other' category
cluster_small <- function(treatment_dict){
  treatment_dict["Other*"] = 0
  other_dict = c()
  remove = c()
  total_removed = 0
  for(i in 1:length(names(treatment_dict))){
    key = names(treatment_dict)[[i]]
    if(treatment_dict[key] < (sum(treatment_dict) / 50)){
      if(key != "Other*"){
        treatment_dict["Other*"] = treatment_dict["Other*"] + treatment_dict[key]
        other_dict[key] = treatment_dict[key]
        total_removed = total_removed + 1
      }
      remove = append(remove, key)
    }
  }
  if(total_removed > 1){
    for(j in 1:length(remove)){
      if(remove[[j]] != "Other*"){
        treatment_dict = treatment_dict[names(treatment_dict) != remove[[j]]]
      }
    }
      if(treatment_dict["Other*"] == 0){
        treatment_dict = treatment_dict[names(treatment_dict) != "Other*"]
      }
  }else{
    treatment_dict = treatment_dict[names(treatment_dict) != "Other*"]
    other_dict = c()
  }
  return(list(treatment_dict, other_dict))
}

# Function to convert treatment dicts to strings
string_format_treatment_list <- function(treatment_list, total_treatment_time, glossary){
  treatment_string = ""
  treatment_labels = c("AZA" = "azathioprine", "ALE" = "alemtuzumab", "CLA" = "cladribine", "DMF" = "dimethylfumarate", "NAT" = "natalizumab", "FIN" = "fingolimod", "GLA" = "glatiramer acetate", "IVIG" = "immunoglobulins", "IFN" = "interferons", "MIT" = "mitoxantron", "MTX" = "methotrexate", "OFA" = "ofatumumab", "OCR" = "ocrelizumab", "RTX" = "rituximab", "PON" = "ponesimod", "ONA" = "ozanimod", "SIP" = "siponimod", "TER" = "teriflunomide")
  if(glossary){
    return(paste(names(treatment_labels), treatment_labels, sep = " = ", collapse = ", "))
  }
  for(i in 1:length(names(treatment_list))){
    key = names(treatment_list)[[i]]
    if(i == 1){
      first = ""
    }else{
      first = ", "
    }
    treatment_string = paste(treatment_string, first, treatment_labels[key], " (", round(treatment_list[key]/total_treatment_time*100, digits=2), "%)", sep="")
  }
  return(treatment_string)
}
