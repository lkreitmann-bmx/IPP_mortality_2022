#--------------------------------------------------------------------------------------#
# R script used in 2022 paper on IPP gene set to predict 30-day mortality in sepsis 
# patients using publicly available microarray data 
# Script author - LK
# 08/06/2022

#--------------------------------------------------------------------------------------#

rm(list = ls())

#--------------------------------------------------------------------------------------#
####                                  LIBRARIES                                     ####
#--------------------------------------------------------------------------------------#

library(pROC)
library(purrr)
library(ggsci)
library(tidyverse)
library(devtools)
library(knitr)
library(COCONUT)    
library(R.utils)
library(caret)
library(glmnet)
library(foreach)
library(yardstick)
library(tableone)
library(cowplot)

#--------------------------------------------------------------------------------------#
####                               DATA DEFINITION                                  ####
#--------------------------------------------------------------------------------------#

GENEXP_LIST = readRDS("GENEXP_LIST.RDS")
PHENO_LIST = readRDS("PHENO_LIST.RDS")
IPP_v6 = c("GNLY", "CTLA4", "ZAP70", "OAS2", "CD177", "CIITA", "TAP2", "BPGM", "C3AR1", 
           "TDRD9", "CD74", "MDC1", "TNF", "ADGRE3", "CCNB1IP1", "CD3D", "IL10", "CX3CR1", 
           "CD274", "IL7R", "IL1RN", "CXCL10", "IFNG", "IL1R2", "S100A9", "ARL14EP", 
           "norm_1", "norm_2", "norm_3")
studies = GENEXP_LIST %>% imap(~.y) %>% unlist %>% unname
modsLabels = c("Ridge", "Lasso", "ElasticNet", "PLS", "RF", "Radial_SVM")

#--------------------------------------------------------------------------------------#
####                                   FUNCTIONS                                    ####
#--------------------------------------------------------------------------------------#

# this function is used to build COCONUT-normalized data set
build_COCONUT_input_tp_2 = function(study, GENEXP_LIST, PHENO_LIST, tp = "all"){
  
  tmp = list()
  
  if(all(tp == "all")){
    tmp$pheno = PHENO_LIST[[study]] %>% 
      column_to_rownames("subject_ID") %>% 
      as.data.frame
    
  } else {
    tmp$pheno = PHENO_LIST[[study]] %>% 
      filter(timepoint %in% tp | non_infected == 1) %>% 
      column_to_rownames("subject_ID") %>% 
      as.data.frame
  }
  
  tmp$genes = GENEXP_LIST[[study]] %>% 
    as.data.frame() %>%
    select(all_of(rownames(tmp$pheno))) %>% 
    as.matrix
  if(all(colnames(tmp[["genes"]]) == rownames(tmp[["pheno"]]))){
    cat(paste("OK !! in", study, "rows and columns match \n"))
  } else {
    cat(paste("in", study, "rows and columns don't match \n"))
  }
  return(tmp)
}
#--------------------------------------------------------------------------------------#

# this function gives a short table to compare clinical variables in train and test
f_eval_train_test = function(df, seed, by_patient = FALSE, prop){
  set.seed(seed)
  
  # define train and test
  ID_unique = df %>% select(ID_rep, y_died) %>% distinct() %>% pull(ID_rep)              
  tmp = rep(FALSE, length(ID_unique))
  tmp[createDataPartition(y = (df %>% select(ID_rep, y_died) %>% distinct() %>% .$y_died), 
                          p = prop, list = FALSE, times = 1)] = TRUE
  tmp = tibble(ID_unique, tmp) %>% dplyr::rename(ID_rep = ID_unique, TRAIN = tmp)
  trainIndex = left_join(x = df,
                         y = tmp,
                         by = "ID_rep")
  df_cases_train = trainIndex %>% filter(TRAIN == TRUE)
  df_cases_test = trainIndex %>% filter(TRAIN == FALSE)
  
  # build table 1
  tmp = PHENO_LIST %>% map_dfr(~.x) %>% 
    filter(subject_ID %in% df$ID) %>% 
    mutate(TRAIN = trainIndex$TRAIN)
  tmp = left_join(tmp, (trainIndex %>% select(ID, y_died) %>% rename(subject_ID = ID)), 
                  by = "subject_ID") %>% 
    rename(survival = y_died)
  platform = sub("\\..*", "", tmp$ID_rep) %>% 
    map(~ df_cohorts %>% filter(name == .x) %>% pull(platform)) %>% 
    unlist %>% unname
  tmp = tmp %>% mutate(platform = platform)
  
  # by patient or not
  if(by_patient){
    tmp = tmp %>% select(timing, ID_rep, age, bact, ethnic, gender, TRAIN, survival, platform) %>% 
      distinct()
  }
  
  # print table 1
  tableone = CreateTableOne(vars = c("age", "gender", "timing", "bact", "ethnic", 
                                     "platform", "survival"),
                            strata = "TRAIN",
                            data = tmp,
                            factorVars = c("timing", "gender", "platform", "bact", 
                                           "ethnic", "survival"),
                            argsNonNormal = list("timing", "age", "gender", "platform", 
                                                 "bact", "ethnic", "survival"), 
                            includeNA = T)
  print(tableone)
}
#--------------------------------------------------------------------------------------#

## plotting model performance
plot_Models_perf_2 = function(Res_models, data_test = df_cases_test_2, title = NULL){
  
  modsLabels = Res_models %>% names()
  
  # train
  train_AUPRC = Res_models %>% map(~ .x$resample) %>% map_dfc(~ .x$AUPRC) %>% 
    pivot_longer(cols = everything(),
                 names_to = "model_name") %>% 
    mutate(metric = rep("AUPRC", nrow(.)), dataset = rep("train", nrow(.))) %>% 
    arrange(model_name) %>% 
    select(value, model_name, metric, dataset)
  
  train_AUROC = Res_models %>% map(~ .x$resample) %>% map_dfc(~ .x$AUROC) %>% 
    pivot_longer(cols = everything(),
                 names_to = "model_name") %>% 
    mutate(metric = rep("AUROC", nrow(.)), dataset = rep("train", nrow(.))) %>% 
    arrange(model_name) %>% 
    select(value, model_name, metric, dataset)
  
  train_all = bind_rows(train_AUPRC, train_AUROC) %>% 
    mutate(metric = factor(metric, levels = c("AUROC", "AUPRC"))) %>% 
    mutate(model_name = factor(model_name, levels = modsLabels))
  
  h_prev = tibble(metric = "AUPRC",
                  prev = data_test$y_died %>% table %>% prop.table() %>% .[1] %>% unname)
  
  test_point = Res_models %>% map(~ predict(.x, newdata = data_test, type = "prob")) %>% 
    map(~ mutate(.x, y_died = data_test$y_died)) %>% 
    map_dfr(function(.x){
      tibble(AUROC = yardstick::roc_auc(data = .x, truth = y_died, dead) %>% .$.estimate,
             AUPRC = yardstick::pr_auc(data = .x, truth = y_died, dead) %>% .$.estimate)
    }) %>% 
    mutate(model_name = modsLabels) %>% 
    pivot_longer(cols = -model_name,
                 values_to = "value",
                 names_to = "metric") %>% 
    mutate(metric = factor(metric, levels = c("AUROC", "AUPRC"))) %>% 
    mutate(dataset = "train")
  
  plot = train_all %>% 
    # mutate(model_name = factor(model_name, levels = modsLabels),
    #        dataset = factor(dataset, levels = c("train", "test")),
    #        metric = factor(metric, levels = c("AUROC", "AUPRC"))) %>% 
    ggplot(aes(x = model_name, y = value)) + 
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    scale_fill_jama(name = "Models") +
    geom_boxplot(aes(fill = model_name)) +
    ylim(c(0,1)) +
    geom_hline(data = h_prev, aes(yintercept = prev), colour = "grey60", size = 2, linetype = "longdash") +
    geom_point(data = test_point, aes(x = model_name, y = value), 
               color = "black", fill = "grey60", size = 5, shape = 23)  +
    facet_grid(.~ fct_inorder(metric)) +
    labs(title = title, size = 10, fontface = 'bold') +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(size=11),
          strip.text.x = element_text(size = 11))
  
  return(plot)
}
#--------------------------------------------------------------------------------------#

# plot  global gene expression across all cohorts
COCO_boxplot = function(COCO_studies){
  COCOin  = imap(GENEXP_LIST, ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST))
  COCOout = COCONUT(COCOin, control.0.col = "sepsis")
  study = COCOin %>% imap(~.y) %>% unlist %>% unname
  
  df_plot = lapply(COCO_studies, function(df){
    df_raw_dis = COCOout$rawDiseaseList[[df]] %>% .$genes
    df_norm_dis = COCOout$COCONUTList[[df]] %>% .$genes 
    df_norm_cont = COCOout$controlList$GSEs[[df]]  %>% .$genes 
    df_raw_cont = COCOin[[df]] %>% .$genes %>% .[, colnames(.) %in% colnames(df_norm_cont)]
    list_df = list(df_raw_dis = df_raw_dis, 
                   df_norm_dis = df_norm_dis, 
                   df_raw_cont = df_raw_cont, 
                   df_norm_cont = df_norm_cont) %>% 
      map(~ as_tibble(., rownames = "SYMBOL")) %>% 
      map(~ pivot_longer(.x, cols = !SYMBOL, 
                         values_to = "genexp_value",
                         names_to = "subject_ID")) %>% 
      imap(~ mutate(.x, norm_status = factor(ifelse(grepl("raw", .y), rep("raw", nrow(.x)), 
                                                    rep("norm", nrow(.x))),
                                             levels = c("raw", "norm"),
                                             labels = c("Not normalized", "After normalization")))) %>% 
      map(~ .x %>% mutate(study_key = as.character(df))) %>% 
      map_dfr(~ .x) 
    return(list_df)
  }) %>%  map_dfr(~.x) %>% mutate(study_key = fct_inorder(study_key)) %>% 
    group_by(study_key) %>% filter(subject_ID %in% sample(unique(subject_ID), 
                                                          size = round(length(unique(subject_ID))/10, 10), 
                                                          replace = F)) %>% 
    ungroup()
  
  study_key = as.character(df_plot$study_key)
  study_key[grepl("Pankla", study_key)] = "Pankla"
  study_key[study_key == "Bermejo_Martin"] = "Bermejo-Martin"
  study_key[study_key == "Burnham_tot"] = "Burnham"
  study_key[study_key == "Smith_tot"] = "Smith"
  study_key = str_replace(study_key, "_", " ")
  df_plot$study_key = fct_inorder(study_key)
  
  boxplot = df_plot %>% 
    ggplot(aes(x = study_key, y = genexp_value, fill = study_key)) +
    facet_grid(as.factor(norm_status) ~ ., scales = "free") +
    # facet_grid(as.factor(norm_status) ~ .) +
    geom_violin(scale = "area") +
    theme(axis.title.y = element_blank()) +
    scale_fill_viridis(name = "Studies", discrete = T, option = "turbo") +
    xlab("Studies") +
    ylab("Gene expression values") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  return(boxplot)
}
#--------------------------------------------------------------------------------------#

# plot density plots by gene before and after co-normalization (sup_fig_2)
COCO_densityplot = function(COCO_studies, genes_select){
  COCOin  = imap(GENEXP_LIST, ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST))
  COCOout = COCONUT(COCOin, control.0.col = "sepsis")
  study = COCOin %>% imap(~.y) %>% unlist %>% unname 
  
  df_plot = lapply(study, function(df){
    df_raw_dis = COCOout$rawDiseaseList[[df]] %>% .$genes
    df_norm_dis = COCOout$COCONUTList[[df]] %>% .$genes 
    df_norm_cont = COCOout$controlList$GSEs[[df]]  %>% .$genes 
    df_raw_cont = COCOin[[df]] %>% .$genes %>% .[, colnames(.) %in% colnames(df_norm_cont)]
    list_df = list(df_raw_dis = df_raw_dis, 
                   df_norm_dis = df_norm_dis, 
                   df_raw_cont = df_raw_cont, 
                   df_norm_cont = df_norm_cont) %>% 
      map(~ as_tibble(., rownames = "SYMBOL")) %>% 
      map(~ filter(., SYMBOL %in% genes_select)) %>% 
      map(~ pivot_longer(.x, cols = !SYMBOL, 
                         values_to = "genexp_value",
                         names_to = "subject_ID")) %>% 
      imap(~ mutate(.x, norm_status = factor(ifelse(grepl("raw", .y), rep("raw", nrow(.x)), 
                                                    rep("norm", nrow(.x))),
                                             levels = c("raw", "norm"),
                                             labels = c("Not normalized", "Normalized")))) %>% 
      imap(~ mutate(.x, dis_status = factor(ifelse(grepl("cont", .y), rep("Controls", nrow(.x)), 
                                                   rep("Cases", nrow(.x))),
                                            levels = c("Controls", "Cases")))) %>% 
      map_dfr(~ .x) %>% 
      mutate(study_key = rep(as.character(df), nrow(.)))
    return(list_df)
  }) %>% map_dfr(~.x) %>% mutate(study_key = factor(study_key, levels = study)) 
  
  study_key = as.character(df_plot$study_key)
  study_key[grepl("Pankla", study_key)] = "Pankla"
  study_key[study_key == "Bermejo_Martin"] = "Bermejo-Martin"
  study_key[study_key == "Burnham_tot"] = "Burnham"
  study_key[study_key == "Smith_tot"] = "Smith"
  study_key = str_replace(study_key, "_", " ")
  df_plot$study_key = fct_inorder(study_key)
  
  df_plot %>% 
    ggplot(aes(x = genexp_value, color = dis_status, fill = norm_status)) +
    theme_bw() +
    facet_wrap(~ as.factor(dis_status), scales = "free") +
    geom_density_ridges(aes(y = fct_inorder(study_key)), alpha = 0.5, panel_scaling = F, size = 0.8) + 
    scale_color_jama(name = "Disease status", 
                     labels = c("Controls", "Cases")) +
    scale_fill_manual(name = "Normalization", 
                      labels = c("Not normalized", "After normalization"),
                      values =  c("white", "grey")) +
    theme(axis.title.y = element_blank()) +
    xlab("Gene expression values")
}
#--------------------------------------------------------------------------------------#

# plot the effect of co-normalization across studies for 2 genes with different function in sepsis
# (sup_fig_3)
COCO_dotplot = function(COCO_studies, genes_select){
  COCOin  = imap(GENEXP_LIST, ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST))
  COCOout = COCONUT(COCOin, control.0.col = "sepsis")
  study = COCOin %>% imap(~.y) %>% unlist %>% unname
  
  df_plot = lapply(study, function(df){
    df_raw_dis = COCOout$rawDiseaseList[[df]] %>% .$genes
    df_norm_dis = COCOout$COCONUTList[[df]] %>% .$genes 
    df_norm_cont = COCOout$controlList$GSEs[[df]]  %>% .$genes 
    df_raw_cont = COCOin[[df]] %>% .$genes %>% .[, colnames(.) %in% colnames(df_norm_cont)]
    list_df = list(df_raw_dis = df_raw_dis, 
                   df_norm_dis = df_norm_dis, 
                   df_raw_cont = df_raw_cont, 
                   df_norm_cont = df_norm_cont) %>% 
      map(~ as_tibble(., rownames = "SYMBOL")) %>% 
      map(~ filter(., SYMBOL %in% genes_select)) %>% 
      map(~ pivot_longer(.x, cols = !SYMBOL, 
                         values_to = "genexp_value",
                         names_to = "subject_ID")) %>% 
      imap(~ mutate(.x, norm_status = factor(ifelse(grepl("raw", .y), rep("raw", nrow(.x)), 
                                                    rep("norm", nrow(.x))),
                                             levels = c("raw", "norm"),
                                             labels = c("Not normalized", "Normalized")))) %>% 
      imap(~ mutate(.x, dis_status = factor(ifelse(grepl("cont", .y), rep("cont", nrow(.x)), 
                                                   rep("dis", nrow(.x))),
                                            levels = c("cont", "dis")))) %>% 
      map_dfr(~ .x) %>% 
      mutate(study_key = rep(as.character(df), nrow(.)))
    return(list_df)
  }) %>% map_dfr(~.x) 
  
  study_key = pull(df_plot, study_key)
  study_key[grepl("Pankla", study_key)] = "Planka"
  study_key[study_key == "Bermejo_Martin"] = "Bermejo-Martin"
  study_key[study_key == "Burnham_tot"] = "Burnham"
  study_key[study_key == "Smith_tot"] = "Smith"
  study_key = str_replace(study_key, "_", " ")
  df_plot$study_key = study_key
  
  df_plot = df_plot %>% 
    mutate(study_key = fct_inorder(study_key)) %>% 
    mutate(SYMBOL = factor(SYMBOL, 
                           levels = (df_plot %>% group_by(SYMBOL) %>% summarise(mean_ge = mean(genexp_value)) %>% 
                                       ungroup %>% arrange(mean_ge) %>% pull(SYMBOL)))) %>% 
    arrange(SYMBOL, dis_status) %>%
    mutate(subject_ID = subject_ID %>% fct_inorder())
  my_colors = wes_palette(name="Moonrise1")[c(2,3)]
  names(my_colors) = levels(df_plot$SYMBOL)
  plot_res = df_plot %>% 
    ggplot(aes(x = subject_ID, y = genexp_value)) +
    theme(panel.background = element_rect(fill = "white")) +
    theme(panel.border = element_rect(linetype = "solid", fill = NA)) +
    geom_point(aes(color = SYMBOL, shape = dis_status)) +
    scale_shape_manual(name = "Disease status", labels = c("Controls", "Cases"),
                       values = c(21, 19)) +
    scale_color_manual(name = "Genes", values=c("#69b3a2", "#404080")) +
    facet_grid(norm_status ~ study_key, scales = "free") +
    ylim(c(0,20)) +
    ylab("Gene Expression Values") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("Subjects") 
  
  plot_res
}
#--------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------#
#####                                   DATA SET 1                                 #####
#--------------------------------------------------------------------------------------#

## run COCONUT
COCOin  = imap(GENEXP_LIST, ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST, tp = 1))
COCOout = COCONUT(COCOin, control.0.col = "sepsis", byPlatform = FALSE)

## build matrix of predictors
df_cases_1 = COCOout$COCONUTList %>% 
  map_dfc(~ .x$genes) %>% # gather all cases after co-normalization
  as.matrix %>% t %>% as_tibble(rownames = "ID")  
tmp = PHENO_LIST %>% map_dfr(~.x) %>% select(subject_ID, ID_rep)
wh = lapply(df_cases_1$ID, function(x){which(tmp$subject_ID == x)}) %>% unlist
tmp = tmp[wh,] %>% mutate(verif_ID = df_cases_1$ID) 
df_cases_1 = df_cases_1 %>% mutate(ID_rep = tmp$ID_rep) %>% select(ID_rep, everything())

# build vector of outcomes
y_died_1 = COCOout$COCONUTList %>% 
  map_dfr(~ .x$pheno) %>% 
  filter(rownames(.) %in% df_cases_1$ID) %>% 
  pull(died)
y_died_1 = ifelse(y_died_1 == 1, "dead", "alive") %>% factor(levels  = c("dead", "alive"))

# build model matrix
df_cases_1 = df_cases_1 %>% mutate(y_died = fct_inorder(y_died_1))

# standardize
tmp = df_cases_1 %>% select(- c(ID_rep, ID, y_died)) %>% lapply(function(vec){
  mean = mean(vec)
  std = sd(vec)
  vec = (vec - mean)/std
})
df_cases_1 = bind_cols(df_cases_1 %>% select(ID_rep, ID),
                       as_tibble(tmp),
                       tibble(y_died = df_cases_1 %>% pull(y_died)))

## train/test split
set.seed(7)
ID_unique = df_cases_1 %>% select(ID_rep, y_died) %>% distinct() %>% pull(ID_rep)              
tmp = rep(FALSE, length(ID_unique))
tmp[caret::createDataPartition(y = (df_cases_1 %>% select(ID_rep, y_died) %>% distinct() %>% .$y_died), 
                               p = .7, list = FALSE, times = 1)] = TRUE
tmp = tibble(ID_unique, tmp) %>% dplyr::rename(ID_rep = ID_unique, TRAIN = tmp)
trainIndex_1 = left_join(x = df_cases_1,
                         y = tmp,
                         by = "ID_rep") %>% 
  mutate(y_died = factor(y_died, levels = c("dead", "alive")))
df_cases_train_1 = trainIndex_1 %>% filter(TRAIN == TRUE)
df_cases_test_1 = trainIndex_1 %>% filter(TRAIN == FALSE)

## compare phenotypes (age, sex, ethnicity, chip, com/noso)
f_eval_train_test(df = df_cases_1, seed = 7, by_patient = F, prop = 0.7)

#--------------------------------------------------------------------------------------#
#####                                  DATA SET 2                                  #####
#--------------------------------------------------------------------------------------#

# remove studies that don't have data on the right tp
tp = c(3:7)
studies_to_remove_tp = studies[!(PHENO_LIST %>% 
                                   map(~ .x$timepoint %in% tp) %>% 
                                   map(~ sum(.x)) %>% 
                                   unlist %>% 
                                   as.logical)]
studies_in = setdiff(studies, studies_to_remove_tp)

## run COCONUT
COCOin  = imap(GENEXP_LIST[studies_in], ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST, tp = tp))
COCOout = COCONUT(COCOin, control.0.col = "sepsis", byPlatform = FALSE)

## build matrix of predictors
df_cases_2 = COCOout$COCONUTList %>%
  map_dfc(~ .x$genes) %>% # gather all cases after co-normalization
  as.matrix %>% t %>% as_tibble(rownames = "ID")
tmp = PHENO_LIST %>% map_dfr(~.x) %>% select(subject_ID, ID_rep)
wh = lapply(df_cases_2$ID, function(x){which(tmp$subject_ID == x)}) %>% unlist
tmp = tmp[wh,] %>% mutate(verif_ID = df_cases_2$ID)
df_cases_2 = df_cases_2 %>% mutate(ID_rep = tmp$ID_rep) %>% select(ID_rep, everything())

# build vector of outcomes
y_died_2 = COCOout$COCONUTList %>%
  map_dfr(~ .x$pheno) %>%
  filter(rownames(.) %in% df_cases_2$ID) %>%
  pull(died)
y_died_2 = ifelse(y_died_2 == 1, "dead", "alive") %>% factor(levels  = c("dead", "alive"))

# build model matrix
df_cases_2 = df_cases_2 %>% mutate(y_died = fct_inorder(y_died_2))

# standardize
tmp = df_cases_2 %>% select(- c(ID_rep, ID, y_died)) %>% lapply(function(vec){
  mean = mean(vec)
  std = sd(vec)
  vec = (vec - mean)/std
})
df_cases_2 = bind_cols(df_cases_2 %>% select(ID_rep, ID),
                       as_tibble(tmp),
                       tibble(y_died = df_cases_2 %>% pull(y_died)))

## train/test split
set.seed(7)
ID_unique = df_cases_2 %>% select(ID_rep, y_died) %>% distinct() %>% pull(ID_rep)
tmp = rep(FALSE, length(ID_unique))
tmp[caret::createDataPartition(y = (df_cases_2 %>% select(ID_rep, y_died) %>% distinct() %>% .$y_died), 
                               p = .7, list = FALSE, times = 1)] = TRUE
tmp = tibble(ID_unique, tmp) %>% dplyr::rename(ID_rep = ID_unique, TRAIN = tmp)
trainIndex_2 = left_join(x = df_cases_2,
                         y = tmp,
                         by = "ID_rep") %>%
  mutate(y_died = factor(y_died, levels = c("dead", "alive")))
df_cases_train_2 = trainIndex_2 %>% filter(TRAIN == TRUE)
df_cases_test_2 = trainIndex_2 %>% filter(TRAIN == FALSE)
df_cases_train_2$y_died %>% table %>% prop.table()

#------------------------------------------------------------------------------------------#
######                                     INITIALIZE TUNING                                  ######
#------------------------------------------------------------------------------------------#

## essential mySummary function
my_Summary = function(data, lev = NULL, model = NULL){
  c(
    AUROC = yardstick::roc_auc(data = data, truth = obs, dead) %>% .$.estimate,
    AUPRC = yardstick::pr_auc(data = data, truth = obs, dead) %>% .$.estimate,
    Se =  yardstick::sens(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    Sp =  yardstick::spec(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    PPV =  yardstick::spec(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    NPV =  yardstick::spec(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    Kappa = yardstick::kap(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    mcc = yardstick::mcc(data = data, truth = obs, estimate = pred) %>% .$.estimate,
    Accuracy = yardstick::accuracy(data = data, truth = obs, estimate = pred) %>% .$.estimate
  )
}

## initialize tuning
metric = "AUPRC"
set.seed(123)
modsLabels = c("Ridge", "Lasso", "Elastic_Net", "PLS", "RF", "Radial_SVM") 
mod = NULL
k_fold = 10
r_repet = 5
preProc = c("center","scale")
samplingOpt = NULL
# samplingOpt = "smote"

## train_Control
train_Control <- trainControl(
  method = "repeatedcv",
  number = k_fold,
  repeats = r_repet,
  sampling = samplingOpt, 
  selectionFunction = "best",
  classProbs = TRUE, 
  summaryFunction = my_Summary,
  verboseIter = TRUE,
  savePredictions = TRUE,
  returnResamp = "final",
  allowParallel = FALSE
)

#------------------------------------------------------------------------------------------#
######                                  MODEL IPPv6 tp = 1                            ######
#------------------------------------------------------------------------------------------#

MOD_IPP_d1 = NULL
genes = IPP_v6
data = df_cases_train_1 %>% select(all_of(genes), y_died)

## Ridge
tune_Grid = expand.grid(.alpha = 0, .lambda = seq(0.05, 2, length = 5))
MOD_IPP_d1$Ridge = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet", 
                         trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_IPP_d1$Ridge)

## Lasso
tune_Grid = expand.grid(alpha = 1, lambda = seq(0.000001, 0.05, length = 40))
MOD_IPP_d1$Lasso = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_IPP_d1$Lasso)

## EN
tune_Grid = expand.grid(alpha = seq(0.1, 0.9, length = 10),
                        lambda = seq(0.00001, 0.005, length = 10))
MOD_IPP_d1$ElasticNet = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet", 
                              trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_IPP_d1$ElasticNet)

##  PLS
MOD_IPP_d1$PLS = caret::train(y_died ~ ., data = data, method = 'pls', metric = 'AUROC',
                       tuneLength  = 6, trControl = train_Control, preProcess = preProc)
plot(MOD_IPP_d1$PLS)

## RF
mtry = sqrt(ncol(data))
tune_Grid = expand.grid(.mtry = 2:sqrt(ncol(data)))
MOD_IPP_d1$RF = caret::train(y_died ~ ., data = data, method = 'rf', metric = 'AUROC',
                      tuneGrid = tune_Grid, trControl = train_Control, preProcess = preProc)
plot(MOD_IPP_d1$RF)

## SVM 
tune_Grid = expand.grid(C = seq(0.001, 10, length = 4),
                        sigma = 10^seq(-3.5, -2, length = 4))
MOD_IPP_d1$Radial_SVM = caret::train(y_died ~., data = data, method = "svmRadial", metric = "AUROC",
                              trControl = train_Control, preProcess = preProc, tuneGrid = tune_Grid)
plot(MOD_IPP_d1$Radial_SVM)

#------------------------------------------------------------------------------------------#
######                              MODEL all genes tp = 1                            ######
#------------------------------------------------------------------------------------------#

MOD_all_d1 = NULL
all_genes = COCOout$controlList$bayesParams$gamma.star %>% colnames
genes = all_genes
data = df_cases_train_1 %>% select(all_of(genes), y_died)

## Ridge
tune_Grid = expand.grid(.alpha = 0, .lambda = 10^seq(-5, 1, length = 10))
MOD_all_d1$Ridge = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_all_d1$Ridge)

## Lasso
tune_Grid = expand.grid(alpha = 1, lambda = seq(0.01, 0.05, length = 10))
MOD_all_d1$Lasso = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid)
plot(MOD_all_d1$Lasso)

## EN
tune_Grid = expand.grid(alpha = seq(0.15, 0.9, length = 7),
                        lambda = seq(0.00001, 0.05, length = 6))
MOD_all_d1$ElasticNet = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                              trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_all_d1$ElasticNet)

##  PLS
MOD_all_d1$PLS = caret::train(y_died ~ ., data = data, method = 'pls', metric = 'AUROC',
                       tuneLength  = 10, trControl = train_Control, preProcess = preProc)
plot(MOD_all_d1$PLS)

## RF
mtry = sqrt(ncol(data))
tune_Grid = expand.grid(.mtry = 2:mtry)
MOD_all_d1$RF = caret::train(y_died ~ ., data = data, method = 'rf', metric = 'AUROC',
                      tuneGrid = tune_Grid, trControl = train_Control, preProcess = preProc)
plot(MOD_all_d1$RF)

## SVM 
tune_Grid = expand.grid(C = seq(0.001, 10, length = 4),
                        sigma = 10^seq(-6, -2, length = 4))
MOD_all_d1$Radial_SVM = caret::train(y_died ~., data = data, method = "svmRadial", metric = "AUROC",
                              trControl = train_Control, preProcess = preProc, tuneGrid = tune_Grid)
plot(MOD_all_d1$Radial_SVM)

# alternatively this can be done in Python (much faster)

#------------------------------------------------------------------------------------------#
######                              MODEL top 29 genes tp = 1                         ######
#------------------------------------------------------------------------------------------#

MOD_top29_d1 = NULL
genes = c("SYT7", "KCNK13", "PCDHA5", "DNAJC5G", "ZNF354A", "PRG2", "RGS1", "MT1F", "TRIB3",
          "DIAPH3", "GSTO2", "NT5DC3", "CCR9", "KPNA6", "CHST13", "CUBN", "HINT3",
          "TMEM128", "NPDC1", "SLC39A9", "RABEPK", "TMEM129", "ARHGAP11A", "OAZ3", "EPC2",
          "OPRK1", "CAMTA2", "CABP2", "HBEGF")
data = df_cases_train_1 %>% select(all_of(genes), y_died)

## Ridge
tune_Grid = expand.grid(.alpha = 0, .lambda = seq(0.05, 2, length = 5))
MOD_top29_d1$Ridge = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_top29_d1$Ridge)

## Lasso
tune_Grid = expand.grid(alpha = 1, lambda = seq(0.000001, 0.05, length = 40))
MOD_top29_d1$Lasso = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_top29_d1$Lasso)

## EN
tune_Grid = expand.grid(alpha = seq(0.1, 0.9, length = 10),
                        lambda = seq(0.00001, 0.005, length = 10))
MOD_top29_d1$ElasticNet = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                              trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_top29_d1$ElasticNet)

##  PLS
MOD_top29_d1$PLS = caret::train(y_died ~ ., data = data, method = 'pls', metric = 'AUROC',
                       tuneLength  = 6, trControl = train_Control, preProcess = preProc)
plot(MOD_top29_d1$PLS)

## RF
mtry = sqrt(ncol(data))
tune_Grid = expand.grid(.mtry = 2:sqrt(ncol(data)))
MOD_top29_d1$RF = caret::train(y_died ~ ., data = data, method = 'rf', metric = 'AUROC',
                      tuneGrid = tune_Grid, trControl = train_Control, preProcess = preProc)
plot(MOD_top29_d1$RF)

## SVM 
tune_Grid = expand.grid(C = seq(0.001, 10, length = 4),
                        sigma = 10^seq(-3.5, -2, length = 4))
MOD_top29_d1$Radial_SVM = caret::train(y_died ~., data = data, method = "svmRadial", metric = "AUROC",
                              trControl = train_Control, preProcess = preProc, tuneGrid = tune_Grid)
plot(MOD_top29_d1$Radial_SVM)

#------------------------------------------------------------------------------------------#
######                                MODEL IPP tp = 3:7                              ######
#------------------------------------------------------------------------------------------#

MOD_IPP_d3_7 = NULL
genes = IPP_v6
data = df_cases_train_2 %>% select(all_of(genes), y_died)

## Ridge
tune_Grid = expand.grid(.alpha = 0, .lambda = seq(0,01, 2, length = 10))
MOD_IPP_d3_7$Ridge = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid)
plot(MOD_IPP_d3_7$Ridge)

## Lasso
tune_Grid = expand.grid(alpha = 1, lambda = seq(0.000001, 0.2, length = 10))
MOD_IPP_d3_7$Lasso = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                         trControl = train_Control, tuneGrid = tune_Grid)
plot(MOD_IPP_d3_7$Lasso)

## EN
tune_Grid = expand.grid(alpha = seq(0.1, 0.9, length = 10),
                        lambda = seq(0.00001, 0.3, length = 10))
MOD_IPP_d3_7$ElasticNet = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                              trControl = train_Control, tuneGrid = tune_Grid)
plot(MOD_IPP_d3_7$ElasticNet)

##  PLS
MOD_IPP_d3_7$PLS = caret::train(y_died ~ ., data = data, method = 'pls', metric = 'AUROC',
                       tuneLength  = 6, trControl = train_Control, preProcess = preProc)
plot(MOD_IPP_d3_7$PLS)

## RF
mtry = sqrt(ncol(data))
tune_Grid = expand.grid(.mtry = 2:mtry)
MOD_IPP_d3_7$RF = caret::train(y_died ~ ., data = data, method = 'rf', metric = 'AUROC', tuneGrid = tune_Grid, 
                      trControl = train_Control, preProcess = preProc)
plot(MOD_IPP_d3_7$RF)

## SVM 
tune_Grid = expand.grid(C = seq(0.001, 10, length = 4),
                        sigma = 10^seq(-5, -3, length = 4))
MOD_IPP_d3_7$Radial_SVM = caret::train(y_died ~., data = data, method = "svmRadial", metric = "AUROC",
                              trControl = train_Control, preProcess = preProc, tuneGrid = tune_Grid)
plot(MOD_IPP_d3_7$Radial_SVM)

#------------------------------------------------------------------------------------------#
######                            MODEL all genes tp = 3:7                            ######
#------------------------------------------------------------------------------------------#

MOD_all_d3_7 = NULL
all_genes = COCOout$controlList$bayesParams$gamma.star %>% colnames
genes = all_genes
data = df_cases_train_1 %>% select(all_of(genes), y_died)

## Ridge
tune_Grid = expand.grid(.alpha = 0, .lambda = 10^seq(-5, 1, length = 10))
MOD_all_d3_7$Ridge = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                                trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_all_d3_7$Ridge)

## Lasso
tune_Grid = expand.grid(alpha = 1, lambda = seq(0.01, 0.05, length = 10))
MOD_all_d3_7$Lasso = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                                trControl = train_Control, tuneGrid = tune_Grid)
plot(MOD_all_d3_7$Lasso)

## EN
tune_Grid = expand.grid(alpha = seq(0.15, 0.9, length = 7),
                        lambda = seq(0.00001, 0.05, length = 6))
MOD_all_d3_7$ElasticNet = caret::train(y_died ~., data = data, metric = "AUROC", method = "glmnet",
                                     trControl = train_Control, tuneGrid = tune_Grid, preProcess = preProc)
plot(MOD_all_d3_7$ElasticNet)

##  PLS
MOD_all_d3_7$PLS = caret::train(y_died ~ ., data = data, method = 'pls', metric = 'AUROC',
                              tuneLength  = 10, trControl = train_Control, preProcess = preProc)
plot(MOD_all_d3_7$PLS)

## RF
mtry = sqrt(ncol(data))
tune_Grid = expand.grid(.mtry = 2:mtry)
MOD_all_d3_7$RF = caret::train(y_died ~ ., data = data, method = 'rf', metric = 'AUROC',
                             tuneGrid = tune_Grid, trControl = train_Control, preProcess = preProc)
plot(MOD_all_d3_7$RF)

## SVM 
tune_Grid = expand.grid(C = seq(0.001, 10, length = 4),
                        sigma = 10^seq(-6, -2, length = 4))
MOD_all_d3_7$Radial_SVM = caret::train(y_died ~., data = data, method = "svmRadial", metric = "AUROC",
                                     trControl = train_Control, preProcess = preProc, tuneGrid = tune_Grid)
plot(MOD_all_d3_7$Radial_SVM)

#--------------------------------------------------------------------------------------#
####                                      Table 2                                   ####
#--------------------------------------------------------------------------------------#

tab_2 = f_eval_train_test(df = df_cases_1, seed = 7, by_patient = F, prop = 0.7)

#--------------------------------------------------------------------------------------#
####                                      Figure 1                                  ####
#--------------------------------------------------------------------------------------#

## run COCONUT
COCOin  = imap(GENEXP_LIST[studies], ~ build_COCONUT_input_tp_2(.y, GENEXP_LIST, PHENO_LIST, tp = 1))
COCOout = COCONUT(COCOin, control.0.col = "sepsis", byPlatform = FALSE)

# norm
df_pca_norm = COCOout$COCONUTList %>% 
  map(~ as.matrix(t(.x$genes))) %>%  
  map(~ .x %>% as_tibble(rownames = "ID")) %>% 
  imap_dfr(~ .x %>% mutate(studies = .y))
pca_norm = prcomp(df_pca_norm[,-c(1, ncol(df_pca_norm))], retx = TRUE, center = TRUE, scale. = TRUE)
pc1 = pca_norm$x[,"PC1"]
pc2 = pca_norm$x[,"PC2"]
studies = pull(df_pca_norm, studies)
studies[grepl("Pankla", studies)] = "Planka"
studies[studies == "Bermejo_Martin"] = "Bermejo-Martin"
studies[studies == "Burnham_tot"] = "Burnham"
studies[studies == "Smith_tot"] = "Smith"
studies = str_replace(studies, "_", " ")
tmp = table(fct_inorder(studies))
colors = rep(viridis::turbo(17), times = unname(tmp)) 
shapes = rep(rep(15:20, 3)[1:17], unname(tmp))
plot_pca_norm = tibble(pc1 = pc1, 
                       pc2 = pc2,
                       studies = studies,
                       colors = colors,
                       shapes = shapes)
colors = plot_pca_norm$colors %>% set_names(nm = plot_pca_norm$studies)
shapes = plot_pca_norm$shapes %>% set_names(nm = plot_pca_norm$studies)
var_exp = pca_norm$sdev^2/sum(pca_norm$sdev^2)*100
var_exp = round(var_exp, digits = 1)
p_norm = ggplot(data = plot_pca_norm, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = studies, shape = studies), size = 3, alpha = 1) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"))

# raw
common_genes = colnames(df_pca_norm)[2:(ncol(df_pca_norm) - 1)]
df_pca_raw = COCOout$rawDiseaseList %>%
  map(~ as.matrix(t(.x$genes))) %>% 
  map(~ .x[, common_genes]) %>% 
  map(~ .x %>% as_tibble(rownames = "ID")) %>% 
  imap_dfr(~ .x %>% mutate(group = .y))
pca_raw = prcomp(df_pca_raw[,-c(1, ncol(df_pca_raw))], retx = TRUE, center = TRUE, scale. = TRUE)
pc1 = pca_raw$x[,"PC1"]
pc2 = pca_raw$x[,"PC2"]
tmp = table(fct_inorder(studies))
colors = rep(viridis::turbo(17), times = unname(tmp)) 
shapes = rep(rep(15:20, 3)[1:17], unname(tmp))
plot_pca_raw = tibble(pc1 = pc1, 
                      pc2 = pc2,
                      studies = studies,
                      colors = colors,
                      shapes = shapes)
colors = plot_pca_raw$colors %>% set_names(nm = plot_pca_raw$studies)
shapes = plot_pca_raw$shapes %>% set_names(nm = plot_pca_raw$studies)
p_raw = ggplot(data = plot_pca_raw, aes(x = pc1, y = pc2)) +
  geom_point(aes(color = studies, shape = studies), size = 3) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme(legend.position="none") +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"))

fig_1 = p_raw + p_norm

#--------------------------------------------------------------------------------------#
####                                      Figure 2                                  ####
#--------------------------------------------------------------------------------------#

fig_2 = plot_Models_perf_2(MOD_IPP_d1, data_test = df_cases_test_1)

#--------------------------------------------------------------------------------------#
####                                      Figure 3                                  ####
#--------------------------------------------------------------------------------------#

data_test = df_cases_test_1

pred1 = predict(MOD_IPP_d1$RF, data_test, type = "prob") %>% mutate(y_died = data_test$y_died)
roc_curve1 = pROC::roc(pred1$y_died, pred1$dead, ci = TRUE, plot = FALSE)
ciobj1 <- pROC::ci.se(roc_curve1, specificities=seq(0, 1, l=25))
dat.ci1 <- data.frame(x = as.numeric(rownames(ciobj1)),
                      lower = ciobj1[, 1],
                      upper = ciobj1[, 3]) %>% as_tibble(rownames = "step") %>% 
  mutate(model = "roc_mod1")

pred2 = predict(MOD_all_d1$RF, data_test, type = "prob") %>% mutate(y_died = data_test$y_died)
roc_curve2 = pROC::roc(pred2, response = truth, predictor = pred, ci = TRUE)
ciobj2 <- pROC::ci.se(roc_curve2, specificities=seq(0, 1, l=25))
dat.ci2 <- data.frame(x = as.numeric(rownames(ciobj2)),
                      lower = ciobj2[, 1],
                      upper = ciobj2[, 3]) %>% as_tibble(rownames = "step") %>% 
  mutate(model = "roc_mod2")

pred3 = predict(MOD_top29_d1$RF, data_test, type = "prob") %>% mutate(y_died = data_test$y_died)
roc_curve3 = pROC::roc(pred3$y_died, pred3$dead, ci = TRUE, plot = FALSE)
ciobj3 <- pROC::ci.se(roc_curve3, specificities=seq(0, 1, l=25))
dat.ci3 <- data.frame(x = as.numeric(rownames(ciobj3)),
                      lower = ciobj3[, 1],
                      upper = ciobj3[, 3]) %>% as_tibble(rownames = "step") %>% 
  mutate(model = "roc_mod3")

dat.ci_tot = list(dat.ci1 = dat.ci1, dat.ci2 = dat.ci2, dat.ci3 = dat.ci3)

colors = pal_jama()(6)[c(4, 1, 2)]

test1 = pROC::roc.test(roc_curve1, roc_curve2)
test2 = pROC::roc.test(roc_curve1, roc_curve3)

fig_3 = pROC::ggroc(list("IPP" = roc_curve1, "all genes" = roc_curve2, "top 29 genes" = roc_curve3), size = 1) + 
  theme_minimal() + 
  scale_color_manual(values = colors, name = "Models") +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha = 0.7, color = "grey") + 
  coord_equal() +
  ggtitle(paste0(gsub("Area under the curve", "AUROC", capture.output(roc_curve1$auc)), ", ", 
                 gsub("[(DeLong)]", "", capture.output(roc_curve1$ci)), " for IPP model\n",
                 gsub("Area under the curve", "AUROC", capture.output(roc_curve2$auc)), ", ", 
                 gsub("[(DeLong)]", "", capture.output(roc_curve2$ci)), " for 'all genes' model\n",
                 gsub("Area under the curve", "AUROC", capture.output(roc_curve3$auc)), ", ", 
                 gsub("[(DeLong)]", "", capture.output(roc_curve3$ci)), " for 'top 29 genes' model\n",
                 "IPP vs. 'all genes', p=", round(test1$p.value, 3), " (DeLong)", "\n",
                 "IPP vs. 'top 29 genes', p=", round(test2$p.value, 3), " (DeLong)", "\n")
  ) +
  theme(plot.title = element_text(size = 11),
        legend.text = element_text(size=11),)

for(i in c(1:3)){
  p = p + 
    geom_ribbon(data = dat.ci_tot[[i]], 
                aes(x = x, ymin = lower, ymax = upper),
                fill = colors[i],
                alpha = 0.15, 
                inherit.aes = F) 
}

#--------------------------------------------------------------------------------------#
####                                      Figure 4                                  ####
#--------------------------------------------------------------------------------------#

data_test = df_cases_test_3

pred1 = predict(MOD_IPP_d3_7$RF, data_test, type = "prob") %>% mutate(y_died = data_test$y_died)
roc_curve1 = pROC::roc(pred1$y_died, pred1$dead, ci = TRUE, plot = FALSE)
ciobj1 <- pROC::ci.se(roc_curve1, specificities=seq(0, 1, l=25))
dat.ci1 <- data.frame(x = as.numeric(rownames(ciobj1)),
                      lower = ciobj1[, 1],
                      upper = ciobj1[, 3]) %>% as_tibble(rownames = "step") %>% 
  mutate(model = "roc_mod1")

pred2 = predict(MOD_all_d3_7$Ridge, data_test, type = "prob") %>% mutate(y_died = data_test$y_died)
roc_curve2 = pROC::roc(pred2$y_died, pred2$dead, ci = TRUE, plot = FALSE)
ciobj2 <- pROC::ci.se(roc_curve2, specificities=seq(0, 1, l=25))
dat.ci2 <- data.frame(x = as.numeric(rownames(ciobj2)),
                      lower = ciobj2[, 1],
                      upper = ciobj2[, 3]) %>% as_tibble(rownames = "step") %>% 
  mutate(model = "roc_mod2")

dat.ci_tot = list(dat.ci1 = dat.ci1, dat.ci2 = dat.ci2)#, dat.ci3 = dat.ci3)

colors = pal_jama()(6)[c(4, 1)]

test1 = pROC::roc.test(roc_curve1, roc_curve2)

fig_4 = pROC::ggroc(list("IPP" = roc_curve1, "all genes" = roc_curve2), size = 1) + 
  theme_minimal() + 
  scale_color_manual(values = colors, name = "Models") +
  geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha = 0.7, color = "grey") + 
  coord_equal() +
  ggtitle(paste0(gsub("Area under the curve", "AUROC", capture.output(roc_curve1$auc)), ", ", 
                 gsub("[(DeLong)]", "", capture.output(roc_curve1$ci)), " for IPP model\n",
                 gsub("Area under the curve", "AUROC", capture.output(roc_curve2$auc)), ", ", 
                 gsub("[(DeLong)]", "", capture.output(roc_curve2$ci)), " for 'all genes' model\n",
                 "IPP vs. 'all genes', p=", round(test1$p.value, 3), " (DeLong)")
  ) +
  theme(plot.title = element_text(size = 11),
        legend.text = element_text(size=11),)

for(i in c(1:2)){
  p = p + 
    geom_ribbon(data = dat.ci_tot[[i]], 
                aes(x = x, ymin = lower, ymax = upper),
                fill = colors[i],
                alpha = 0.15, 
                inherit.aes = F) 
}

#--------------------------------------------------------------------------------------#
####                                      Figure 5                                  ####
#--------------------------------------------------------------------------------------#

f_threshold = function(mod, df){
  roc_result = f_ROC(mod, df)
  threshold = coords(roc_result, 
                     ret = c('th', "sen", "sp", "ppv", "npv", "precision", "recall")) %>% 
    filter(sensitivity > 0.8) %>% 
    filter(specificity == max(specificity)) %>% 
    filter(sensitivity == max(sensitivity)) %>% 
    pull(threshold)
  return(threshold)
}

# DAY 1
df_pred = df_cases_test_1

table_enrichment_d1 = MOD2_final$RF %>% predict(newdata = df_pred, type = "prob") %>% 
  select(dead) %>% 
  mutate(y_died = df_pred$y_died) %>% 
  pROC::roc(data = .,
            response = y_died,
            predictor = dead,
            levels = c("alive", "dead"),
            direction = "<",
            ci = TRUE, 
            plot = FALSE)

threshold = coords(table_enrichment_d1, 
                   ret = c('th', "sen", "sp", "ppv", "npv", "precision", "recall")) %>% 
  filter(sensitivity > 0.8) %>% 
  filter(specificity == max(specificity)) %>% 
  filter(sensitivity == max(sensitivity)) %>% 
  pull(threshold)

threshold2 = coords(table_enrichment_d1, x = "best", best.method = "closest.topleft") %>% 
  pull(threshold)

table_enrichment_d1 = MOD2_final$RF %>% predict(newdata = df_pred, type = "prob") %>% 
  select(dead) %>% 
  mutate(y_died = df_pred$y_died) %>% 
  mutate(predict = ifelse(dead > threshold2, "High risk", "Low risk")) %>% 
  group_by(predict) %>% 
  summarise(n = n(), 
            n_dead = sum(ifelse(y_died == "alive", 0, 1))) 
prop_dead = table_enrichment_d1$n_dead/table_enrichment_d1$n
CI_half = NULL
for(i in 1:nrow(table_enrichment_d1)){
  CI_half[i] = binom.test(x = table_enrichment_d1$n_dead[i], 
                          n = table_enrichment_d1$n[i]) %>% 
    .$conf.int %>% as.numeric() %>% .[2] - prop_dead[i]
}

table_enrichment_d1 = table_enrichment_d1 %>% mutate(prop_dead = prop_dead,
                                                     CI_half = CI_half) %>% 
  mutate(timepoint = "day 1")

# DAY >2
df_pred = df_cases_test_3

table_enrichment_d3 = MOD10_def$RF %>% predict(newdata = df_pred, type = "prob") %>% 
  select(dead) %>% 
  mutate(y_died = df_pred$y_died) %>% 
  pROC::roc(data = .,
            response = y_died,
            predictor = dead,
            levels = c("alive", "dead"),
            direction = "<",
            ci = TRUE, 
            plot = FALSE)

threshold2 = coords(table_enrichment_d3, x = "best", best.method = "closest.topleft") %>% 
  pull(threshold)

table_enrichment_d3 = MOD10_def$RF %>% predict(newdata = df_pred, type = "prob") %>% 
  select(dead) %>% 
  mutate(y_died = df_pred$y_died) %>% 
  mutate(predict = ifelse(dead > threshold2, "High risk", "Low risk")) %>% 
  group_by(predict) %>% 
  summarise(n = n(), 
            n_dead = sum(ifelse(y_died == "alive", 0, 1))) 
prop_dead = table_enrichment_d3$n_dead/table_enrichment_d3$n
CI_half = NULL
for(i in 1:nrow(table_enrichment_d3)){
  CI_half[i] = binom.test(x = table_enrichment_d3$n_dead[i], 
                          n = table_enrichment_d3$n[i]) %>% 
    .$conf.int %>% as.numeric() %>% .[2] - prop_dead[i]
}

table_enrichment_d3 = table_enrichment_d3 %>% mutate(prop_dead = prop_dead,
                                                     CI_half = CI_half) %>% 
  mutate(timepoint = "day >2")

table_enrichment_tot = bind_rows(table_enrichment_d1, table_enrichment_d3) %>% 
  mutate(timepoint = fct_inorder(timepoint)) %>% 
  mutate(predict = factor(predict, levels = c("Low risk", "High risk"))) %>% 
  mutate(prop_dead = 100*prop_dead) %>% 
  mutate(CI_half = 100*CI_half)

my_plot = table_enrichment_tot %>% 
  arrange(timepoint, predict) %>% 
  mutate(
    newXtick = as.character(1:n())
  ) %>% 
  mutate(annotations = c(rep(p1, 2), rep(p2, 2)))

p1 = chisq.test(table_enrichment_tot %>%
                  filter(timepoint == "day 1") %>%
                  select(n, n_dead) %>%
                  mutate(n_alive = n - n_dead) %>%
                  select(n_dead, n_alive)) %>%
  .$p.value

p2 = fisher.test(table_enrichment_tot %>%
                   filter(timepoint == "day >2") %>%
                   select(n, n_dead) %>%
                   mutate(n_alive = n - n_dead) %>%
                   select(n_dead, n_alive)) %>%
  .$p.value

my_plot_1 = my_plot %>% 
  filter(timepoint == "day 1") %>% 
  ggplot(aes(x = newXtick, y = prop_dead, fill = predict)) +
  geom_bar(position = dodge, stat = "identity", color = "black") +
  scale_fill_manual(values = c("#3c9845", "#00427f")) +
  ylim(c(0,100)) +
  geom_errorbar(aes(ymin = prop_dead, 
                    ymax = prop_dead + CI_half),
                color = "black",
                alpha = 1,
                width =.2,
                position = dodge) +
  ylab("Proportion of dead patients") +
  facet_wrap(~ "Day 1", scales = "free_x") +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = "None",
        
        strip.text.x = element_text(size = 11, face = "bold")) +
  xlab("") +
  ggsignif::geom_signif(annotations = paste0("p=", signif(p1, digits = 3)), 
                        comparisons = list(c("1", "2")),
                        y_position = 40)

my_plot_2  = my_plot %>% 
  filter(timepoint == "day >2") %>% 
  ggplot(aes(x = newXtick, y = prop_dead, fill = predict)) +
  geom_bar(position = dodge, stat = "identity", color = "black") +
  scale_fill_manual(values = c("#3c9845", "#00427f")) +
  ylim(c(0,100)) +
  geom_errorbar(aes(ymin = prop_dead, 
                    ymax = prop_dead + CI_half),
                color = "black",
                alpha = 1,
                width =.2,
                position = dodge) +
  xlab("") +
  ylab("") +
  facet_wrap(~ "Days >2", scales = "free_x") +
  guides(fill = guide_legend(title="Risk groups based on IPP")) + 
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 11, face = "bold"),
        legend.text = element_text(size=11)) +
  ggsignif::geom_signif(comparisons = list(c("3", "4")),
                        annotations = paste0("p=",signif(p2, digits=3)),
                        y_position = 95)

fig_5 = my_plot_1 + my_plot_2


#--------------------------------------------------------------------------------------#
####                                SUPPLEMENTARY MATERIAL                          ####
#--------------------------------------------------------------------------------------#

# Tables
sup_tab_3 = f_eval_train_test(df = df_cases_2, seed = 7, by_patient = T, prop = 0.7)

# Figures
sup_fig_1 = COCO_boxplot(studies)
sup_fig_2 = COCO_densityplot(studies, "CD3D")
sup_fig_3 = COCO_dotplot(studies, c("CEACAM1", "CLDN8"))

#--------------------------------------------------------------------------------------#
####                                      END                                       ####
#--------------------------------------------------------------------------------------#
