#####################################################################
## All of the commands used to generate the figures of this paper. ##
## Using the included Patient_Data_Results.xlsx file, this file    ##
## can be ran to reproduce each figure.                            ##
#####################################################################

library(tidyverse)
library(readxl)
library(survminer)
require("survival")

data = data.frame(read_excel("Patient_Data_Results.xlsx", sheet = 1))#, row.names = 1)
acala = filter(data, treatment == "acalabrutinib")
ibru = filter(data, treatment == "ibrutinib")
no_53 = filter(data, tp53_mutation == "no")
#################################################  Figure 1 ######################################################
## Figure 1c
df = data
table(df$response, df$evolving_CLL_gene)
fisher.test(df$response, df$evolving_CLL_gene)

## Same as figure 1c, but excluding patients with baseline TP53. 
## Not in figure, but used in the same paragraph of the manuscript.
table(no_53$response, no_53$evolving_CLL_gene)
fisher.test(no_53$response, no_53$evolving_CLL_gene)

#################################################  Figure 2 ######################################################
## KM Survival of All patients, Landmark set at 2 years. 
df = data
lm_dat = 
  df %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)

fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time after year 2 landmark", ylab = "Relapse-free probability",
                legend.labs = c("No CLL-driver evolution", "TP53 mutation", "CLL-driver evolution"),
                palette = c("navy", "red3", "black"),
                title = "All patients",
                ggtheme = theme_light())
p

## Replicate excluding patients with p53 for P value ##
lm_dat = 
  no_53 %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)
fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                pval = TRUE,pval.method = TRUE,
                break.time.by = 0.5,
                risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years",ylab = "Relapse-free probability",
                legend.labs = c("No Cll-driver evolution", "CLL-driver evolution"),
                palette = c("navy", "red3"),
                ggtheme = theme_light())
p

######################################## Supplemental Figure 1b #################################################
df = data
makeGraphs = function(df = NA, title = "", name = ""){
  fixed_column = c()
  for (sample_timepoints in df$timepoints){
    final_timepoint_vector = c("0")
    numerical_times = grep("Y", strsplit(sample_timepoints, ", ")[[1]], value = TRUE)
    final_timepoint_vector = c(final_timepoint_vector, str_sub(numerical_times, 2, -1))
    fixed_column = c(fixed_column, paste(final_timepoint_vector, collapse = ","))
  }
  df$timepoints = fixed_column
  df = subset(df, select=c("patient", "timepoints", "response_richters"))
  df = separate_rows(df, timepoints, sep=",")
  if (name == "ibrutinib"){
    clrs = c("firebrick2", "darkred", "royalblue3")
  }
  else{clrs = c("firebrick2", "royalblue3")}
  pl = ggplot(df, aes(x=as.numeric(timepoints), y=patient)) + 
    geom_point(size=3, shape=16, aes(color=response_richters)) +
    geom_line(aes(color=response_richters)) +
    xlab("Time at sampling (Years)") + 
    ylab("Patient ID") +
    xlim(0,8) +
    ggtitle(title) + 
    scale_color_manual(values = clrs) +
    labs(color='Response') +
    theme_bw() + 
    theme(axis.text = element_text(size = 9.2)) + 
    theme(legend.text = element_text(size = 9.2))
  print(pl)
}

ibrutinib_order=c("1","2","4","3","7",
                  "6","8","9","21",
                  "13","14","12","11",
                  "19","5","16","17",
                  "18","20","10","15")
acalabrutinib_order=c("22","30","31","35",
                      "24","25","36","27","28","29",
                      "32","33","34","23","26","37","38")

ibru$patient = factor(ibru$patient, levels=ibrutinib_order)
acala$patient = factor(acala$patient, levels=acalabrutinib_order)

makeGraphs(acala, "Acalabrutinib Cohort Sampling Time Points", "acalabrutinib")
makeGraphs(ibru, "Ibrutinib Cohort Sampling Time Points", "ibrutinib")

######################################## Supplemental Figure 2 #################################################
## Comparing the two drugs ##
## Supplemental Figure 2a ##
df = data
fit <- survfit(Surv(treatment_time, response_status) ~ treatment, data = df)
p = ggsurvplot(fit,data = df,
               pval = TRUE,pval.method = TRUE,
               break.time.by = 0.5,
               risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
               xlab = "Time in years",ylab = "Relapse-free probability",
               legend.labs = c("Acalabrutinib", "Ibrutinib"),
               palette = c("navy", "red3"),
               ggtheme = theme_light())
p

## Supplemental Figure 2b ##
df = no_53
fit <- survfit(Surv(treatment_time, response_status) ~ treatment, data = df)
p = ggsurvplot(fit,data = df,
               pval = TRUE,pval.method = TRUE,
               break.time.by = 0.5,
               risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
               xlab = "Time in years",ylab = "Relapse-free probability",
               legend.labs = c("Acalabrutinib", "Ibrutinib"),
               palette = c("navy", "red3"),
               ggtheme = theme_light())
p

############################################ Supplemental Figure 3 #################################################
## A curve for TP53 as strata. ##
## Supplemental Figure 3a
df = data
p53 = filter(df, tp53_mutation == "yes")
fit <- survfit(Surv(treatment_time, response_status) ~ tp53_mutation, data = p53)
p <- ggsurvplot(fit,
                data = p53,
                conf.int = FALSE,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 3.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years",ylab = "Relapse-free probability",
                legend = "none",
                palette = c("red3"),
                title = "Patients with baseline TP53 mutations",
                ggtheme = theme_light())
p

## Supplemental Figure 3b
table(df$response, df$tp53_mutation)
fisher.test(df$response, df$tp53_mutation)

############################################ Supplemental Figure 4 #################################################
## KM Survival of IBRUTINIB patients ONLY, Landmark set at 2 years. ##
## SUPPLEMENTAL FIGURE 4A. ##
df = ibru
lm_dat = 
  df %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)

fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time after year 2 landmark", ylab = "Relapse-free probability",
                legend.labs = c("No CLL-driver evolution", "TP53 mutation", "CLL-driver evolution"),
                palette = c("navy", "red3", "black"),
                title = "Ibrutinib patients only",
                ggtheme = theme_light())
p

## Replicate IBRUTINIB excluding patients with p53 for P value ##
df = filter(ibru, tp53_mutation == "no")
lm_dat = 
  df %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)
fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                pval = TRUE,pval.method = TRUE,
                break.time.by = 0.5,
                risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years",ylab = "Relapse-free probability",
                legend.labs = c("No Cll-driver evolution", "CLL-driver evolution"),
                palette = c("navy", "red3"),
                title = "Ibrutinib patients only",
                ggtheme = theme_light())
p

## KM Survival of ACALABRUTINIB patients ONLY, Landmark set at 2 years. ##
## SUPPLEMENTAL FIGURE 4B. ##
df = acala
lm_dat = 
  df %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)

fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time after year 2 landmark", ylab = "Relapse-free probability",
                legend.labs = c("No CLL-driver evolution", "TP53 mutation", "CLL-driver evolution"),
                palette = c("navy", "red3", "black"),
                title = "Acalabrutinib patients only",
                ggtheme = theme_light())
p

## Replicate ACALABRUTINIB excluding patients with p53 for P value ##
df = filter(acala, tp53_mutation == "no")
lm_dat = 
  df %>%
  filter(treatment_time >= 2)
lm_dat = 
  lm_dat %>%
  mutate(lm_time = treatment_time - 2)
fit <- survfit(Surv(lm_time, response_status) ~ evol_CLL_y2, data = lm_dat)
p <- ggsurvplot(fit,
                pval = TRUE,pval.method = TRUE,
                break.time.by = 0.5,
                risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years",ylab = "Relapse-free probability",
                legend.labs = c("No Cll-driver evolution", "CLL-driver evolution"),
                palette = c("navy", "red3"),
                title = "Acalabrutinib patients only",
                ggtheme = theme_light())
p