library(readxl)
library(tidyverse)
library(reshape)

## Figure 1b. Summary of timepoints for each patient
data = data.frame(read_excel("CLL_Patient_Data.xlsx", sheet = 1))

makeGraphs = function(df = NA, title = "", name = ""){
  fixed_column = c()
  for (sample_timepoints in df$timepoints){
    final_timepoint_vector = c("0")
    numerical_times = grep("Y", strsplit(sample_timepoints, ", ")[[1]], value = TRUE)
    final_timepoint_vector = c(final_timepoint_vector, str_sub(numerical_times, 2, -1))
    fixed_column = c(fixed_column, paste(final_timepoint_vector, collapse = ","))
  }
  df$timepoints = fixed_column
  df = subset(df, select=c("patient", "timepoints", "Response_Richters"))
  df = separate_rows(df, timepoints, sep=",")
  if (name == "ibrutinib"){
    clrs = c("firebrick2", "darkred", "royalblue3")
  }
  else{clrs = c("firebrick2", "royalblue3")}
  ggplot(df, aes(x=as.numeric(timepoints), y=patient)) + 
    geom_point(size=3, shape=16, aes(color=Response_Richters)) +
    geom_line(aes(color=Response_Richters)) +
    xlab("Time at sampling (Years)") + 
    ylab("Patient ID") +
    xlim(0,8) +
    ggtitle(title) + 
    scale_color_manual(values = clrs) +
    labs(color='Response') +
    theme_bw()
  print(pl)
  if (name == "ibrutinib"){
    ggsave(paste(name,"_timepoints.png", sep = ""), plot=pl, device = "png", width = 7.45, height = 6.5)
  }
  else{
    ggsave(paste(name,"_timepoints.png", sep = ""), plot=pl, device = "png", width = 7, height = 6.5)
  }
}

ibrutinib_order=c("250422","250733","251150","250801","251853",
                  "251646","251940","252336","EJW0719-1","EJW1940-4",
                  "111330026","111330031","111330025","111330019",
                  "111330111","251218","111330087","111330103",
                  "111330105","111330126","111330009","111330082")
acalabrutinib_order=c("250635","251623","251742","251936",
                      "250859","250890","251972","251363","251497","251572",
                      "251786","251863","251882","250738","251185","EJW1856-2","EJW1983-3")

ibru = filter(data, treatment=='ibrutinib')
acala = filter(data, treatment=="acalabrutinib")
ibru$patient = factor(ibru$patient, levels=ibrutinib_order)
acala$patient = factor(acala$patient, levels=acalabrutinib_order)

makeGraphs(acala, "Acalabrutinib Cohort Sampling Time Points", "acalabrutinib")
makeGraphs(ibru, "Ibrutinib Cohort Sampling Time Points", "ibrutinib")


###############################################
## Treatment vs response (Table 1)
table(df$treatment, df$response)
fisher.test(df$treatment, df$response)

###############################################
###  Kaplan Meier survival(relapse) curves  ###
###############################################
library(survminer)
require("survival")
library(readxl)
library(tidyverse)


data = data.frame(read_excel("CLL_Patient_Data.xlsx", sheet = 1))
acala = filter(data, treatment == "acalabrutinib")
ibru = filter(data, treatment == "ibrutinib")
no53 = filter(data, tp53_mut == "no")
noRich = filter(data, Richters == "no")
df = no53
df = data
df = ibru
df = acala

## Evolution of subclones that include any mutations ##
## Figure 3a ##
df = data
df = ibru
df = acala 
fit <- survfit(Surv(treatment_time, status) ~ evolution_no53, data = df)
p <- ggsurvplot(fit,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years", ylab = "% Relapse free",
                legend.labs = c("Evolution", "No Evolution", "TP53 mutation"),
                palette = c("black", "navy", "red3"),
                title = "Acalabrutinib patients only",
                ggtheme = theme_light())

## Exclude patients with TP53 mutations, then compare those with and without evolving subclones to get Log-Rank P value) ##
no_53 = filter(df, evolution_no53 != "TP53")
fit <- survfit(Surv(treatment_time, status) ~ evolution_no53, data = no_53)
p <- ggsurvplot(fit,
                data = no_53,
                pval = TRUE,
                pval.method = TRUE,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years", ylab = "% Relapse free",
                ggtheme = theme_light())

## CLL-driver evolution (Figure 3b) ##
fit <- survfit(Surv(treatment_time, status) ~ evolving_CLL_no53, data = df)
p <- ggsurvplot(fit,
           break.time.by = 0.5,
           risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
           xlab = "Time in years", ylab = "% Relapse free",
           legend.labs = c("No CLL-driver evolution", "TP53 mutation", "CLL-driver evolution"),
           palette = c("navy", "red3", "black"),
           title = "All patients",
           ggtheme = theme_light())

## Exclude patients with TP53 mutations, then compare those with and without CLL-driver genes in evolving subclones to get Log-Rank P value) ##
no_53 = filter(df, evolution_no53 != "TP53")
fit <- survfit(Surv(treatment_time, status) ~ evolving_CLL_no53, data = no_53)
p <- ggsurvplot(fit,
                data = no_53,
                pval = TRUE,
                pval.method = TRUE,
                break.time.by = 0.5,
                risk.table = TRUE,risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
                xlab = "Time in years", ylab = "% Relapse free",
                ggtheme = theme_light())

## Comparing the two drugs ##
## Supplemental figure 1
df = data
fit <- survfit(Surv(treatment_time, status) ~ treatment, data = df)
p = ggsurvplot(fit,data = df,
               pval = TRUE,pval.method = TRUE,
               break.time.by = 0.5,
               risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
               xlab = "Time in years",ylab = "% Relapse free",
               legend.labs = c("Acalabrutinib", "Ibrutinib"),
               palette = c("navy", "red3"),
               ggtheme = theme_light())

## Exclude patients with baseline TP53 mutations ##
df = no53
fit <- survfit(Surv(treatment_time, status) ~ treatment, data = df)
p = ggsurvplot(fit,data = df,
               pval = TRUE,pval.method = TRUE,
               break.time.by = 0.5,
               risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
               xlab = "Time in years",ylab = "% Relapse free",
               legend.labs = c("Acalabrutinib", "Ibrutinib"),
               palette = c("navy", "red3"),
               ggtheme = theme_light())

## Exclude patients who developed Richter's transformation ##
df = noRich
fit <- survfit(Surv(treatment_time, status) ~ treatment, data = df)
p = ggsurvplot(fit,data = df,
               pval = TRUE,pval.method = TRUE,
               break.time.by = 0.5,
               risk.table = TRUE, risk.table.fontsize = 2.75,risk.table.height = .2,risk.table.y.text = FALSE,
               xlab = "Time in years",ylab = "% Relapse free",
               legend.labs = c("Acalabrutinib", "Ibrutinib"),
               palette = c("navy", "red3"),
               ggtheme = theme_light())

###################################################
## Comparing CNVs to TP53 mutation
table(df$CNV_1, df$tp53_mut)
fisher.test(df$CNV_1, df$tp53_mut)

table(df$tp53_mut, df$evolution_status)
fisher.test(df$tp53_mut, df$evolution_status)

use = filter(df, evolution_type == "no_effect")
df = use
table(df$response, df$tp53_mut)
fisher.test(df$response, df$tp53_mut)

####################################
### Other analyses used in paper ###

## CLL specific evolution vs outcome in patients with evolution
df = data.frame(read_excel("CLL_Patient_Data.xlsx", sheet = 1))
no_53 = filter(df, tp53_mut == "no")
ibru = filter(no_53, evolution_status == "evolution")
acala = filter(no_53, treatment == "acalabrutinib")

# Evolving CLL genes vs treatment outcomes.
table(df$response, df$evolving_CLL_gene)
fisher.test(df$response, df$evolving_CLL_gene)
# Filtering out TP53
table(no_53$response, no_53$evolving_CLL_gene)
fisher.test(no_53$response, no_53$evolving_CLL_gene)
# Evolving CLL genes vs outcomes within each treatment.
table(ibru$response, ibru$evolving_CLL_gene)
fisher.test(ibru$response, ibru$evolving_CLL_gene)
table(acala$response, acala$evolving_CLL_gene)
fisher.test(acala$response, acala$evolving_CLL_gene)

# Evolving CLL Genes regardless of genes mutated.
table(df$response, df$evolution_status)
fisher.test(df$response, df$evolution_status)
table(no_53$response, no_53$evolution_status)
fisher.test(no_53$response, no_53$evolution_status)
table(use$treatment, use$evolving_CLL_gene)
fisher.test(use$treatment, use$evolving_CLL_gene)

###################################################################################################
### Make graphs for number of CNVs on the X and patient on the Y
sorted = df[order(df$num_cnvs),]
tmp = as.character(factor(sorted$patient))
sorted$patient = factor(sorted$patient, levels=tmp)
ggplot(sorted, aes(x=num_cnvs, y=patient)) +
  geom_point(size=4, shape=16, aes(color=response)) +
  xlab("Number of CNVs") + 
  ylab("Patient ID") +
  xlim(0,5) +
  theme_bw()

ggplot(sorted, aes(x=num_cnvs, y=patient)) +
  geom_point(size=4, shape=16, aes(color=evolution_type)) +
  xlab("Number of CNVs") + 
  ylab("Patient ID") +
  xlim(0,14) +
  theme_bw()

###################################################################################################
## Plotting the first evolution event

data = data.frame(read_excel("CLL_Patient_Data.xlsx", sheet = 1))
ibru = filter(data, treatment=='ibrutinib')
acala = filter(data, treatment=="acalabrutinib")

df = filter(data, evolution_status=="evolution")
df = filter(ibru, evolution_status=="evolution")
df = filter(acala, evolution_status=="evolution")
df = filter(data, evolving_CLL_gene=="yes")

ggplot(df, aes(x=as.numeric(year_at_evolution), y=patient)) + 
  geom_point(size=3, shape=16, aes(color=response)) +
  xlab("Time at sampling (Years)") + 
  ylab("Patient ID") +
  xlim(0,3) +
  ggtitle("Year at first evolution") + 
  theme_bw()

ggplot(df, aes(x=as.numeric(year_at_evolution))) + 
  geom_bar(width = 0.2, aes(fill=response)) +
  guides(fill=FALSE) +
  facet_grid(cols = vars(response)) +
  scale_fill_manual(values = c("red3", "navy")) +
  scale_fill_manual(values = c("firebrick3", "royalblue4")) +
  xlab("Years after treatment began") + 
  ylab("Number of patients having first evolution event") +
  xlim(0,3.5) +
  ggtitle("First evolving clone with CLL-driver") +  
  theme_bw()