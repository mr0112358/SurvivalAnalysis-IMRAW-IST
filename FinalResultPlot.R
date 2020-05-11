library(survival)
library(survminer)

IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))
## Change unit of time to year
IMRAWIST$survival = IMRAWIST$survival / 365.25
IMRAWIST$AMLTIME  = IMRAWIST$AMLTIME  / 365.25
IMRAWIST$TREATMENT = factor(IMRAWIST$NIH, levels=c(0,1), labels= c("IMRAW", "IST"))

#######################################################
#                     Tables                          #
#######################################################
# Summary Table
SumStat = function(DF){
  No = nrow(DF)
  DthRate = sum(DF$DIED == 1) / No * 100
  ANC = paste0(sprintf("%.1f",mean(DF$NEUTRO, na.rm = T)), "pm" , sprintf("%.1f",1.96 * sd(DF$NEUTRO, na.rm = T)/sqrt(length(complete.cases(DF$NEUTRO)))))
  Platelets = paste0(sprintf("%.1f",mean(DF$PLATE, na.rm = T)), "pm" , sprintf("%.1f",1.96 * sd(DF$PLATE, na.rm = T)/sqrt(length(complete.cases(DF$PLATE)))))
  data.frame(No, DthRate, ANC, Platelets)
}

SumStat2 = function(DF){
  Tab = rbind(
  SumStat(DF),
  SumStat(subset(DF,GENDER=="F")),
  SumStat(subset(DF,GENDER=="M")),
  SumStat(subset(DF,AGE <= 60)),
  SumStat(subset(DF,AGE > 60))
  )
  rownames(Tab) = c("Total","Female","Male","Age leq 60", "Age >60")
  Tab
}

Total_IMRAW_Sum = SumStat2(subset(IMRAWIST, NIH == "IMRAW"))
Total_IST_Sum = SumStat2(subset(IMRAWIST, NIH == "IST"))
Total_Sum = cbind(Total_IMRAW_Sum,Total_IST_Sum)
AML_IMRAW_Sum = SumStat2(subset(IMRAWIST, NIH == "IMRAW" & AML == 1))
AML_IST_Sum = SumStat2(subset(IMRAWIST, NIH == "IST" & AML == 1))
AML_Sum = cbind(AML_IMRAW_Sum,AML_IST_Sum)

xtable::xtable(rbind(Total_Sum,AML_Sum), digits=1)

# CoxPH Total death
coxfit <- coxph(Surv(time = survival, event = DIED, type = "right") 
                ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST)
coxsum = summary(coxfit)

coxTab = cbind(coxsum$coefficients[,c(1,2)], coxsum$conf.int[,c(3,4)], coxsum$coefficients[,5])

xtable::xtable(coxTab, digits = 3)


#######################################################
#                    Figures                          #
#######################################################
# KM Survival Curve
DEATH_KM = survfit(Surv(time = survival, event = DIED) ~ TREATMENT, 
                    type = "kaplan-meier",
                    data = IMRAWIST)
AML_KM = survfit(Surv(time = AMLTIME, event = AML) ~ TREATMENT, 
                    type = "kaplan-meier",
                    data = IMRAWIST)

KMplots <- list()
KMplots[[1]] <- ggsurvplot(DEATH_KM, data = IMRAWIST,
                           xlim = c(0,12),
                           break.x.by=3,
                           xlab = "Time (years)",
                           font.title = c(10, "bold", "black"),
                           title = "Survival Curves of Time to Death",
                           legend.title = "Treatment",
                           legend.labs = c("IMRAW","IST"),
                           #surv.median.line = "hv",
                          risk.table = TRUE,
                          risk.table.fontsize = 3,
                          tables.y.text = FALSE,
                          conf.int = TRUE,
                          ggtheme = theme_light())
KMplots[[2]] <- ggsurvplot(AML_KM, data = IMRAWIST,
                           xlim = c(0,12),
                           break.x.by=3,
                           xlab = "Time (years)",
                           font.title = c(10, "bold", "black"),
                           title = "Survival Curves of Time to AML",
                           legend.title = "Treatment",
                           legend.labs = c("IMRAW","IST"),
                           #surv.median.line = "hv",
                          risk.table = TRUE,
                          risk.table.fontsize = 3,
                          tables.y.text = FALSE,
                          conf.int = TRUE,
                          ggtheme = theme_light())
KM = arrange_ggsurvplots(KMplots, print = FALSE,
                    ncol = 2, nrow = 1, risk.table.height = 0.25)
ggsave(KM, width = 6.2, height = 4, device = "pdf", filename = "KM.pdf")

# Competing Risk Curves
library(cmprsk)
IMRAWIST$status = 0 #Event 0 Censored (survival)
IMRAWIST$status[which(IMRAWIST$AML==1 & IMRAWIST$DIED==1)] = 1 #Event 1 AML Death
IMRAWIST$status[which(IMRAWIST$AML==0 & IMRAWIST$DIED==1)] = 2 #Event 2 Non-AML Death
IMRAWIST$status = factor(IMRAWIST$status, levels=c(0,1,2), labels= c("Censored", "AML Death", "Non-AML Death"))
IMRAWIST$NIH = factor(IMRAWIST$NIH, levels=c(0,1), labels= c("IMRAW", "IST"))

cmpfit = with(IMRAWIST, cuminc(ftime = survival, fstatus = status, group = NIH, cencode = "Censored"))
cmp = ggcompetingrisks(cmpfit,
                 xlim = c(0.5,11.7),
                 xlab = "Time (years)",
                 #legend.title = "Event",
                 #legend.labs = c("AML Death","Non-AML Death"),
                 conf.int = TRUE,
                 title = "Cumulative Incidence Functions of AML and Non-AML Deaths",
                 legend = "top",
                 ggtheme = theme_bw())
ggsave(cmp, width = 6.2, height = 3.5, device = "pdf", filename = "cmp.pdf")
