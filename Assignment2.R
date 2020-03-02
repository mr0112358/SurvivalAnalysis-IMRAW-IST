# Read Data
library(doParallel)
library(survival)
library(xtable)
IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))
IMRAW = subset(IMRAWIST, NIH == 0)
IST = subset(IMRAWIST, NIH == 1)
IMRAWIST = list(IMRAW, IST)

# function to find indices of time for years 1, 3 and 5
findTime = function(survObj){
  year1 = sum(survObj$time <= 365)
  year3 = sum(survObj$time <= 1095)
  year5 = sum(survObj$time <= 1826) # At least one year is leap year
  c(year1, year3, year5)
}

registerDoParallel(detectCores())

# KM estimator for survival functions
KMtables = foreach(d = 1:2) %dopar% {
dataName = ifelse(d==1, "IMRAW", "IST")

KM_death = survfit(Surv(time = survival, event = DIED) ~ 1,
                    type = "kaplan-meier",
                    data = IMRAWIST[[d]])

SurvBoot = foreach(B = 1:500, .combine = 'rbind') %dopar% {
n = nrow(IMRAWIST[[d]])
set.seed(B)
bootIndex = sample.int(n, replace = TRUE)
KM_death_boot = survfit(Surv(time = survival, event = DIED) ~ 1,
                    type = "kaplan-meier",
                    data = IMRAWIST[[d]][bootIndex,])

yearsBoot = findTime(KM_death_boot)
KM_death_boot$surv[yearsBoot]
}

years = findTime(KM_death)

TableKM = t(rbind(KM_death$surv[years],
  KM_death$lower[years],
  KM_death$upper[years],
  bootCI = sapply(1:3, function(y) {quantile(SurvBoot[,y], probs = c(0.025, 0.975))})))
rownames(TableKM) = c("year 1","year 3","year 5")
colnames(TableKM) = c("hat S(t)", "2.5%", "97.5%", "2.5% Boot", "97.5% Boot")
xtable(TableKM,digits = 4, 
             caption = paste0("KM estimates of $S(1)$, $S(2)$, $S(5)$
             with 95% CIs by Greenwood formula and bootstrap for ", dataName, " patients"))
}


## NA estimator for cumulative hazard function
NAtables = foreach(d = 1:2) %dopar% {
dataName = ifelse(d==1, "IMRAW", "IST")

NA_death = survfit(Surv(time = survival, event = DIED) ~ 1,
                    type = "fleming-harrington",
                    data = IMRAWIST[[d]])

CumHazardBoot = foreach(B = 1:500, .combine = 'rbind') %dopar% {
n = nrow(IMRAWIST[[d]])
set.seed(B)
bootIndex = sample.int(n, replace = TRUE)
NA_death_boot = survfit(Surv(time = survival, event = DIED) ~ 1,
                    type = "fleming-harrington",
                    data = IMRAWIST[[d]][bootIndex,])

yearsBoot = findTime(NA_death_boot)
NA_death_boot$cumhaz[yearsBoot]
}

years = findTime(NA_death)

TableNA = t(rbind(NA_death$cumhaz[years],
  -log(NA_death$upper[years]),
  -log(NA_death$lower[years]),
  bootCI = sapply(1:3, function(y) {quantile(CumHazardBoot[,y], probs = c(0.025, 0.975))})))
rownames(TableNA) = c("year 1","year 3","year 5")
colnames(TableNA) = c("hat S(t)", "2.5%", "97.5%", "2.5% Boot", "97.5% Boot")
xtable(TableNA,digits = 4, 
             caption = paste0("NA estimates of $\\Lambda(1)$, $\\Lambda(2)$, $\\Lambda(5)$ 
             with 95% CIs by Greenwood formula and bootstrap for ", dataName, " patients"))
}

print(KMtables[[1]])
print(KMtables[[2]])
print(NAtables[[1]])
print(NAtables[[2]])
