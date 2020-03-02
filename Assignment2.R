# Read Data
library(doParallel)
library(survival)
library(xtable)
IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))
IMRAW = subset(IMRAWIST, NIH == 0)
IST = subset(IMRAWIST, NIH == 1)
IMRAWIST = list(IMRAW, IST)

registerDoParallel(detectCores())

# KM estimator for survival functions
KMtables = foreach(d = 1:2) %dopar% {
  dataName = ifelse(d == 1, "IMRAW", "IST")
  
  KM_death = survfit(
    Surv(time = survival, event = DIED) ~ 1,
    type = "kaplan-meier",
    error = "greenwood",
    data = IMRAWIST[[d]]
  )
  
  SurvBoot = foreach(B = 1:500, .combine = 'rbind') %dopar% {
    n = nrow(IMRAWIST[[d]])
    set.seed(B)
    bootIndex = sample.int(n, replace = TRUE)
    KM_death_boot = survfit(Surv(time = survival, event = DIED) ~ 1,
                            type = "kaplan-meier",
                            data = IMRAWIST[[d]][bootIndex, ])
    KM_death_boot = summary(KM_death_boot, times = c(1, 3, 5), scale = 365.25)
    KM_death_boot$surv
  }
  
  KM_death = summary(KM_death, times = c(1, 3, 5), scale = 365.25)
  
  TableKM = t(rbind(
    KM_death$surv,
    KM_death$lower,
    KM_death$upper,
    bootCI = sapply(1:3, function(y) {
      quantile(SurvBoot[, y], probs = c(0.025, 0.975))
    })
  ))
  rownames(TableKM) = c("year 1", "year 3", "year 5")
  colnames(TableKM) = c("$\\hS(t)$", "2.5%", "97.5%", "2.5% Boot", "97.5% Boot")
  xtable(
    TableKM,
    digits = 4,
    caption = paste0(
      "KM estimates of $S(1)$, $S(2)$, $S(5)$
             with 95% CIs by Greenwood formula and bootstrap for ",
      dataName,
      " patients"
    )
  )
}


# NA estimator for cumulative hazard function
NAtables = foreach(d = 1:2) %dopar% {
  dataName = ifelse(d == 1, "IMRAW", "IST")
  
  NA_death = survfit(Surv(time = survival, event = DIED) ~ 1,
                     type = "fleming-harrington",
                     data = IMRAWIST[[d]])
  
  CumHazardBoot = foreach(B = 1:500, .combine = 'rbind') %dopar% {
    n = nrow(IMRAWIST[[d]])
    set.seed(B)
    bootIndex = sample.int(n, replace = TRUE)
    NA_death_boot = survfit(Surv(time = survival, event = DIED) ~ 1,
                            type = "fleming-harrington",
                            data = IMRAWIST[[d]][bootIndex, ])
    NA_death_boot = summary(NA_death_boot, times = c(1, 3, 5), scale = 365.25)
    NA_death_boot$cumhaz
  }
  
  NA_death = summary(NA_death, times = c(1, 3, 5), scale = 365.25)
  
  TableNA = t(rbind(
    NA_death$cumhaz,
    -log(NA_death$upper),
    -log(NA_death$lower),
    bootCI = sapply(1:3, function(y) {
      quantile(CumHazardBoot[, y], probs = c(0.025, 0.975))
    })
  ))
  rownames(TableNA) = c("year 1", "year 3", "year 5")
  colnames(TableNA) = c("$\\hLam(t)$", "2.5%", "97.5%", "2.5% Boot", "97.5% Boot")
  xtable(
    TableNA,
    digits = 4,
    caption = paste0(
      "NA estimates of $\\Lambda(1)$, $\\Lambda(2)$, $\\Lambda(5)$
             with 95% CIs by asymptotics and bootstrap for ",
      dataName,
      " patients"
    )
  )
}

print(KMtables[[1]])
print(KMtables[[2]])
print(NAtables[[1]])
print(NAtables[[2]])
