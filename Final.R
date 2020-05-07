# Read and recenter Data
library(survival)
IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))
## Centering
IMRAWIST$AGE = scale(IMRAWIST$AGE, center = TRUE, scale = FALSE)
IMRAWIST$NEUTRO = scale(IMRAWIST$NEUTRO, center = TRUE, scale = FALSE)
IMRAWIST$PLATE = scale(IMRAWIST$PLATE, center = TRUE, scale = FALSE)
## Change unit of time to year
IMRAWIST$survival = IMRAWIST$survival / 365.25
IMRAWIST$AMLTIME  = IMRAWIST$AMLTIME  / 365.25
## Correct AML status for subjects 826 and 936
IMRAWIST$GENDER = as.integer(IMRAWIST$GENDER=="M")
#IMRAWIST$AML[826] = 1; IMRAWIST$AML[936] = 1
#IMRAWIST$BLASTS <- IMRAWIST$DIAGNOS <- IMRAWIST$IPSS <- IMRAWIST$cyto.s <- NULL

# Problem 1
library(survC1)
DataDeath = with(IMRAWIST, cbind(survival, DIED, AGE,GENDER, NEUTRO, PLATE, NIH))
DataDeath = CompCase(DataDeath) # Eliminate rows with NAs
# 1.1 & 1.2
Inf.Cval.Delta(DataDeath[,c(1,2)], covs0 = DataDeath[,-c(1,2,7)], covs1 = DataDeath[,-c(1,2)],
               tau = 3, itr = 200)
# 1.3 & 1.4
Inf.Cval.Delta(DataDeath[,c(1,2)], covs0 = DataDeath[,-c(1,2,7)], covs1 = DataDeath[,-c(1,2)],
               tau = 5, itr = 200)

# Problem 2
library(nricens)
cox1<-coxph(Surv(survival,DIED)~AGE+GENDER+NEUTRO+PLATE,
            data=IMRAWIST,x=TRUE) 
cox2<-coxph(Surv(survival,DIED)~AGE+GENDER+NEUTRO+PLATE+NIH,
            data=IMRAWIST,x=TRUE)
nri1 = nricens(mdl.std = cox1, mdl.new = cox2, t0 = 3,
        updown = "category", cut = c(0.1, 0.3))
nri2 = nricens(mdl.std = cox1, mdl.new = cox2, t0 = 3,
        updown = "category", cut = c(0.1, 0.3, 0.5))
nri3 = nricens(mdl.std = cox1, mdl.new = cox2, t0 = 5,
        updown = "category", cut = c(0.2, 0.4))
nri4 = nricens(mdl.std = cox1, mdl.new = cox2, t0 = 5,
        updown = "category", cut = c(0.2, 0.4, 0.6))

# Problem 3
library(cmprsk)
IMRAWIST$status = 0 #Event 0 Censored (survival)
IMRAWIST$status[which(IMRAWIST$AML==1 & IMRAWIST$DIED==1)] = 1 #Event 1 AML Death
IMRAWIST$status[which(IMRAWIST$AML==0 & IMRAWIST$DIED==1)] = 2 #Event 2 Non-AML Death
# IMRAWIST$status = factor(IMRAWIST$status, levels=c(0,1,2), labels= c("Censored", "Non-AML", "AML"))
IMRAWIST$NIH = factor(IMRAWIST$NIH, levels=c(0,1), labels= c("IMRAW", "IST"))

cmpfit = with(IMRAWIST, cuminc(ftime = survival, fstatus = status, group = NIH))
## 3.1
CumInc3 = with(timepoints(cmpfit, times = 3),{
   CritVal = qnorm(0.975)
   CI.dev = CritVal * sqrt(var/(log(est)*est)^2)
   est.loglog = log(-log(est))
   CI.u = exp(-exp(est.loglog - CI.dev))
   CI.l = exp(-exp(est.loglog + CI.dev))
   table = cbind(est, sqrt(var), CI.l, CI.u)
   colnames(table) = c("est","se","0.025","0.975") 
   rownames(table) = c("IMRAW:AML", "IST:AML", "IMRAW:NonAML", "IST:NonAML")
   table
})
CumInc3
## 3.2
CumInc5 = with(timepoints(cmpfit, times = 5),{
   CritVal = qnorm(0.975)
   CI.dev = CritVal * sqrt(var/(log(est)*est)^2)
   est.loglog = log(-log(est))
   CI.u = exp(-exp(est.loglog - CI.dev))
   CI.l = exp(-exp(est.loglog + CI.dev))
   table = cbind(est, sqrt(var), CI.l, CI.u)
   colnames(table) = c("est","se","0.025","0.975") 
   rownames(table) = c("IMRAW:AML", "IST:AML", "IMRAW:NonAML", "IST:NonAML")
   table
})
CumInc5

## 3.3
IMRAWIST = IMRAWIST[complete.cases(IMRAWIST[ , c("survival","AGE", "GENDER", "NEUTRO", "PLATE", "NIH")]),]
crrfit1 = with(IMRAWIST, crr(ftime = survival, fstatus = status, cov1 = cbind(AGE, GENDER, NEUTRO, PLATE, NIH),
                   failcode = 1))
coxfit3 = coxph(Surv(time = survival, event = as.numeric(status==1), type = "right") 
               ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST)
summary(coxfit3)
## 3.4
crrfit2 = with(IMRAWIST, crr(ftime = survival, fstatus = status, cov1 = cbind(AGE, GENDER, NEUTRO, PLATE, NIH),
                   failcode = 2))
coxfit4 = coxph(Surv(time = survival, event = as.numeric(status==2), type = "right") 
               ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST)
summary(coxfit4)
