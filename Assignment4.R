# Read Data
library(survival)
IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))

# Cox Model without interaction of treatment and sex
coxfit1 <- coxph(Surv(time = survival, event = DIED, type = "right") 
                ~ AGE + I(GENDER=="M") + NEUTRO + PLATE + NIH,
                data = IMRAWIST)
summary(coxfit1)

# Cox Model with interaction of treatment and sex
coxfit2 <- coxph(Surv(time = survival, event = DIED, type = "right") 
                 ~ AGE + I(GENDER=="M") + NEUTRO + PLATE + NIH + I(GENDER=="M"):NIH,
                 data = IMRAWIST)
summary(coxfit2)
