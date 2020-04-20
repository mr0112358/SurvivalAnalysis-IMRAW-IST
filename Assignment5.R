# Read and centralize data
library(survival)
IMRAWIST = data.frame(readxl::read_excel("IMRAWandISTnov.xls"))
IMRAWIST$AGE = scale(IMRAWIST$AGE, center = TRUE, scale = FALSE)
IMRAWIST$NEUTRO = scale(IMRAWIST$NEUTRO, center = TRUE, scale = FALSE)
IMRAWIST$PLATE = scale(IMRAWIST$PLATE, center = TRUE, scale = FALSE)

# Function for summary table
sumtable = function(survregobj){
  survregSum = summary(survregobj)
  table = cbind(survregSum$table[2:6,1], survregSum$table[2:6,4], 
    exp(survregSum$table[2:6,1]/survregobj$scale), 
    exp((survregSum$table[2:6,1] - 1.96*survregSum$table[2:6,2])/survregobj$scale), 
    exp((survregSum$table[2:6,1] + 1.96*survregSum$table[2:6,2])/survregobj$scale))
  colnames(table) = c("estimate", "p-value", "hazard ratio", "lower CI", "upper CI")
  table
}

# Cox Model for exponential hazard baseline
regfit1 <- survreg(Surv(time = survival, event = DIED, type = "right") 
                ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST,
                dist = "exponential")
sumtable(regfit1)

# Cox Model for Weibull hazard baseline
regfit2 <- survreg(Surv(time = survival, event = DIED, type = "right") 
                ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST,
                dist = "weibull")
sumtable(regfit2)

# Cox Model for log-Normal hazard baseline
regfit3 <- survreg(Surv(time = survival, event = DIED, type = "right") 
                ~ AGE + GENDER + NEUTRO + PLATE + NIH,
                data = IMRAWIST,
                dist = "lognormal")
sumtable(regfit3)
