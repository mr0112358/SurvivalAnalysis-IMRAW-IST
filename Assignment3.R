library(survival)
library(ggplot2)
library(data.table)

# Settings
n = 200
theta = 1 / 10
gamma = 1 / 12

# Simulation
S_all = vapply(1:500, function(B) {
  ## generate data
  set.seed(B)
  T = rexp(n, rate = theta)
  set.seed(B + 500)
  C = rexp(n, rate = gamma)
  
  delta = as.integer(T <= C)
  X = T * delta + C * (1 - delta)
  
  ## MLE
  theta_ML = sum(delta) / sum(X)
  S_ML = c(exp(-theta_ML * 8), exp(-theta_ML * 12), exp(-theta_ML * 16))
  ## KM
  S_KM = summary(survfit(Surv(
    time = X,
    event = delta,
    type = "right"
  ) ~ 1),
  type = "kaplan-meier",
  times = c(8, 12, 16))$surv
  ## combine MLE and KM into a vector
  c(S_ML, S_KM)
}, numeric(6))

# 1. Histograms

## Add tags
tags = cbind(
  c("S(8)", "S(12)", "S(16)", "S(8)", "S(12)", "S(16)"),
  c("MLE", "MLE", "MLE", "KM", "KM", "KM")
)
colnames(tags) = c("SurvFunc", "Method")
S_est = data.table(cbind(tags, S_all))
S_est_long = data.frame(melt(
  S_est,
  measure.vars = 3:502,
  variable.name = "COLS",
  value.name = "Estimate"
)[, COLS := NULL])
S_est_long$Estimate = as.numeric(S_est_long$Estimate)

ggplot(S_est_long, aes(x = Estimate, fill = Method)) +
  geom_histogram(color = "black", alpha = 0.5) +
  facet_grid(factor(Method, 
    levels = c("MLE", "KM")) ~ factor(SurvFunc, levels = c("S(8)", "S(12)", "S(16)")))

# 2. MSE
## Theoretic result
S = c(exp(-theta * 8), exp(-theta * 12), exp(-theta * 16))
S = c(S, S)
MSE = t(matrix(vapply(1:6, function(k) {
  mean((S_all[k, ] - S[k]) ^ 2)
}, numeric(1)), ncol = 2))
colnames(MSE) = c("S(8)", "S(12)", "S(16)")
rownames(MSE) = c("MLE", "KM")

print(xtable::xtable(MSE, digits = 5))
