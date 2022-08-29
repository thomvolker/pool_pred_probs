
###############################################################################
## Required packages
###############################################################################

library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library(furrr)
library(mice)

###############################################################################
## Simulation functions
###############################################################################

# Create the correlation matrix for all variables
get_cormat <- function(cor, vars) {
  cormat <- diag(length(vars))
  cormat[cormat == 0] <- cor
  as.matrix(cormat)
}

# Make the regression coefficients given a prespecified effect size (specified in
# data simulation function)
get_coefs <- function(var_yhat, rel_weight, cormat) {
  sqrt(var_yhat / sum(rel_weight %*% t(rel_weight) * cormat)) * rel_weight
}

# Calculate an intercept based on a known prevalence (to vary class imbalance)
get_intercept <- function(betas, cormat, prevalence) {
  
  sigma <- c(sqrt(t(betas) %*% cormat %*% betas))
  
  integral <- function(x, b0) (1/sqrt(2 * pi)) * exp(-x^2 / 2) / (1 + exp(-b0 - sigma * x))
  
  f_integrate <- function(g, prevalence, g_lower, g_upper, ...) {
    prevalence - integrate(f = g, lower = g_lower, upper = g_upper, ...)$value
  }
  
  uniroot(f = f_integrate,
          interval = c(-99999999, 99999999),
          prevalence = prevalence,
          g = integral,
          g_lower = -Inf,
          g_upper = Inf)$root
}

# Simulate the complete samples
sim_dat <- function(n, r2, prevalence = 0.5, rel_weight, cormat) {
  var_yhat <- (r2 * pi^2/3) / (1-r2)
  coefs <- get_coefs(var_yhat, rel_weight, cormat)
  b0 <- get_intercept(coefs, cormat, prevalence = prevalence)
  
  X <- MASS::mvrnorm(n, mu = rep(0, length(coefs)), Sigma = cormat)
  Y <- rbinom(n, 1, 1 / (1 + exp(-(b0 + X %*% coefs))))
  
  bind_cols(X = data.frame(X), Y = Y)
}

###############################################################################
## Test simulation functions
###############################################################################

set.seed(123)
relative.strength <- 1:5
cormat            <- get_cormat(0.3, relative.strength)
df                <- sim_dat(n = 100000,                     # sample size
                             r2 = 0.25,                      # McKelvey and Zavoina's R^2
                             prevalence = 0.2,               # prevalence of the outcome
                             rel_weight = relative.strength, # relative strength of the coefficients
                             cormat = cormat)

mean(df$Y) # about .20

fit <- glm(Y ~ ., data = df, family = binomial) 

summary(fit) # relative strength of coefficients makes sense
performance::r2_mckelvey(fit) # R^2 about 0.25


###############################################################################
## Run simulation
###############################################################################

nsim <- 10

n    <- c(50, 250, 500)
r2   <- c(0.15, 0.5, 0.85)
prev <- c(0.20, 0.35, 0.50)

plan(multisession)

sim <- future_map_dfr(1:nsim, 
                      ~expand_grid(N = n, R2 = r2, Prevalence = prev, 
                                   relative.strength = list(relative.strength),
                                   cormat = list(cormat)) %>%
                        mutate(dat = pmap(list(N, R2, Prevalence, relative.strength, cormat), 
                                          function(N, R2, Prevalence, relative.strength, cormat) {
                                            sim_dat(n = N, r2 = R2, prevalence = Prevalence,
                                                    rel_weight = relative.strength, cormat = cormat)
                                          }),
                               amp = map(dat, ~ampute(.x)), # probably specify other arguments here
                               imp = map(amp, ~mice(.x$amp, m = 5, print=F))), # and specify other arguments here as well
                      .id = "NSIM",
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)

saveRDS(sim, file = "results//pool_pred_probs.RDS")


