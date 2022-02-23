# Example of Bayesian Mediation Analysis via Nimble.

# This is a proof-of-concept; my actual analysis requires Bayesian analysis,
# so I won't have a frequentist model to compare with.

library(dplyr)
library(ggplot2)
library(tidyr)
library(nimble)

# Data Simulation -----------------------------------------
x <- runif(100, 0, 10)
m <- 2 + x + rnorm(length(x), 0, 2)
y <- 5 + 2*x + 0.5*m + rnorm(length(x), 0, 2)
# Indirect effect: 0.5
# Total effect: 2.5?

ggplot() + aes(x = x, y = y) + geom_point()

# Frequentist version -------------------------------------
mx <- lm(m ~ x)
yx <- lm(y ~ x)
ymx <- lm(y ~ m + x)

ab <- coef(mx)["x"]*coef(ymx)["m"]
print(paste0("Indirect effect (a*b):         ", round(ab, 3)))
bbprime <- coef(yx)["x"] - coef(ymx)["x"]
print(paste0("Indirect effect (b - b prime): ", round(bbprime, 3)))
print(paste0("total effect:                  ", round(ab + coef(ymx)["x"], 3)))

# Bayesian version ----------------------------------------

code <- nimbleCode({
    # Model
    for (i in 1:N) {
        m[i] ~ dnorm(m_mean[i], m_prec)
        m_mean[i] <- beta3 + alpha*x[i]
        
        y[i] ~ dnorm(y_mean[i], y_prec)
        y_mean[i] <- beta1 + beta*m[i] + tau_prime*x[i]
    }
    
    # Priors
    beta ~ dnorm(0, 1/100)
    beta1 ~ dnorm(0, 1/100)
    beta3 ~ dnorm(0, 1/100)
    alpha ~ dnorm(0, 1/100)
    tau_prime ~ dnorm(0, 1/100)
    m_prec ~ dgamma(0.001, 0.001)
    y_prec ~ dgamma(0.001, 0.001)
    
    # Monitor effects
    AlphaBeta <- alpha*beta
    m_sig <- 1/sqrt(m_prec)
    y_sig <- 1/sqrt(y_prec)
})

consts <- list(x = x, N = length(y))
data <- list(m = m, y = y)
inits <- list(beta = 0, beta1 = 0, beta3 = 0, alpha = 0,
    tau_prime = 0, m_prec = 1, y_prec = 1)

nimbleMCMC_samples <- nimbleMCMC(code = code, 
    constants = consts, 
    data = data, 
    inits = inits,
    nburnin = 1000, niter = 10000, thin = 10, nchains = 3,
    monitors = c("beta", "beta1", "beta3", "alpha", "tau_prime", 
        "y_sig", "m_sig",
        "AlphaBeta"))

samples_df <- lapply(1:length(nimbleMCMC_samples), function(x) {
    df = as.data.frame(nimbleMCMC_samples[[x]])
    df$chain <- x
    df$iter <- 1:nrow(df)
    df
}) %>% bind_rows() 

samples_long <- pivot_longer(samples_df, -c("iter", "chain"))
ggplot(samples_long) +
    aes(x = iter, y = value, colour = factor(chain)) +
    geom_line() + 
    facet_wrap(~ name) +
    labs(title = "Convergence is Excellent.")

ggplot(samples_df) +
    aes(x = AlphaBeta, colour = factor(chain)) +
    geom_density() +
    geom_vline(xintercept = ab) +
    labs(title = "Indirect Effect",
        subtitle = "Black line is MLE from Frequentist analysis")

# In conclusion, Nimble does this just fine and Bayesian estimates are
# equivalent to frequentist results.

# Hooray
