library(nimble)
library(codatools)
library(ggplot2)
library(dplyr)
library(patchwork)

log_growth <- function(t, asymptote, midpoint, scale) {
    asymptote / (1 + exp((midpoint - t) * scale))
}

t <- seq(1, 30, 1)
p <- log_growth(t, 0.9, 10, 0.3)
n <- t^2 - 3 * t
n <- n - min(n) + 5
y <- rbinom(length(t), prob = p, size = n)
y / n


log_growth_code <- nimbleCode({
    for (i in 1:N) {
        p[i] <- asymptote / (1 + exp((midpoint - t[i]) * scale))
        Y[i] ~ dbinom(prob = p[i], size = size[i])
        Y2[i] ~ dbinom(prob = p[i], size = size[i])
    }

    # Priors
    # spike and slab so asymptote can be exactly 1
    asymptote01 ~ dbeta(6, 2) # Left skewed, [0, 1]
    asymptote1 ~ dbinom(prob = 0.75, size = 1) # 75% chance asy==1
    asymptote <- asymptote1*asymptote01 + 1 - asymptote1
    midpoint ~ dunif(0, maxt)
    scale ~ dgamma(2, 2) # [0, Inf]
})

log_growth_data <- list(Y = y, 
    Y2 = rep(NA, length(y)))
log_growth_cons <- list(
    t = t,
    N = length(y),
    size = n)
log_growth_cons$maxt <- max(t)
log_growth_init <- list(
    asymptote = c(0.5, 0.75, 0.99), 
    midpoint = mean(log_growth_cons$t) + c(-10, 0, 10),
    scale = c(0.25, 0.75, 1.5)
    )
log_growth_init <- function() {
    list(
        asymptote01 = runif(1, 0.5, 1),
        asymptote1 = rbinom(1, 1, 0.75),
        midpoint = runif(1, 
            mean(log_growth_cons$t) - 10, 
            mean(log_growth_cons$t) - 10),
        scale = runif(1, 0.25, 1.5)
    )
}

log_growth_coda <- nimbleMCMC(
    code = log_growth_code,
    constants = log_growth_cons,
    data = log_growth_data,
    inits = log_growth_init,
    nburnin = 1000, niter = 5000, nchains = 3,
    monitors = c("asymptote", "scale", "midpoint", 
        "Y2", "p"),
    samplesAsCodaMCMC = TRUE)


codf <- coda_df(log_growth_coda)

# Parameters
gasy <- ggplot(codf) +
    aes(x = iter, y = asymptote, colour = factor(chain)) +
    geom_line()
gmid <- ggplot(codf) +
    aes(x = iter, y = midpoint, colour = factor(chain)) +
    geom_line()
gsca <- ggplot(codf) +
    aes(x = iter, y = scale, colour = factor(chain)) +
    geom_line()
gasy + gmid + gsca + plot_layout(guides = "collect")
