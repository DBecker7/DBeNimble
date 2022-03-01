library(nimble)
library(MASS)
rho <- 0.2
sigma <- matrix(data = c(1, rho, rho, 2), ncol = 2)

xy <- mvrnorm(1000, c(0, 4), Sigma = sigma)

nimble_data <- list(xy = xy)
nimble_const <- list(N = nrow(xy))
nimble_inits <- list(rho = -0.2, mu1 = 4, mu2 = 0, vx = 2, vy = 1)
nimble_model <- nimbleCode({
    for (i in 1:N) {
        xy[i, 1:2] ~ dmnorm(mu[1:2], cov = Cov[1:2, 1:2])
    }

    Cov[1, 1] <- vx
    Cov[1, 2] <- rho
    Cov[2, 1] <- rho
    Cov[2, 2] <- vy
    vx ~ dgamma(1, 1)
    vy ~ dgamma(1, 1)
    rho ~ dnorm(0, 1 / 10)
    mu[1:2] <- c(mu1, mu2)
    mu1 ~ dnorm(0, 1 / 10)
    mu2 ~ dnorm(0, 1 / 10)
})

nimble_samples <- nimbleMCMC(
    code = nimble_model,
    data = nimble_data,
    constants = nimble_const,
    inits = nimble_inits,
    nburnin = 1000, niter = 2000, nchains = 2,
    samplesAsCodaMCMC = TRUE,
    monitors = c("mu1", "mu2", "vx", "vy", "rho")
)

summary(nimble_samples)
