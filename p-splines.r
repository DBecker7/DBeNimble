library(nimble)
library(splines)

plot(mpg ~ hp, data = mtcars)

# https://bragqut.github.io/2016/05/24/samclifford-splines/
bspline <- function(x, K, bdeg = 3, cyclic = FALSE,
                    xl = min(x), xr = max(x)) {
    x <- as.matrix(x, ncol = 1)
    ndx <- K - bdeg

    # as outlined in Eilers and Marx (1996)
    dx <- (xr - xl) / ndx
    t <- xl + dx * (-bdeg:(ndx + bdeg))
    T <- (0 * x + 1) %*% t
    X <- x %*% (0 * t + 1)
    P <- (X - T) / dx
    B <- (T <= X) & (X < (T + dx))
    r <- c(2:length(t), 1)

    for (k in 1:bdeg){
        B <- (P * B + (k + 1 - P) * B[, r]) / k
    }

    B <- B[, 1:(ndx + bdeg)]

    if (cyclic == 1) {
        for (i in 1:bdeg){
        B[, i] <- B[, i] + B[, K - bdeg + i]
        }
        B <- B[, 1:(K - bdeg)]
    }
    
    return(B)
}
x <- seq(0, 400, 0.1)
xb <- bspline(x, 10)
matplot(xb, type = "l")
b <- bspline(mtcars$hp, 10)

a <- solve(t(b) %*% b) %*% t(b) %*% mtcars$mpg
newmpg <- b %*% a
plot(x = mtcars$hp, y = mtcars$mpg)
lines(x = sort(mtcars$hp), y = newmpg[order(mtcars$hp)])
# Too wiggly! Let's penalise.

d1_model <- nimbleCode({
    for (i in 1:N) {
        mu[i] <- inprod(beta[1:K, i], alpha[1:K])
        y[i] ~ dnorm(mean = mu[i],
            sd = sigma_y)
    }

    alpha[1] ~ dnorm(0, sd = sigma_a)
    for (k in 2:K) {
        alpha[k] ~ dnorm(alpha[k - 1], sd = sigma_a)
    }

    sigma_y ~ dgamma(0.1, 0.1)
    sigma_a ~ dgamma(0.1, 0.1)
})

d1_const <- list(beta = t(b), N = nrow(mtcars), K = ncol(b))
d1_data <- list(y = mtcars$mpg)

d1_samples <- nimbleMCMC(
    code = d1_model, const = d1_const,
    data = d1_data,
    niter = 5000,
    monitors = c("alpha")
)
alpha_est <- apply(d1_samples, 2, median)
d1_est <- b %*% alpha_est
lines(sort(mtcars$hp), y = d1_est[order(mtcars$hp)],
    col = 2)

d2_model <- nimbleCode({
    for (i in 1:N) {
        mu[i] <- inprod(beta[1:K, i], alpha[1:K])
        y[i] ~ dnorm(mean = mu[i],
            sd = sigma_y)
    }

    alpha[1] ~ dnorm(0, sd = sigma_a)
    alpha[2] ~ dnorm(alpha[1], sd = sigma_a)
    for (k in 3:K) {
        alpha[k] ~ dnorm(
            mean = -2 * alpha[k - 1] + alpha[k - 2],
            sd = sigma_a)
    }

    sigma_y ~ dgamma(0.1, 0.1)
    sigma_a ~ dgamma(0.1, 0.1)
})

d2_const <- list(beta = t(b), N = nrow(mtcars), K = ncol(b))
d2_data <- list(y = mtcars$mpg)

d2_samples <- nimbleMCMC(
    code = d2_model, const = d2_const,
    data = d2_data,
    niter = 5000,
    monitors = c("alpha")
)
alpha_est <- apply(d2_samples, 2, median)
d2_est <- b %*% alpha_est
lines(sort(mtcars$hp), y = d2_est[order(mtcars$hp)],
    col = 3)


qvar_model <- nimbleCode({
    for (i in 1:N) {
        mu[i] <- inprod(beta[1:K, i], alpha[1:K])
        y[i] ~ dnorm(mean = mu[i], sd = sigma_y)
    }

    alpha[1:K] ~ dmnorm(zeros[1:K], prec[1:K, 1:K])

    sigma_y ~ dgamma(0.1, 0.1)
    lambda ~ dgamma(0.1, 0.1)

    prec[1:K, 1:K] <- lambda * Q[1:K, 1:K]
})

make_q <- function(degree, k, epsilon = 1e-3) {
  x <- diag(k)
  e <- diff(x, differences = degree)
  return(t(e) %*% e + x * epsilon)
}

q <- make_q(2, ncol(b))

qvar_const <- list(beta = t(b), N = nrow(mtcars), K = ncol(b),
    Q = q, zeros = rep(0, ncol(b)))
qvar_data <- list(y = mtcars$mpg)

qvar_samples <- nimbleMCMC(
    code = qvar_model, const = qvar_const,
    data = qvar_data,
    niter = 15000,
    monitors = c("alpha")
)
alpha_est <- apply(qvar_samples, 2, median)
q_est <- b %*% alpha_est
lines(sort(mtcars$hp), y = q_est[order(mtcars$hp)],
    col = 4)
legend("topright", legend = c("None", "1", "2", "Q"), col = 1:4, lty = 1)
