# D. Be Nimble

"DBeNimble" is short for "Devan BEcker Nimble", in case that wasn't obvious.

My collection of [NIMBLE](https://r-nimble.org/) examples/snippets for personal reference.

The official documentation is detailed and complete, but I just want the minimal examples.

# D. Be Quick

[From the BUGS example,](https://r-nimble.org/nimbleExamples/nimble_build_a_model.html) which is based on [this model](https://www.multibugs.org/examples/latest/Pumps.html).

Number of failures in a power plant pump is modelled as theta times the length of operating time. 

```r
library(nimble, warn.conflicts = FALSE)

pumpCode <- nimbleCode({ 
    for (i in 1:N){
        theta[i] ~ dgamma(alpha, beta)
        lambda[i] <- theta[i]*t[i]
        x[i] ~ dpois(lambda[i])
    }
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1, 1.0)
})

# Note that covars are supplied as constants
pumpConsts <- list(N = 10,
    t = c(94.3, 15.7, 62.9, 126, 5.24,
        31.4, 1.05, 1.05, 2.1, 10.5))

# data refers to response
pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
    theta = rep(0.1, pumpConsts$N))
# Could also specify inits as a function for arbitrary number of chains
pumpInitsFun <- function() {
    alpha = runif(1, 0.5, 1.5),
    beta = runif(1, 0.5, 1.5),
    theta = runif(pumpConsts$N, 0, 0.2)
}

nimbleMCMC_samples <- nimbleMCMC(code = pumpCode, 
    constants = pumpConsts, 
    data = pumpData, 
    inits = pumpInits, # inits = pumpInitsFun
    nburnin = 1000, niter = 10000,
    monitors = c("alpha", "beta", "theta"),
    # monitors = names(pumpInitsFun()),
    samplesAsCodaMCMC = TRUE)

samples_df <- lapply(1:length(nimbleMCMC_samples), function(x) {
    df = as.data.frame(nimbleMCMC_samples[[x]])
    df$chain <- x
    df$iter <- 1:nrow(df)
    df
}) %>% bind_rows() 

ggplot(nimbleMCMC_samples) +
    aes(x = iter, y = value, colour = factor(chain)) +
    geom_line() + 
    facet_wrap(~ name)
```


# D. Be Confused

- `X[1:n, 1:p]` and `beta[1:p]`
    - must specify the dimensions and length
- `data` is the target, `constants` are the features (incl., *e.g.*, `n = length(response)`)
- Matrix multiplication retains matrix class, must extract into vector.
    - `(beta[1:p] %*% x[i, 1:p])[1,1]`
- `tau2 ~ dgamma(1,1)`; `Y[i] ~ dnorm(0, sd = 1/sqrt(tau2))`
