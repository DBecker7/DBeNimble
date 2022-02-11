# D. Be Nimble

"DBeNimble" is short for "Devan BEcker Nimble", in case that wasn't obvious.

My collection of [NIMBLE](https://r-nimble.org/) examples/snippets for personal reference.

The official documentation is detailed and complete, but I just want the minimal examples.

# D. Be Quick

[From the BUGS example.](https://r-nimble.org/nimbleExamples/nimble_build_a_model.html)

```r
library(nimble, warn.conflicts = FALSE)

pumpCode <- nimbleCode({ 
# Define relationships between nodes
    for (i in 1:N){
        theta[i] ~ dgamma(alpha,beta)
        lambda[i] <- theta[i]*t[i]
        x[i] ~ dpois(lambda[i])
    }
    
    # Set priors
    alpha ~ dexp(1.0)
    beta ~ dgamma(0.1,1.0)
})

# Create some contrants, data, and initial values to pass to the model builder
pumpConsts <- list(N = 10,
    t = c(94.3, 15.7, 62.9, 126, 5.24,
        31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
    theta = rep(0.1, pumpConsts$N))

nimbleMCMC_samples <- nimbleMCMC(code = pumpCode, 
    constants = pumpConsts, 
    data = pumpData, 
    inits = pumpInits,
    nburnin = 1000, niter = 10000)
```

