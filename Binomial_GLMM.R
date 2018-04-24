# JAGS model

model
{
    # Priors
    for (i in 1:n.groups){ # n.groups is the number of animals
      alpha[i] ~ dnorm(mu.int, tau.int)  # Intercepts
      #beta[i] ~ dnorm(mu.beta, tau.beta) # Slopes
    }
    
    mu.int ~ dnorm(0, 0.001) # Hyperparameter for random intercepts
    tau.int <- 1 / (sigma.int * sigma.int)
    sigma.int ~ dunif(0,10)
    
    beta1 ~ dnorm(0,0.001) # Individual group beta (post collaring)
    beta2 ~ dnorm(0,0.001) # Individual group beta (post collaring + 3 days)
    
    # Beta for time past 7
    beta3 ~ dnorm(0,0.001)
    
    #mu.beta ~ dnorm(0, 0.001)
    #tau.beta <- 1 / (sigma.beta * sigma.beta)
    #sigma.beta ~ dunif(0,10)
  
    # Binomial Likelihood
    for (i in 1:n){
      C[i] ~ dbin(p[i], N[i])
      logit(p[i]) <- alpha[ID[i]] + beta1 * equals(PERIOD[i],2) + beta2 * equals(PERIOD[i],3) + beta3 * TimePast[i]
    }
    
    # Derived parameters
    Contrast1v2 <- beta1
    Contrast1v3 <- beta2
    Contrast2v3 <- beta2-beta1  # Group comparison
}