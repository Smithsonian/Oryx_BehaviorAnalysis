model{
  
  #  PRIORS FOR GLOBAL PARAMETERS (alphas represent intercept for relative probs. of each outcome at the pre-treatment period on the log scale)
  alpha[1] <- 0  # Here we are fixing the relative probability of the reference outcome (HU) to zero on the log scale
  for (j in 2:n.outcomes) {   # loop over response outcomes
    alpha[j] ~ dnorm(0, 0.001)  # assign diffuse priors to the relative probs. of all outcomes except the reference (for period 1)
  }
  
  # betas represent change from period one on the log scale
  for (j in 1:n.outcomes){   # loop over response outcomes
    beta[1, j] <- 0   # beta[1,] are fixed to zero because there is no period-adjustment needed for period 1 since it is the reference
  }
  
  for (i in 2:3) {    # loop over time periods 2 and 3 (post-collaring)
    beta[i, 1] <- 0   # as for period 1, we have to fix the relative prob. of the reference outcome (HU) to zero on the log scale
    for (j in 2:n.outcomes){   # loop over response outcomes
      beta[i, j] ~ dnorm(0, 0.001)  # assign diffuse priors to change (periods 1-2 and periods 1-3) in rel. probs. of outcomes
    } 
  }
  
  # tau.j parameters represent inter-individual variation in relative probs. of outcomes
  tau.j[1] <- 1   # this is a placeholder zero - we always fix the rel. prob. of the reference outcome to zero so there is no inter-indiv. variation for j=1
  for (j in 2:n.outcomes){  # loop over response outcomes
    tau.j[j] ~ dgamma(1,0.001)  # diffuse gamma prior for precisions (1/variance)
  }
  
  #  DEFINE INDIVIDUAL-LEVEL PARAMETERS
  for (idx in 1:nind){   # Loop over individuals to define individual-level random effects
    eps[idx,1] <- 0      # Rel. prob. of reference outcome fixed to zero, so there is no adjustment among indiviuals
    for (j in 2:n.outcomes){  # loop over outcomes
      eps[idx,j] ~ dnorm(0,tau.j[j])  # Adjustment to rel. probs. of each response outcome to account for individual variation
    }
  }
  
  # LIKELIHOOD  
  for (i in 1:n) {     # loop over observations
    # Multinomial response
    Y[i, ] ~ dmulti(p[i, ] , N[i])

    for (j in 1:n.outcomes) {     # loop around
      p[i,j] <- phi[i,j] / sum(phi[i, ])
      log(phi[i,j]) <- alpha[j] + beta[PERIOD[i], j] + eps[ind[i], j]
    }
  }
  
  # Derived quantities
  for (j in 1:n.outcomes){
    PROBS[1,j] <- PHI[1,j] / sum(PHI[1,])
    log(PHI[1,j]) <- alpha[j] + beta[1,j]
    PROBS[2,j] <- PHI[2,j] / sum(PHI[2,])
    log(PHI[2,j]) <- alpha[j] + beta[2,j]
    PROBS[3,j] <- PHI[3,j] / sum(PHI[3,])
    log(PHI[3,j]) <- alpha[j] + beta[3,j]
    sigma2[j] <- 1/tau.j[j] 
  }
}
