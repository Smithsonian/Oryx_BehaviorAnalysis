model{
  
  # PRIORS FOR GLOBAL PARAMETERS
  # Alphas represent intercept for relative probabilities of each outcome at the pre-treatment period on the log scale
  # Fix the relative probability of the reference outcome (HU) to zero on the log scale
  alpha[1] <- 0  
  
  # Loop over response outcomes
  for (j in 2:n.outcomes) {   
    # Assign diffuse priors to the relative probabilities of all outcomes except the reference (for period 1)
    alpha[j] ~ dnorm(0, 0.001)  
  }
  
  # Betas represent change from period one on the log scale
  # Loop over response outcomes
  for (j in 1:n.outcomes){
    # beta[1,] are fixed to zero because there is no period-adjustment needed for period 1 since it is the reference
    beta[1, j] <- 0   
  }
  
  # Loop over time periods 2 (treatment) and 3 (post-treatment)
  for (i in 2:3) { 
    # As for period 1, we have to fix the relative probabilities of the reference outcome (HU) to zero on the log scale
    beta[i, 1] <- 0
    # Loop over response outcomes
    for (j in 2:n.outcomes){
      # Assign diffuse priors to change (periods 1-2 and periods 1-3) in rel. probs. of outcomes
      beta[i, j] ~ dnorm(0, 0.001)   
      } 
  }
  
  # tau.j parameters represent inter-individual variation in relative probs. of outcomes
  # Loop over response outcomes
  for (j in 1:n.outcomes){  
    # Mean of individual random effects is 0
    mu.re[j] <- 0           
  }
  
  # PRIORS FOR ELEMENTS OF PRECISION MATRIX
  # df set to j+1
  prec[1:6,1:6] ~ dwish(R[,],7)  
  # Convert precision to covariance matrix
  sigma[1:6,1:6] <- inverse(prec[,])   
  # Correlation between outcome 1 and 2
  # rho <- sigma[1,2]/sqrt(sigma[1,1]*sigma[2,2])  
  
  # DEFINE INDIVIDUAL-LEVEL PARAMETERS
  # Loop over individuals to define individual-level random effects
  for (idx in 1:nind){  
    # Rel. prob. of reference outcome fixed to zero, so there is no adjustment among indiviuals
    eps[idx,1:6] ~ dmnorm(mu.re[], prec[,])      
  }
  
  # LIKELIHOOD 
  # Loop over observations
  for (i in 1:n) {     
    # Multinomial response
    Y[i, ] ~ dmulti(p[i, ] , N[i])

    # Loop through outcomes
    for (j in 1:n.outcomes) {     
      p[i,j] <- phi[i,j] / sum(phi[i, ])
      log(phi[i,j]) <- alpha[j] + beta[PERIOD[i], j] + eps[ind[i], j]
    }
  }
  
  # DERIVED QUANTITIES
  for (j in 1:n.outcomes){
    PROBS[1,j] <- PHI[1,j] / sum(PHI[1,])
    log(PHI[1,j]) <- alpha[j] + beta[1,j]
    PROBS[2,j] <- PHI[2,j] / sum(PHI[2,])
    log(PHI[2,j]) <- alpha[j] + beta[2,j]
    PROBS[3,j] <- PHI[3,j] / sum(PHI[3,])
    log(PHI[3,j]) <- alpha[j] + beta[3,j]
  }
}
