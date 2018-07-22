model{
  
  #  PRIORS
  for (i in 1:3){
    # Fix the relative probability of the baseline behavior to zero for all time periods (i.e., pre/post collar)
    alpha[i,1] <- 0;   
    
    # Loop over response behaviors (HU, HD, Laying, Headshake, etc...)
    for (j in 2:n.outcomes) { 
      # Assign diffuse normal priors to the relative probabilities of all behaviors except the reference
      alpha[i,j] ~ dnorm(0, 0.001)
    }
  }
  
  # LIKELIHOOD 
  # *************************
  # Loop over all observations
  for (i in 1:n) {     
    # Create multinomial response
    # Per activity, counts (Y) are distributed multi-nomially with per activity (vector of probabilities - p) and N number of trials
    Y[i, ] ~ dmulti(p[i, ] , N[i])

    # Internal loop through all the outcomes/behaviors at each observation
    for (j in 1:n.outcomes) {
      # Convert relative probabilities into actual probabilities (i.e., divide by the sum of all outcomes)
      p[i,j] <- phi[i,j] / sum(phi[i, ])
      # Log link function
      log(phi[i,j]) <- alpha[period[i],j]
    }
  }
  
  # Derived quantities
  for (j in 1:n.outcomes){
    PROBS[1,j] <- PHI[1,j] / sum(PHI[1,])
    log(PHI[1,j]) <- alpha[1,j]
    PROBS[2,j] <- PHI[2,j] / sum(PHI[2,])
    log(PHI[2,j]) <- alpha[2,j]
    PROBS[3,j] <- PHI[3,j] / sum(PHI[3,])
    log(PHI[3,j]) <- alpha[3,j]
  }
}
