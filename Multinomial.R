model{
  
  #  PRIORS
  alpha[1] <- 0;       # zero contrast for baseline activity
  for (j in 2:n.outcomes) {   # loop over activities
    alpha[j] ~ dnorm(0, 0.001)
  }
  
  for (j in 1:n.outcomes){  
    beta[1, j] <- 0   # zero contrast for control treatment
  }
  for (i in 2:3) {    # loop over treatments
    beta[i, 1] <- 0   # zero contrast for baseline activity
    for (j in 2:n.outcomes){  
      beta[i, j] ~ dnorm(0, 0.001)
    } 
  }
  
  # LIKELIHOOD  
  for (i in 1:n) {     # loop over observations
    # Multinomial response
    Y[i, ] ~ dmulti(p[i, ] , N[i])

    for (j in 1:n.outcomes) {     # loop around
      p[i,j] <- phi[i,j] / sum(phi[i, ])
      log(phi[i,j]) <- alpha[j] + beta[PERIOD[i], j]
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
  }
}
