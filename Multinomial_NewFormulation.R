model{
  
  #  PRIORS
  for (i in 1:3){
    alpha[i,1] <- 0;   # Here we are fixing the relative probability of the baseline behavior to zero for all time periods (pre/post collar)

    for (j in 2:n.outcomes) {   # loop over response behaviors (HU, HD, Laying, Headshake, etc...)
      alpha[i,j] ~ dnorm(0, 0.001)   # Here we are assigning diffuse normal priors to the relative probabilities of all behaviors except the reference
    }
  }
  
  # LIKELIHOOD  
  for (i in 1:n) {     # loop over observations
    # Multinomial response
    Y[i, ] ~ dmulti(p[i, ] , N[i])

    for (j in 1:n.outcomes) {     # loop around
      p[i,j] <- phi[i,j] / sum(phi[i, ])
      log(phi[i,j]) <- alpha[period[i],j]
    }
  }
  
  # Derived quantities
  PROBS <- p   # We actually don't need to re-calculate these, we could just add p to the parameters to monitor or just save this as 'PROBS'
}
