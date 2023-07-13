model{
  for(i in 1:length(mgca.d)){
    mgca.d[i] ~ dnorm(mgca.f[i], 1 / 0.15^2)
    d18O.d[i] ~ dnorm(d18O.f[i], 1 / 0.15^2)
  }
  
  for(i in 1:length(mgca.d)){
    mgca.f[i] = mgca.b * exp(mgca.m * sst[i])
    d18O.f[i] = d18O.sw[i] - sst[i] * d18O.m + d18O.b
    
    sst[i] ~ dunif(15, 35)
    d18O.sw[i] ~ dnorm(0.5, 1/0.5^2)
  }
  
  mgca.b ~ dnorm(0.4, 1/0.02^2)
  mgca.m ~ dnorm(0.96, 1/0.004^2)
  d18O.b ~ dnorm(2.9, 1/0.1^2)
  d18O.m ~ dnorm(0.213, 1/0.005^2)
}