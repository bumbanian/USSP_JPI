
#initialize the chain
init = function(){
  sst = runif(1, 15, 35)
  d18O.sw = rnorm(1, 0.5, 0.5)
  mgca.b = rnorm(1, 0.4, 0.02)
  mgca.m = rnorm(1, 0.096, 0.004)
  d18O.b = rnorm(1, 2.9, 0.075)
  d18O.m = rnorm(1, 0.213, 0.005)
  
  return(list("sst" = sst, "d18O.sw" = d18O.sw,
              "mgca.b" = mgca.b, "mgca.m" = mgca.m,
              "d18O.b" = d18O.b, "d18O.m" = d18O.m))
}

#proposal
prop = function(parm){
  sst = rnorm(1, parm$sst, 0.4)
  d18O.sw = rnorm(1, parm$d18O.sw, 0.08)
  mgca.b = rnorm(1, parm$mgca.b, 0.02)
  mgca.m = rnorm(1, parm$mgca.m, 0.0012)
  d18O.b = rnorm(1, parm$d18O.b, 0.03)
  d18O.m = rnorm(1, parm$d18O.m, 0.002)
  return(list("sst" = sst, "d18O.sw" = d18O.sw,
              "mgca.b" = mgca.b, "mgca.m" = mgca.m,
              "d18O.b" = d18O.b, "d18O.m" = d18O.m))
}

#likelihood
like = function(mod, obs, parm){
  small = 1e-9
  #obs
  mgca.f.sd = 0.1
  d18O.f.sd = 0.15
  mgca.f.p = log(max(dnorm(obs$mgca.f, mod$mgca.f, mgca.f.sd), small))
  d18O.f.p = log(max(dnorm(obs$d18O.f, mod$d18O.f, d18O.f.sd), small))
  
  #parm
  sst.p = log(max(dunif(parm$sst, 15, 35), small))
  d18O.sw.p = log(max(dnorm(parm$d18O.sw, 0.5, 0.5), small))
  mgca.b.p = log(max(dnorm(parm$mgca.b, 0.4, 0.02), small))
  mgca.m.p = log(max(dnorm(parm$mgca.m, 0.096, 0.004), small))
  d18O.b.p = log(max(dnorm(parm$d18O.b, 2.9, 0.1), small))
  d18O.m.p = log(max(dnorm(parm$d18O.m, 0.213, 0.005), small))
  
  return(sum(c(mgca.f.p, d18O.f.p, sst.p, d18O.sw.p, mgca.b.p,
               mgca.m.p, d18O.b.p, d18O.m.p)))
}

#PSM calcs
psm = function(parm){
  mgca.f = parm$mgca.b * exp(parm$mgca.m * parm$sst)
  d18O.f = parm$d18O.sw - parm$sst * parm$d18O.m + parm$d18O.b
  return(list("mgca.f" = mgca.f, "d18O.f" = d18O.f))
}

#run inversion
invert = function(obs, n = 1000, burnin = n/5){
  
  #space for results
  parms = data.frame("sst" = double(n), "d18O.sw" = double(n), 
                     "mgca.b" = double(n), "mgca.m" = double(n),
                     "d18O.b" = double(n), "d18O.m" = double(n))
  mods = data.frame("mgca.f" = double(n), "d18O.f" = double(n))
  
  #start with initial condition
  parm = init()
  mod = psm(parm)
  l = like(mod, obs, parm)
  
  #counters
  i = 1; ti = 1
  
  #progress bar
  pb = txtProgressBar(0, n)
  
  #iterate until we have the desired number of samples
  while(i <= n){
    #proposal
    parm.n = prop(parm)
    mod.n = psm(parm.n)
    l.n = like(mod.n, obs, parm.n)
    
    #Hastings ratio
    H = exp(l.n - l)
    
    #accept?
    if(runif(1) < H){
      #save values
      parms[i,] = (unlist(parm.n))
      mods[i,] = unlist(mod)
      
      #update current state
      parm = parm.n
      mod = mod.n
      l = l.n
      i = i + 1
      setTxtProgressBar(pb, i)
    }
    ti = ti + 1
  }
  
  #acceptance rate
  cat("\nAcceptance rate:", signif(n / ti, 2), "\n")
  
  #package results and trim the burnin period
  result = cbind(parms, mods)
  result = result[-c(1:burnin),]
  
  return(result)
}

#Full dataset####

#initialize the chain
init.c = function(){
  mgca.b = rnorm(1, 0.4, 0.02)
  mgca.m = rnorm(1, 0.096, 0.004)
  d18O.b = rnorm(1, 2.9, 0.1)
  d18O.m = rnorm(1, 0.213, 0.005)
  
  return(list("mgca.b" = mgca.b, "mgca.m" = mgca.m,
              "d18O.b" = d18O.b, "d18O.m" = d18O.m))
}

init.v = function(){
  sst = runif(1, 28, 30)
  d18O.sw = runif(1, 0, 0.5)

  return(list("sst" = sst, "d18O.sw" = d18O.sw))
}

#proposal
prop.c = function(parm){
  mgca.b = rnorm(1, parm$mgca.b, 0.003)
  mgca.m = rnorm(1, parm$mgca.m, 0.0004)
  d18O.b = rnorm(1, parm$d18O.b, 0.02)
  d18O.m = rnorm(1, parm$d18O.m, 0.0007)
  return(list("mgca.b" = mgca.b, "mgca.m" = mgca.m,
              "d18O.b" = d18O.b, "d18O.m" = d18O.m))
}

prop.v = function(parm){
  sst = rnorm(1, parm$sst, 0.1)
  d18O.sw = rnorm(1, parm$d18O.sw, 0.02)
  return(list("sst" = sst, "d18O.sw" = d18O.sw))
}

#likelihood
like.c = function(parm){
  small = 1e-15
  #parm
  mgca.b.p = log(max(dnorm(parm$mgca.b, 0.4, 0.02), small))
  mgca.m.p = log(max(dnorm(parm$mgca.m, 0.096, 0.002), small))
  d18O.b.p = log(max(dnorm(parm$d18O.b, 2.9, 0.1), small))
  d18O.m.p = log(max(dnorm(parm$d18O.m, 0.213, 0.005), small))
  
  return(sum(c(mgca.b.p, mgca.m.p, d18O.b.p, d18O.m.p)))
}

like.v = function(mod, obs, parm){
  small = 1e-15
  #obs
  mgca.f.sd = 0.15
  d18O.f.sd = 0.15
  mgca.f.p = log(max(dnorm(obs$mgca.f, mod$mgca.f, mgca.f.sd), small))
  d18O.f.p = log(max(dnorm(obs$d18O.f, mod$d18O.f, d18O.f.sd), small))
  
  #parm
  sst.p = log(max(dunif(parm$sst, 15, 35), small))
  d18O.sw.p = log(max(dnorm(parm$d18O.sw, 0.5, 0.5), small))

  return(sum(c(mgca.f.p, d18O.f.p, sst.p, d18O.sw.p)))
}

#run inversion
invert.all = function(obs, n = 1000, burnin = n/5){
  
  #number of obs
  nobs = length(obs$mgca.f)
  
  #space for results
  parms = data.frame("mgca.b" = double(n), "mgca.m" = double(n),
                     "d18O.b" = double(n), "d18O.m" = double(n))
  sst = d18O.sw = matrix(nrow = n, ncol = nobs)

  #space for initial condition
  parm.c = parm.v = mod = parm.c.n = parm.v.n = mod.n = list()
  l.v = l.v.n = double(nobs)
  
  #get initial conditions
  parm.c = init.c()
  l.c = like.c(parm.c)
  for(j in seq(nobs)){
    parm.v[[j]] = init.v()
    mod[[j]] = psm(append(parm.v[[j]], parm.c))
    l.v[j] = like.v(mod[[j]], 
             list("mgca.f" = obs$mgca.f[j], "d18O.f" = obs$d18O.f[j]), 
             parm.v[[j]])
  }
  
  #counters
  i = 1; ti = 1
  
  #progress bar
  pb = txtProgressBar(0, n)
  
  #iterate until we have the desired number of samples
  while(i <= n){
    #fixed parms
    parm.c.n = prop.c(parm.c)
    l.c.n = like.c(parm.c.n)
    
    #cycle through all obs
    for(j in seq(nobs)){
      parm.v.n[[j]] = prop.v(parm.v[[j]])
      mod.n[[j]] = psm(append(parm.v.n[[j]], parm.c.n))
      l.v.n[j] = like.v(mod.n[[j]], 
                    list("mgca.f" = obs$mgca.f[j], "d18O.f" = obs$d18O.f[j]), 
                    parm.v.n[[j]])
    }
    
    #Hastings ratio
    H = exp(sum(c(l.c.n, l.v.n)) - sum(c(l.c, l.v)))
    
    #accept?
    if(runif(1) < H){
      #save values
      parms[i,] = unlist(parm.c.n)
      for(j in seq(nobs)){
        sst[i, j] = parm.v.n[[j]]$sst
        d18O.sw[i, j] = parm.v.n[[j]]$d18O.sw
      }

      #update current state
      parm.c = parm.c.n
      parm.v = parm.v.n
      mod = mod.n
      l.c = l.c.n
      l.v = l.v.n
      i = i + 1
      setTxtProgressBar(pb, i)
    }
    ti = ti + 1
  }
  
  #acceptance rate
  cat("\nAcceptance rate:", signif(n / ti, 2), "\n")
  
  #package results and trim the burnin period
  sst = sst[-c(1:burnin),]
  d18O.sw = d18O.sw[-c(1:burnin),]
  parms = parms[-c(1:burnin),]
  result = list("sst" = sst, "d18O.sw" = d18O.sw, "parms" = parms)

  return(result)
}
