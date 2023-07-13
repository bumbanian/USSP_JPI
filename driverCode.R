
#measurement
obs = list("mgca.f" = 5.3, "d18O.f" = -1.7)

post = invert(obs, 1e4)
post2 = invert(obs, 1e4)

#plots

plot(post$sst, type = "l")
lines(post2$sst, col = 2)

plot(post$d18O.sw, type = "l")
lines(post2$d18O.sw, col = 2)

plot(density(post$sst))
lines(density(post2$sst), col = 2)

plot(density(post$d18O.sw))
lines(density(post2$d18O.sw), col = 2)

plot(post$mgca.b, type = "l")
plot(post$mgca.m, type = "l")
plot(post$d18O.b, type = "l")
plot(post$d18O.m, type = "l")

#all data
d = data.frame("Depths" = c(0.890,0.940,0.990,1.040,1.090,1.140,1.190,1.200,1.210,1.220,1.230,1.240,1.250,1.260,1.270,1.280,1.290,1.300,1.31,1.32,1.33,1.34,1.39,1.44,1.48),
                  "mgca.f" = c(3.919,4.092,4.242,4.158,4.309,4.708,5.401,5.406,5.517,5.523,5.571,5.282,5.592,5.369,5.267,5.419,5.224,5.188,5.151,5.176,4.989,4.732,3.706,3.618,3.724),
                  "d18O.f" = c(-1.43,-1.54,-1.47,-1.67,-1.56,-1.65,-1.67,-1.72,-1.63,-1.68,-1.75,-1.76,-1.89,-1.99,-1.95,-1.93,-2.13,-1.98,-1.99,-1.93,-1.94,-1.9,-1.5,-1.43,-1.42))
obs.all = list("mgca.f" = d$mgca.f, "d18O.f" = d$d18O.f)

post = invert.all(obs.all, 1e4)

plot(post$sst[,1], type = "l")
plot(post$d18O.sw[,1], type = "l")
plot(post$parms$mgca.b, type = "l")
plot(post$parms$mgca.m, type = "l")
plot(post$parms$d18O.b, type = "l")
plot(post$parms$d18O.m, type = "l")
