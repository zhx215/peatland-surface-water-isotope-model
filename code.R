rm(list=ls())

height <- 10 # m
distance <- 1000 # m
dHdL_ratio <- height/distance # topographic gradient
K <- 900/86400 # m/s, hydraulic conductivity, 900 m/day
v <- K*dHdL_ratio # m/s, flow velocity
dx <- 1 # m, distance step
dt <- dx/v # s, time step
PET <- 112/1000 # m, potential ET
AETPET_ratio <- 0.6 # actual ET/potential ET ratio
ET <- PET*AETPET_ratio/(2.592e+6) # m/s, actual ET
S <- 0.81*30*10/1000 # m, thickness of surface water, wet bulk density 0.9 g/cm3, 90%, 30 cm
ET_S_ratio <- ET/S # /s, turnover timescale
f <- 1-ET_S_ratio*dt # residual liquid fraction
nd <- 1/f-1 # Damkohler number
tet <- 0.75 # transpiration fraction
rh <- 0.85 # relative humidity
d18O_G <- -12 # permil, groundwater
d2H_G <- 8*d18O_G+10 # permil, groundwater
d18O_S0 <- -12 # permil, initial surface water
d2H_S0 <- 8*d18O_S0+10 # permil, initial surface water
a18_e_vl <- function(temp){
  1/exp(1137/temp^2-0.4156/temp-0.00207)
}
a2_e_vl <- function(temp){
  1/exp(24844/temp^2-76.248/temp+0.05261)
}
a18 <- a18_e_vl(273.15+18) # 18 degree C, temperature
a2 <- a2_e_vl(273.15+18) # 18 degree C, temperature
k18 <- 0.9723^0.5 # wet soil, aerodynamic exponent
k2 <- 0.9755^0.5 # wet soil, aerodynamic exponent
smow18 <- 0.0020052
smow2 <- 155.76*10^(-6)
rayleigh_ET <- function(nd,t,h,delta_initial,alpha,k,smow){ # isotopic composition of ET flux, see Xia and Winnick, 2021
  frac <-  1/(1/nd+1)
  r_initial <- (delta_initial/1000+1)*smow
  ae <- ((1-t)*((alpha*k*r_initial)/(1-h+(1-t)*k*h))+t*r_initial*(1/(1+(1-t)*k*(h/(1-h)))))/r_initial
  delta_ET <- ((delta_initial+1000)*(1-frac)^ae-1000*(1-frac)-delta_initial)/(1-frac-1)*(1-t)+delta_initial*t
  return(c((delta_initial-frac*delta_ET)/(1-frac),delta_ET))
}
d18O_S <- d18O_S0 # surface water vector
d2H_S <- d2H_S0 # surface water vector
x <- seq(1,distance,dx) # distance vector
for (i in 1:(length(x)-1)){
  d18O_S[i+1] <- d18O_S[i]+dx*ET_S_ratio/v*(d18O_G-rayleigh_ET(nd,tet,rh,d18O_S[i],a18,k18,smow18)[2])
  d2H_S[i+1] <- d2H_S[i]+dx*ET_S_ratio/v*(d2H_G-rayleigh_ET(nd,tet,rh,d2H_S[i],a2,k2,smow2)[2])
}

plot(x,d18O_S)