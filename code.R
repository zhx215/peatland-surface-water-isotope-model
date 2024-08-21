rm(list=ls())

# ΔH, vertical gradient of the transect
height <- 10 # m

# ΔL, horizontal gradient of the transect
distance <- 1000 # m

# ΔH/ΔL, topographic gradient
dHdL_ratio <- height/distance # dimensionless

# saturated hydraulic conductivity
Ks <- 900/86400 # m/s, original value is 900 m/day (converted to m/s)

# horizontal flow velocity
v <- Ks*dHdL_ratio # m/s

# distance in a numerical step
dx <- 1 # m

# time step
dt <- dx/v # s

# potential ET, derived from the Thornthwaite method (see http://ponce.sdsu.edu/onlinethornthwaite.php)
PET <- 112/1000 # m/month, original value is 112 mm/month (converted to m) 

# actual ET/potential ET ratio, estimated from literature data
AETPET_ratio <- 0.6 # dimensionless

# actual ET rate
ET <- PET*AETPET_ratio/(2.592e+6) # m/s (converted from m/month to m/s)

# thickness of surface water, given wet bulk density 0.9 g/cm3, 90% water content, and 30 cm depth
S <- 0.9*0.9*30/100 # m (converted from mm to m)

# ET/S, relative evapotranspiration rate (water turnover timescale)
ET_S_ratio <- ET/S # /s

# residual liquid fraction, used for δET calculation (see Xia and Winnick, 2021; https://doi.org/10.1016/j.epsl.2021.117120)
f <- 1-ET_S_ratio*dt # dimensionless

# Damkohler number (see Xia et al., 2021; https://doi.org/10.1016/j.epsl.2021.117120)
nd <- 1/f-1 # dimensionless

# T/ET, transpiration/evapotranspiration ratio
tet <- 0.75 # dimensionless

# relative humidity
rh <- 0.85 # dimensionless

# δG, the oxygen and hydrogen isotopic compositions of groundwater flux
d18O_G <- -12 # permil
d2H_G <- 8*d18O_G+10 # permil

# δS(x=0), the initial oxygen and hydrogen isotopic compositions of surface water
d18O_S0 <- -12 # permil
d2H_S0 <- 8*d18O_S0+10 # permil

# vapor-liquidequilibrium fractionation equation (see Xia and Winnick, 2021; https://doi.org/10.1016/j.epsl.2021.117120)
a18_e_vl <- function(temp){
  1/exp(1137/temp^2-0.4156/temp-0.00207) 
}
a2_e_vl <- function(temp){
  1/exp(24844/temp^2-76.248/temp+0.05261) 
}

# vapor-liquidequilibrium fractionation factor, assuming air temperature of 18 degree C
a18 <- a18_e_vl(273.15+18) # dimensionless
a2 <- a2_e_vl(273.15+18) # dimensionless

# kinetic fractionation factor, assuming wet soil (aerodynamic exponent=0.5)
k18 <- 0.9723^0.5 # dimensionless
k2 <- 0.9755^0.5 # dimensionless

# isotopic ratios of Vienna Standard Mean Ocean Water (VSMOW)
smow18 <- 0.0020052
smow2 <- 155.76*10^(-6)

# equation to derive the isotopic composition of evapotranspiration flux (see Xia and Winnick, 2021; https://doi.org/10.1016/j.epsl.2021.117120)
rayleigh_ET <- function(nd,t,h,delta_initial,alpha,k,smow){ # inputs as follows: Damkohler number, transpiration/evapotranspiration ratio, relative humidity, the isotopic composition of surface water, equilibrium fractionation factor, kinetic fractionation factor, isotopic ratio of VSMOW
  frac <-  1/(1/nd+1)
  r_initial <- (delta_initial/1000+1)*smow
  ae <- ((1-t)*((alpha*k*r_initial)/(1-h+(1-t)*k*h))+t*r_initial*(1/(1+(1-t)*k*(h/(1-h)))))/r_initial
  delta_ET <- ((delta_initial+1000)*(1-frac)^ae-1000*(1-frac)-delta_initial)/(1-frac-1)*(1-t)+delta_initial*t
  return(c((delta_initial-frac*delta_ET)/(1-frac),delta_ET)) # permil
}

# initialize the δS, the oxygen and hydrogen isotopic compositions of surface water
d18O_S <- d18O_S0 # permil
d2H_S <- d2H_S0 # permil

# initialize the distance vector
x <- seq(1,distance,dx) # m

# calculate the oxygen and hydrogen isotopic compositions of surface water across the 1-D transect
for (i in 1:(length(x)-1)){
  d18O_S[i+1] <- d18O_S[i]+dx*ET_S_ratio/v*(d18O_G-rayleigh_ET(nd,tet,rh,d18O_S[i],a18,k18,smow18)[2]) # permil
  d2H_S[i+1] <- d2H_S[i]+dx*ET_S_ratio/v*(d2H_G-rayleigh_ET(nd,tet,rh,d2H_S[i],a2,k2,smow2)[2]) # permil
}

# plot the oxygen isotopic composition of surface water versus horizontal distance (the control run in Figure 6 of the paper)
plot(x,d18O_S)
