##### SCRIPT INFORMATION

# XXX.R

# XXX

##########################

##### DEVELOPER NOTES

# N/A

##########################

##### SCRIPT USER INPUT

### Environment

rm(list = ls())
#try(dev.off())
setwd('~/Desktop/')

# Accessory packages
library(deSolve)

### Key variables

# Model rate constants
k.1=1e6 # monomer-1 association rate constant, 1/M/s
k.n1=16e-4 # monomer-1 dissociation rate constant, 1/s
k.2=1e-5 # 1:1 complex forward isomerization rate constant, 1/s
k.n2=32e-4 # 1:1 complex backward isomerization rate constant, 1/s
k.3=1e9 # monomer-2 association rate constant, 1/M/s
k.n3=2e-2 # monomer-2 dissociation rate constant, 1/s
k.alpha=20e6 # rate constant for isomerized complex catalysis of complex isomerization, 1/M/s
k.beta=2e6 # rate constant for dimer catalysis of complex isomerization, 1/M/s
k.3b=1e4 # monomer-2 alternate association rate constant, 1/M/s

# Reactant concentrations
E.T=c(1e3/2^c(0:12),0)*1e-9 # total protein concentration, M
L.T=5e-9 # total prey concentration, M

### Extraneous variables

# Simulation parameters
ass.time=60*60*1.5 # total association time, s
diss.time=60*60*1 # total dissociation time, s
dt=25e-3 # integration time-step, s
#
integrator.method='lsoda'
#
j.fold=1e6 # fold jump-dilution

# Data features
save.id='NULL' # identifier for save file(s)

##########################

show(date())

##### BEGIN SCRIPT

# Save simulation parameters
parameters=mget(ls())

# Record equations
Equations=function(t,initials,constants){
  with(as.list(c(initials,constants)),{
    
    dE = k.n1*EL+k.n3*E2L-E*(k.1*L+k.3*ELast)
    dL = k.n1*EL-k.1*E*L
    dEL = k.1*E*L+k.n2*ELast-EL*(k.n1+k.2+k.alpha*ELast+k.beta*E2L+k.3b*E)
    dELast = (k.2+k.alpha*ELast+k.beta*E2L)*EL+k.n3*E2L-ELast*(k.n2+k.3*E)
    dE2L = k.3*E*ELast+k.3b*EL*E-k.n3*E2L

    return(list(c(dE,dL,dEL,dELast,dE2L)))
  })
}

# Additional functions
#none#

### Reaction Simulator

ALL.spr=list(NULL)

for (i in 1:length(E.T)){

# Record constants
constants=c(k.1=k.1,k.n1=k.n1,k.2=k.2,k.n2=k.n2,k.3=k.3,k.n3=k.n3,k.3b=k.3b,k.alpha=k.alpha,k.beta=k.beta)

# Record initial values
initials=c(E=E.T[i],L=L.T,EL=0,ELast=0,E2L=0)

# Reaction
sim.sim1=as.data.frame(ode(y=initials,times=seq(0,ass.time,dt),func=Equations,parms=constants,method = integrator.method))

# Record initial values
initials=c(E=sim.sim1$E[ass.time/dt+1]/j.fold,L=sim.sim1$L[ass.time/dt+1]/j.fold,EL=sim.sim1$EL[ass.time/dt+1]/j.fold,ELast=sim.sim1$ELast[ass.time/dt+1]/j.fold,E2L=sim.sim1$E2L[ass.time/dt+1]/j.fold)

# Reaction
sim.sim2=as.data.frame(ode(y=initials,times=seq(0,diss.time,dt),func=Equations,parms=constants,method = integrator.method))

# FA calculations
SPR.data=rbind(sim.sim1,sim.sim2[-1,])
SPR.data$time=seq(0,ass.time+diss.time,dt)
SPR.data$FA=(SPR.data$E2L+0.5*(SPR.data$ELast+SPR.data$EL))/(SPR.data$E2L+SPR.data$ELast+SPR.data$EL+SPR.data$L)
ALL.spr[[i]]=SPR.data

}

# Plotting
par(fig=c(0,0.5,0.5,1))
plot(NULL,NULL,ylim=c(0,1),xlim=c(0,ass.time/60),main='Association',xlab='Time (min)',ylab='Relative Anisotropy',new=TRUE)
for (i in 1:length(ALL.spr)){
  SPR.data=ALL.spr[[i]]
  lines(SPR.data$time/60,SPR.data$FA,col=i)
}
#
par(fig=c(0.5,1,0.5,1),new=T)
plot(NULL,NULL,ylim=c(0,1),xlim=c(0,diss.time/60),main='Dissociation',xlab='Time (min)',ylab='Relative Anisotropy',new=TRUE)
for (i in 1:length(ALL.spr)){
  SPR.data=ALL.spr[[i]]
  lines(SPR.data$time/60-ass.time/60,SPR.data$FA,col=i)
}
#
par(fig=c(0.2,0.8,0,0.5),new=T)
FPeq=rep(NA,times=length(E.T))
for (i in 1:length(ALL.spr)){
  FPeq[i]=(ALL.spr[[i]])$FA[ass.time/dt]
}
plot(log10(E.T),FPeq,type='b',ylim = c(0,1),main='Equilibrium Binding Curve',xlab = expression('log'[10]*'[E'['T']*'] (M)'),ylab = 'Relative Anisotropy')

##### END SCRIPT

show(date())








