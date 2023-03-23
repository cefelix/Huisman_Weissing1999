#initializing required packages

rm(list=ls())
library("deSolve")
library("dplyr")


####initializing model parameters####

K <- c(1, 0.9, 0.3, 1.04, 0.34, 0.77)
K <- rbind(K, c(0.3, 1, 0.9, 0.71, 1.02, 0.76)) 
K <- rbind(K, c(0.9, 0.3, 1, 0.46, 0.34, 1.07))
K                 #half saturation constants of for each species on each resource 
#(from figure 1c in Huisman Weissing 1999)

C <- c(0.04, 0.07, 0.04, 0.1, 0.03, 0.02)
C <- rbind(C, c(0.08, 0.08, 0.1, 0.1, 0.05, 0.17))
C <- rbind(C, c(0.14, 0.1, 0.1, 0.16, 0.06, 0.14))
C                 #resource content matrix for each species on each resource
#(from figure 1c in Huisman Weissing 1999)

S <- c(6, 10, 14) #Supply concentration of resources 1 to 3
r = 1             #maximum specific growth (same for all species)
m = 0.25          #specific mortality rate (same for all species)
D = m             #System's turnover rate

i=6               #number of plankton species 
j=3               #number of resources the plankton depends on

u=NULL
N=NULL
R=NULL




####building the Huisman-Weissing-model from figure 1c####

#blueprint for similar models where i species depend j resources
#simply add phytoplankton growth terms, resource depletion terms, and update C and K matrices

HW.comp.i_on_j <- function(t, V, parms){
  with(as.list(parms),{
    #initial densities of plankton species (N) and resources (R)
    N[1] = V[1]
    N[2] = V[2]
    N[3] = V[3]
    N[4] = V[4]
    N[5] = V[5]
    N[6] = V[6]
    
    R[1] = V[7]
    R[2] = V[8]
    R[3] = V[9]
    
    
    #minimal growth rates for n species
    #equation (3) in Huisman-Weissing 1999
    for (n in 1:i) {
      
      u.in=NULL
      for (k in 1:j) {
        u.in = rbind(u.in, (r*R[k])/(K[k,n] + R[k]))
        
      }
      u[n] = min(u.in)
    }
    
    
    #sum of c*u*N for j resources 
    #sum term of equation (2) 
    sum.cuN = vector(, length=j)
    for (k in 1:j) {
      for (n in 1:i) {
        sum.cuN[k] = sum.cuN[k] + C[k,n]*N[n] * u[n]
      }
    }
    
    
    #phytoplankton growth equations
    #equation (1) 
    
    dN1dt = N[1] * (u[1] - m)
    dN2dt = N[2] * (u[2] - m)
    dN3dt = N[3] * (u[3] - m)
    dN4dt = N[4] * (u[4] - m)
    dN5dt = N[5] * (u[5] - m)
    dN6dt = N[6] * (u[6] - m)
    
    
    #resource depletion equations
    #equation (2)
    
    dR1dt = D * (S[1] - R[1]) - sum.cuN[1]
    dR2dt = D * (S[2] - R[2]) - sum.cuN[2]
    dR3dt = D * (S[3] - R[3]) - sum.cuN[3]
    
    return(list( c(dN1dt, dN2dt, dN3dt, dN4dt, dN5dt, dN6dt, dR1dt, dR2dt, dR3dt) ))
  })
}




####applying the model from fig. 1c with the parameters given in Huisman-Weissing 1999####

#species 1-3 starting at t=0

outHW.comp.i_on_j <- as.data.frame(
  lsoda( y=     c(N1=0.11,
                  N2=0.12,
                  N3=0.13,
                  N4=0,
                  N5=0,
                  N6=0,
                  R1=6,
                  R2=10,
                  R3=14) ,
         times= seq(0, 1000, by=1),
         func = HW.comp.i_on_j,
         parms= c(r=r,
                  K=K,
                  C=C,
                  m=m,
                  D=D,
                  i=i,
                  j=j)
  )
)

#add species 4 at t=1000

#always leaving out the last row of previous data frame while using rbind()
#as it would occur double otherwise

outHW.comp.i_on_j <- rbind(outHW.comp.i_on_j[1:1000,],
                           as.data.frame(
                             lsoda( y=     c(N1=outHW.comp.i_on_j$N1[1001],
                                             N2=outHW.comp.i_on_j$N2[1001],
                                             N3=outHW.comp.i_on_j$N3[1001],
                                             N4=0.1,
                                             N5=0,
                                             N6=0,
                                             R1=outHW.comp.i_on_j$R1[1001],
                                             R2=outHW.comp.i_on_j$R2[1001],
                                             R3=outHW.comp.i_on_j$R3[1001]) ,
                                    times= seq(1000, 2000, by=1),
                                    func = HW.comp.i_on_j,
                                    parms= c(r=r,
                                             K=K,
                                             C=C,
                                             m=m,
                                             D=D,
                                             i=i,
                                             j=j) 
                             )
                           )
)

#add species 5 at t=2000
outHW.comp.i_on_j <- rbind(outHW.comp.i_on_j[1:2000,],
                           as.data.frame(
                             lsoda( y=     c(N1=outHW.comp.i_on_j$N1[2001],
                                             N2=outHW.comp.i_on_j$N2[2001],
                                             N3=outHW.comp.i_on_j$N3[2001],
                                             N4=outHW.comp.i_on_j$N4[2001],
                                             N5=0.1,
                                             N6=0,
                                             R1=outHW.comp.i_on_j$R1[2001],
                                             R2=outHW.comp.i_on_j$R2[2001],
                                             R3=outHW.comp.i_on_j$R3[2001]) ,
                                    times= seq(2000, 5000, by=1),
                                    func = HW.comp.i_on_j,
                                    parms= c(r=r,
                                             K=K,
                                             C=C,
                                             m=m,
                                             D=D,
                                             i=i,
                                             j=j) 
                             )
                           )
)

#add species 6 at t=5000
outHW.comp.i_on_j <- rbind(outHW.comp.i_on_j[1:5000,],
                           as.data.frame(
                             lsoda( y=     c(N1=outHW.comp.i_on_j$N1[5001],
                                             N2=outHW.comp.i_on_j$N2[5001],
                                             N3=outHW.comp.i_on_j$N3[5001],
                                             N4=outHW.comp.i_on_j$N4[5001],
                                             N5=outHW.comp.i_on_j$N5[5001],
                                             N6=0.1,
                                             R1=outHW.comp.i_on_j$R1[5001],
                                             R2=outHW.comp.i_on_j$R2[5001],
                                             R3=outHW.comp.i_on_j$R3[5001]) ,
                                    times= seq(5000, 15000, by=1),
                                    func = HW.comp.i_on_j,
                                    parms= c(r=r,
                                             K=K,
                                             C=C,
                                             m=m,
                                             D=D,
                                             i=i,
                                             j=j) 
                             )
                           )
)

####3.0 plot creation####
#shortcut to get data
outHW.comp.i_on_j <- read.csv("./output6on3.csv")


####3.1 abundances of the plankton species over time (fig. 1c)####
par(mfrow=c(1,1))
cex=1.8
plot(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N1, type="l", col="red",
     xlab="time (days)", ylab="species abundance", xlim =c(0, 15000),
     cex.lab = cex, cex.axis=cex)
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N2, col="darkgreen")
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N3, col="orange")
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N4, col="darkblue")
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N5, col="purple")
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$N6, col="lightblue")


####3.2 resource densities between t=10000 and t=10100 ####
plot(outHW.comp.i_on_j$time, outHW.comp.i_on_j$R1, type="l", col="orange",
     xlab="time (days)", ylab="resource content", xlim =c(10000, 10100), ylim = c(0,1))
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$R2, col="cyan")
lines(outHW.comp.i_on_j$time, outHW.comp.i_on_j$R3, col="violet")



####3.3 Limit cycles of species 1-6 on Resource 3####
par(mfrow=c(3,2))
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N1, type="l", col="red",
     xlab="Resource 3", ylab="abundance spec 1",
     cex.lab=cex)
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N2, type="l", col="darkgreen",
     xlab="Resource 3", ylab="abundance spec 2",
     cex.lab=cex)
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N3, type="l", col="orange",
     xlab="Resource 3", ylab="abundance spec 3",
     cex.lab=cex)
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N4, type="l", col="darkblue",
     xlab="Resource 3", ylab="abundance spec 4",
     cex.lab=cex)
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N4, type="l", col="purple",
     xlab="Resource 3", ylab="abundance spec 5",
     cex.lab=cex)
plot(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N6, type="l", col="lightblue",
     xlab="Resource 3", ylab="abundance spec 6",
     cex.lab=cex)


####3.4 Limit cycles of species 4,5,6 on Resource 2####
par(mfrow=c(1,1))

plot(outHW.comp.i_on_j$R2, outHW.comp.i_on_j$N5, type="l", col="darkblue",
     xlab="Resource 2", ylab="Species abundance", ylim = c(0,0.2), xlim = c(0,0.75))
#lines(outHW.comp.i_on_j$R3, outHW.comp.i_on_j$N4, col="darkblue")
lines(outHW.comp.i_on_j$R2, outHW.comp.i_on_j$N4, col="purple")
lines(outHW.comp.i_on_j$R2, outHW.comp.i_on_j$N6, col="lightblue")

