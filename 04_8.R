#####################################################
# Figure 4.8. Sediment mineralization model
#####################################################

#General settings
library(deSolve)

# Function to calculate model cost
costf <- function(params){
  with(as.list(params), {
    Carbon   <- meanDepo * mult / k
    outtimes <- as.vector(oxcon$time)
    outmin   <- ode(Carbon, outtimes, minmod, params)
    sqdiff   <- (outmin[,3] - oxcon$cons)^2
    
    sum(sqdiff)
    })
}

# Function to calculate derivative of state variable
minmod <- function(t,Carbon,parameters){ 
  with (as.list(c(Carbon,parameters)),{
    
    minrate <- k*Carbon
    Depo <- approx(Flux[,1],Flux[,2], xout=t)$y
    dCarbon <- mult*Depo - minrate
    
    list(dCarbon,minrate)
  })
}

# initial par estimates
# minimal parameter values
# maximal parameter values
# function to minimise
# nr elements in population
# number of iterations
# number of points in centroid
# relative variation upon stopping
pricefit <- function (par, minpar=rep(-1e8,length(par)), 
                      maxpar=rep(1e8,length(par)), func, npop=max(5*length(par),50), 
                      numiter=10000, centroid = 3, varleft = 1e-8, ...){

  # Initialization
cost <- function (par) func(par,...)
npar <- length(par)
tiny <- 1e-8
varleft <- max(tiny,varleft)
populationpar <- matrix(nrow=npop,ncol=npar,byrow=TRUE, 
                        data= minpar+runif(npar*npop)*rep((maxpar-minpar),npop)) 

colnames(populationpar)<-names(par)
populationpar[1,]<-par
populationcost <- apply(populationpar,FUN=cost,MARGIN=1)
iworst <- which.max(populationcost)
worstcost <- populationcost[iworst]

# Hybridization phase
iter<-0
while(iter < numiter & 
       (max(populationcost)-min(populationcost)) > (min(populationcost)*varleft)){

  iter<-iter+1
  selectpar <- sample(1:npop,size=centroid)
 # for cross-fertilization
  mirrorpar <- sample(1:npop,size=1)
 # for mirroring
  newpar <- colMeans(populationpar[selectpar,]) # centroid
  newpar <- 2*newpar - populationpar[mirrorpar,] # mirroring
  newpar <- pmin( pmax(newpar,minpar) ,maxpar)
  newcost <- cost(newpar)

  if(newcost < worstcost){
      populationcost[iworst] <-newcost
      populationpar [iworst,]<-newpar
      iworst <- which.max(populationcost) # new worst member
      worstcost <- populationcost[iworst]
    }
} # end j loop

ibest    <- which.min(populationcost)
bestpar  <- populationpar[ibest,]
bestcost <- populationcost[ibest]

return (list(par = bestpar, cost = bestcost, 
             poppar = populationpar, popcost=populationcost))
}

# Define problem and data
Flux <- matrix(ncol=2,byrow=TRUE,data=c(
1,
 0.654, 11, 0.167, 21, 0.060, 41, 0.070,
73,
 0.277, 83, 0.186, 93, 0.140,103, 0.255,
113,
 0.231,123, 0.309,133, 1.127,143, 1.923,
153,
 1.091,163, 1.001,173, 1.691,183, 1.404,
194,
 1.226,204, 0.767,214, 0.893,224, 0.737,
234,
 0.772,244, 0.726,254, 0.624,264, 0.439,
274,
 0.168,284, 0.280,294, 0.202,304, 0.193,
315,
 0.286,325, 0.599,335, 1.889,345, 0.996,
355,
 0.681,365, 1.135))
meanDepo <- mean(approx(Flux[,1],Flux[,2], xout=seq(1,365,by=1))$y)
oxcon<-as.data.frame(matrix(ncol=2,byrow=TRUE,data=c(
68, 0.387, 69, 0.447, 71, 0.473, 72, 0.515,
189, 1.210,190, 1.056,192, 0.953,193, 1.133,
220, 1.259,221, 1.291,222, 1.204,230, 1.272,
231, 1.168,232, 1.168,311, 0.963,312, 1.075,
313, 1.023)))
names(oxcon)<-c("time","cons")

plot(oxcon)

multser <- seq(1,1.5,by=.05)
numms   <- length(multser)
kseries <- seq(0.001,0.05,by=0.002)
numks   <- length(kseries)

outcost <- matrix(nrow = numms, ncol = numks)
for (m in 1:numms){
  for (i in 1:numks){
    pars         <- c(k = kseries[i], mult = multser[m])
    outcost[m,i] <- costf(pars)
  }
}

minpos <- which(outcost==min(outcost),arr.ind=TRUE)
multm  <- multser[minpos[1]]
ki     <- kseries[minpos[2]]

optpar <- pricefit(par=c(k=ki,mult=multm),minpar=c(0.001,1), 
                   maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500, 
                   centroid=3,varleft=1e-8)

optpar20 <- pricefit(par=optpar$par,minpar=c(0.001,1),
maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500,
centroid=3,varleft=0.2)

optpar25 <- pricefit(par=optpar$par,minpar=c(0.001,1),
maxpar=c(0.05,1.5),func=costf,npop=50,numiter=500,
centroid=3,varleft=0.025)


outtimes <- seq(1,365,by=1)
Carbon <- meanDepo*optpar$par[2]/optpar$par[1]
names(Carbon) <-"Carbon"
out <- as.data.frame(ode(Carbon,outtimes,minmod, optpar$par))
names(out) <- c("time","Carbon","minrate")
par (oma=c(0,0,0,2))
plot(Flux,type="l",xlab="daynr",ylab="mmol/m2/d",
main="Sediment-detritus model",lwd=2)

lines(out$time,out$minrate,lwd=2,col="darkgrey")
points(oxcon$time,oxcon$cons,pch=25,col="black", bg="darkgray",cex=2)
par(new=TRUE)
plot(out$time,out$Carbon,axes=FALSE,xlab="",ylab="",
type="l",lty=2)
axis(4)
mtext(side=4,"mmolC/m2",outer=TRUE)
legend("topleft",col=c("black","darkgrey","black"),
leg=c("C flux","C mineralization","C concentration"),
lwd=c(2,2,1),lty=c(1,1,2))
