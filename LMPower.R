# Code from Dan Kehler at Quebec-Atlantic Bio-Region, Parks Canada
power.lm<-function(B0=1.5,perc.change=0.05,nsim=1,alpha1=0.05,alpha2=0.1,alpha3=0.2,years=1:5,nreps=1,sd=1,print.p=F){
  
  # The factors to vary are 
  # B0 - the intercept, or starting point (time 0) of the time series 
  # rate - annual change (%)/100 
  # nyears - number of years of data
  # nreps - the number of replicates per year
  # sd - standard deviation of the residuals of the a regression, or of the distribution 
  #   of the normally distributed response
  # alpha = level at which a null hypothesis is rejected.  (1 - confidence) 
  # nsim - the number of simlations used to estimate power.  500 is a minimum.  
  
  # Linear models do *not* have a constant annual percent change.  Instead
  # we look at the percent change (from the starting) point over the duration of 
  # study.  For example, power to detect a 50% change over 5 years. 
  # naturally, the slope needed to effect such a change will differ, depending
  # upon the intercept, and the duration of the study in the following manner:

    B1<- ifelse(perc.change == 0, B0, B0*perc.change / max(years))
  
  
  # It is critical that the right B0 be used!!!!!
    
  p.value<-NULL
  mu <- B0 + B1*years
  year<-as.vector(mapply(rep,years,nreps)) 
  
  for (i in 1:nsim) {

    # Use the mapply function to avoid a loop over the mu vector. 
    # mapply creats a matrix, each column corresponding to a year. 
    # and the number of rows = nreps.  The as.matrix command sticks 
    # each column one after the other in a vector.
    
    response<-as.vector(mapply(rnorm,nreps,mu,MoreArgs=list(sd=sd) ) )
    year<-as.vector(mapply(rep,1:length(mu),nreps)) 
    
    model<-lm(response~year)
    p.value[i]<- anova(model)$"Pr(>F)"[1]
  } # Close nsim loop
  
  # Power is the P(rejecting H0 | H0 is false ) which 
  # translates into the frequency with which a p.value <= alpha is obtained
  power1<- length(p.value[p.value<=alpha1]) / nsim 
  power2<- length(p.value[p.value<=alpha2]) / nsim 
  power3<- length(p.value[p.value<=alpha3]) / nsim 
  if(print.p) print(p.value) 
  return(c(power1,power2,power3))  
}