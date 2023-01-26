#####
#get_dist function
#Author: N. Thompson Hobbs
#Last revision date 1/26/2023
#Purpose: Implements equations S10 - S15 to obtain d vector approximating out of sample data as described in Appendix S1. The kernal density of d approximates the nutrient distribution [x | \theta].
#Dependencies: Output of fit_moments function
#Inputs########################
#[1] "mu.b.all" A matrix of MCMC estimates of the mean of the distribution of biomass for each species across sites with number of columns equal to number of iterations sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[2] "mu.e.all" A matrix of MCMC estimates of the mean of the distribution of metabolizable energy concentration for each species across sites with number of columns equal to number of iternations sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[3] "sigma.b.all" A matrix of MCMC estimates of the standard deviation of the distribution of biomass for each species across sites with number of columns equal to number of iterations sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S9 in Appendix S1.
#[4] "sigma.e.all" A matrix of MCMC estimates of the standard deviation of the distribution of metabolizable energy concentration for each species across sites with number of columns equal to number of iterations sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[5] A matrix of MCMC estimates of the probability that a species is present at a site with number of columns equal to number of iterations sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. The equation in Appendix S1 is not numbered but can be found in paragraph 2.
#[6] tilde.K The number of random draws of draws of metabolizable energy concentration, biomass, and probability of occurrence for each plant type for each of the sample of MCMC chains. See paragraph two on page 3 of Appendix S1.

#[7] MEflag = TRUE is a boolean value that determines whether plant nutritional quality is metabolizable energy or digestible nitrogen.  It is set at TRUE for this analysis because it includes only metabolizable energy

#Hobbs_switch a boolean value that is set at FALSE for this analysis

#Ouputs########################
#me.dist The vector d for each iteration. See Appendix S1 paragraph 1 page 4.
#total.biomass TA derived quantity approximating the marginal posterior distribution of the sum of species biomass.
#n.species The number of species for a given vegetation type
library(tidyverse)
get_dist = function(mu.b.all, sigma.b.all, mu.e.all, sigma.e.all, p, tilde.K = 1000, MEflag = TRUE, Hobbs_switch = FALSE){
 library(LaplacesDemon)
 n.species = dim(mu.b.all)[2] #number of plant species
 n.chains = dim(mu.b.all)[1]
 #set upper and lower limits on of possible values for metabolizable energy and digestible nitrogen
 if(MEflag){
     qual.limit.high = 15
     qual.limit.low = 2}
  else 
  {qual.limit.high = .05  
  qual.limit.low = -.04}

 
 if(MEflag) dist.type = "gamma"
 else dist.type = "norm"

  mean.mass = present = p.biomass = species.biomass = draw.me = array(dim=c(n.species, tilde.K, n.chains))
  p.test = rep(.00001, n.species)
  total.biomass = me.index = matrix(nrow = tilde.K,ncol = n.chains)
  me.dist = matrix(ncol = tilde.K, nrow = n.chains)
  #use this to prevent all zero probabilities in dcat vector
  for(t in 1:n.chains){
    for(k in 1:tilde.K){
      for(i in 1:n.species){
        if(MEflag){
          param1 = mu.e.all[t,i]^2 / (sigma.e.all[t,i]^2 + .00001)
          param2 = mu.e.all[t,i] / (sigma.e.all[t,i]^2 + .00001)
        } else {
          param1 = mu.e.all[t,i]
          param2 = sigma.e.all[t,i]
        }
        mean.mass[i,k,t] = rgamma(1,mu.b.all[t,i]^2 / (sigma.b.all[t,i]^2 + .0001), mu.b.all[t,i]/ (sigma.b.all[t,i]^2 + .0001)) #draw of unobserved average mass of species i at site K conditional on being presentmean.mass[i,k,t]. Eq. S11 in appendix S1
        if(Hobbs_switch) present[i,k,t] = 1 else{
        present[i,k,t] = rbinom(1, prob=p[t,i], size=1)}
      species.biomass[i,k,t] <- present[i,k,t] * mean.mass[i,k,t] 
      #random draw from out of sample site k for species i from chain t.  #truncate gamma distribution (dist.type) to prevent unreasonable values for nutrient concentration. Eq S13 in Appendix S1
      draw.me[i,k,t] = rtrunc(1, a=qual.limit.low, b = qual.limit.high, spec = dist.type, param1, param2) #Eq. S12 in Appendix S1
      
       } #end of i
  #
  total.biomass[k,t] = sum(species.biomass[,k,t]) 
  p.biomass[1:n.species,k,t] = (species.biomass[1:n.species,k,t] / (total.biomass[k,t] + .0000001)) + p.test[] ##Eq. S14 in Appendix S1. p.test prevents all zero probabilities
  #get random draw of index of a species at site k for chain t where the probability is the proportion of the species in the site biomass
  me.index[k,t] = extraDistr::rcat(1, prob = p.biomass[1:n.species,k,t]) #Eq S15 in Appendix S1
  me.dist[t,k] = draw.me[me.index[k,t],k,t] #put chains in rows, out of sample predictions of  ME in columns. One tildeK length vector for each iteration. So for each  iteration there will be tilde.K draws of ME concentration in proportion to the contribution of the species to the total biomass.

     } #end k
    
} # end of t
out = list(me.dist = me.dist, total.biomass = colMeans(total.biomass), n.species = n.species)
return(out)
 } #end of function

#test function
#nut_dist = get_dist(mu.b.all = d[[i]]$mu.b.all, sigma.b.all = d[[i]]$sigma.b.all, mu.e.all = d[[i]]$mu.e.all, sigma.e.all = d[[i]]$sigma.e.all, p = d[[i]]$p, tilde.K = 10, MEflag = FALSE, Hobbs_switch = FALSE)
# 
