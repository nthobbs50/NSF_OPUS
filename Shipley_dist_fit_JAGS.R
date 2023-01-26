###JAGS code implementing Equations S1 - S9 in Appendix S1.

model{
  #species hyperparameters
  for(i in 1:n.species){
    mu.b.all[i] ~ dgamma(.001,.001)
    sigma.b.all[i] ~ dt(0,1/y.scale.b[i]^2,1)T(0,)
    
    mu.e.all[i] ~ dgamma(.001,.001)T(.00001, )
    sigma.e.all[i] ~ dt(0,1/y.scale.e[i]^2,1)T(0,)
    
    p[i] ~ dunif(0, 1)
   
    for(j in 1:y.n.site){
           mu.e[i,j] ~ dgamma(mu.e.all[i]^2/sigma.e.all[i]^2, mu.e.all[i]/sigma.e.all[i]^2)
     } #end of j
  } #end of i
  ############################
  
  #likelihood for presence absence data
  for(i in 1:length(z)){
    z[i] ~ dbinom(p[i], y.n.site)
  }
  #likelihood for biomass data conditional on being present Eq. S7 in appendix S1
  for(i in 1:length(y)){
    y[i] ~ dgamma(mu.b.all[y.species.index[i]]^2 / (sigma.b.all[y.species.index[i]]^2 + .000001), mu.b.all[y.species.index[i]] / (sigma.b.all[y.species.index[i]]^2 + .00001))
  }
  #likelihood for metabolizable energy data conditional on being present. Eq. S2 in appendix S1
    for(i in 1:length(w)){
        w[i] ~ dgamma(mu.e[w.species.index[i], y.site.index[i]]^2 / (y.sigma.e[w.species.index[i]] + .0001 )^2, mu.e[w.species.index[i], y.site.index[i]] / (y.sigma.e[w.species.index[i]])^2)
    }
  
} #end model

