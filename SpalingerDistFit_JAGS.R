
      model{
      #species hyperparameters
      for(i in 1:n.species){
       mu.b.all[i] ~ dgamma(.001,.001)
       sigma.b.all[i] ~ dt(0,1/y.scale.b[i]^2,1)T(0,)
       mu.sigma.b.all[i] ~ dt(0,1/y.scale.b[i]^2,1)T(0,)
       sigma.sigma.b.all[i] ~ dt(0,1/y.scale.b[i]^2,1)T(0,)
       
       mu.e.all[i] ~ dunif(0,25)
       sigma.e.all[i] ~ dt(0,1/y.scale.e[i]^2,1)T(0,)
       mu.sigma.e.all[i] ~ dt(0,1/y.scale.e[i]^2,1)T(0,)
       sigma.sigma.e.all[i] ~ dt(0,1/y.scale.e[i]^2,1)T(0,)
       
       p[i] ~ dunif(0, 1)
       for(j in 1:y.n.site){
       mu.b[i,j] ~ dgamma(mu.b.all[i]^2 / sigma.b.all[i]^2, mu.b.all[i] / sigma.b.all[i]^2)
        sigma.b[i,j] ~ dgamma(mu.sigma.b.all[i]^2 / sigma.sigma.b.all[i]^2, mu.sigma.b.all[i] / sigma.sigma.b.all[i]^2)
       mu.e[i,j] ~ dgamma(mu.e.all[i]^2 / sigma.e.all[i]^2, mu.e.all[i] / sigma.e.all[i]^2)
       sigma.e[i,j] ~ dgamma(mu.sigma.e.all[i]^2 / sigma.sigma.e.all[i]^2, mu.sigma.e.all[i] / sigma.sigma.e.all[i]^2)
      } #end of j
      } #end of i
      #must assume z and y arise from separate processes
      #likelihood for presence absence data
      for(i in 1:length(z)){
        z[i] ~ dbinom(p[i], y.n.site)
      }
      #likelihood for biomass data conditional on being present
      for(i in 1:length(y)){
      y[i] ~ dgamma(mu.b[y.species.index[i], y.site.index[i]]^2 / (sigma.b[y.species.index[i], y.site.index[i]]^2 + .000001), mu.b[y.species.index[i], y.site.index[i]] / (sigma.b[y.species.index[i], y.site.index[i]]^2 + .00001))
      }
      #likelihood for metabolilzable energy data
      for(i in 1:length(w)){
       
        w[i] ~ dgamma(mu.e[w.species.index[i], w.site.index[i]]^2 / (sigma.e[w.species.index[i], w.site.index[i]]^2 + .000001), mu.e[w.species.index[i], w.site.index[i]] / (sigma.e[w.species.index[i], w.site.index[i]])^2 + .00001)
      }
      for(i in 1:n.species){
          p.test[i] = .00001
      }

  } #end model
      
      
      
