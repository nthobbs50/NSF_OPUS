###############################
#fit_moments function
#Author: N. Thompson Hobbs
#Last revision date 1/26/2023
#Purpose: Fits means and variances for metabolizable energy concentration (kJ / g dry matter) and biomass (kg dry matter / ha) of plant types as described in equations S1 - S9 of  Appendix S1
#dependencies: y.all tibble containing biomass and metabolizable energy data
###############################
#Inputs########################
#ME switch controls wheter the analysis is done for metabolizable energy or nitrogen. It is always set to TRUE for the analysis in the paper, which does not include nitrogren.

#y.all are the data with fields:

# [1] "Veg_type"  The vegetation type analyzed           
# [2] "Species"  Code for plant tissue type       
# [3] "site"  Codes for sites, which are later converted to integer indices                   
# [4] "KgDM_Ha"  Observation of biomass of species at each site in kg dry matter / ha               
# [5] "ME"      Observation of mean metabolizable energy concentration(kJ / gram dry matter) for each species at each site                 
# [6] "ME.SD"   Standard deviation of the distribution of the mean metabolizable energy concentration                 
# [7] "N"   Mean digestible nitrogen content.  Not used in this analysis.                     
# [8] "N.SD" Standard deviation of mean digestible nitrogen content. Not used in this analysis
# [9] "presence"   indicator variable describing whether a species was present at a site = 1 if present, 0 otherwise.              
# [10] "consumer_species"   The herbivore species studied      
# [11] "Veg_type_consumer_species" A code combining vegetation type and herbivore species
# [12] "source"     A code for origin of data              
# [13] "cite"     Citation for data                
# [14] "body_mass"  Body mass (kg) of herbviore studied.    


#n.adapt = 1000, n.update = 100000, n.iter = 100000, n.thin = 10 control the number of iterations (with default values) in the MCMC sampler

#n.chains is the number of chains randomly sampled from the converged coda object


###############################
#Outputs#######################
# A R list object (d) with the following components
#[1] "Veg_type"  Vegetation type

#[2]"body_mass" The body mass of the species studied
#[3] "mu.b.all" A matrix of MCMC estimates of the mean of the distribution of biomass for each species across sites with number of columns equal to number of chains sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[4] "mu.e.all" A matrix of MCMC estimates of the mean of the distribution of metabolizable energy concentration for each species across sites with number of columns equal to number of chains sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[5] "sigma.b.all" A matrix of MCMC estimates of the standard deviation of the distribution of biomass for each species across sites with number of columns equal to number of chains sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S9 in Appendix S1.
#[6] "sigma.e.all" A matrix of MCMC estimates of the standard deviation of the distribution of metabolizable energy concentration for each species across sites with number of columns equal to number of chains sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. See equation S7 in Appendix S1.
#[7] A matrix of MCMC estimates of the probability that a species is present at a site with number of columns equal to number of chains sampled from the converged MCMC output and number of rows equal to the number of species for a particular vegetation type. The equation in Appendix S1 is not numbered but can be found in paragraph 2.
#[8] "max_Rhat"  The maximum R_hat value (S. P. Brooks and A. Gelman. General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, 7:434–455, 1997.) for species in a given vegetation type
#[9] "MCMCsummary" A summary of the coda output for each vegetation type produced by the MCMC summary function of the MCMCvis package
#[10] "out.n.Rhat"  MCMC summary output for species with Rhat values > 1.1
#[11] "n.high.Rhat" The number of species for which Rhat values > 1.1
#[12] "Hobbs_switch" a boolean swith not used in this analysis
#[13] "n.species" The number of species in a vegetation type.
  


fit_moments = function(ME.switch = TRUE, y.all, n.adapt = 1000, n.update = 100000, n.iter = 100000, n.thin = 10, n.sample.chains = 1000 ){
  library(tidyverse)
  library(rjags)
  library(MCMCvis) 
  library(stringr)
  library(dclone)
  library(doParallel)
  library(parallel)  
  set.seed(2)
  
#

# 
# ###fit Shipley distributions############################################
#
Shipley_types = unique( y.all %>% filter(source == "Shipley") %>% pull(Veg_type))
Wagoner_types = unique( y.all %>% filter(source == "Wagoner") %>% pull(Veg_type))
Shipley_types = c(Shipley_types, Wagoner_types)
d = vector("list", length(Shipley_types))
Veg_types = unique(Shipley_types)
for(i in 1:length(Veg_types)){
#   #######
  select_type = Veg_types[i]
  print(select_type)
  ##all data for selected type for presence - absence
  z = filter(y.all, Veg_type == select_type) #all data
  y = z %>% filter(presence == 1) #subset data for sites where a species was present
  #create site and species indices needed by JAGS
  y = y %>% mutate(site.index = as.integer(as.factor(site)), species.index = as.integer(as.factor(Species)))
  z = z%>% mutate(site.index = as.integer(as.factor(site)),  species.index = as.integer(as.factor(Species)))

  #check to be sure that species indexes are aligned
  sort(unique(y$species.index))
  sort(unique(z$species.index))

  sort(unique(y$Species))
  sort(unique(z$Species))
  n.species = length(unique(y$species.index))
  n.site = length(unique(y$site.index))

#   ##get values for inits
  if(ME.switch)  species.means_and_sd = y %>% group_by(Species) %>% summarize(mean.b = mean(KgDM_Ha, na.rm = TRUE), sd.b = sd(KgDM_Ha, na.rm = TRUE), mean.e = mean(ME, na.rm = TRUE), sd.e = sd(ME, na.rm = TRUE), n = n()) %>% ungroup() else species.means_and_sd = y %>% group_by(Species) %>% summarize(mean.b = mean(KgDM_Ha, na.rm = TRUE), sd.b = sd(KgDM_Ha, na.rm = TRUE), mean.e = mean(N, na.rm = TRUE), sd.e = sd(N, na.rm = TRUE), n = n()) %>% ungroup()


  species.site.means_and_sd = y %>% group_by(Species,site) %>% summarize(mean = mean(N, na.rm = TRUE), sd = sd(N, na.rm = TRUE), n = n()) %>% ungroup()

  n.sites.present = species.site.means_and_sd %>% group_by(Species) %>% summarize(n = n())
if(ME.switch){ y.quality = y$ME
y.quality.SD = y$ME.SD
} else {
  y.quality = y$N
  y.quality.SD = y$N.SD
}
  data = list(y = y$KgDM_Ha,  #biomass data
              w = y.quality,  #nutritional quality data, ME in this analysis, referred to as x in Appendix S1
              y.sigma.e = y.quality.SD,
              z = n.sites.present$n, # number of sites present
              y.n.site = length(unique(y$site.index)), #number of sites
              y.species.index = y$species.index,
              z.species.index = z$species.index,
              w.species.index = y$species.index,
              w.site.index = y$site.index,
              y.site.index = y$site.index,
              n.species = length(unique(y$species.index)),
              y.scale.b = species.means_and_sd$mean.b * 1.5, #Scale parameters for half-Cauchy
              y.scale.e = species.means_and_sd$mean.e * 1.5

  )
  base_init= list(mu.b.all = species.means_and_sd$mean.b[1:n.species],
                  sigma.b.all = species.means_and_sd$sd.b[1:n.species],  #assume sd  = mean for inits
                  mu.e.all = species.means_and_sd$mean.e[1:n.species],
                  sigma.e.all = abs(species.means_and_sd$mean.e[1:n.species]) * .5,
                  p = rep(.5, n.species)
  )

  
  make_init = function(multiplier, base = base_init){
    inits = list()
    inits = lapply(base, FUN = function(x) x * multiplier)
    inits$.RNG.name = "base::Mersenne-Twister"
    inits$.RNG.seed= round(runif(1,1,200))
    return(inits)
  }
  inits = list(
    make_init(1.2),
    make_init(1),
    make_init(.8)
  )
  
jags_model_file = "Shipley_dist_fit_JAGS.R" 
  variables = c("p", "mu.b.all", "mu.e.all", "sigma.b.all", "sigma.e.all")


  start.time = Sys.time()
  cl=makeCluster(length(inits))
  registerDoParallel(cl)
  zc=jags.parfit(cl,model=jags_model_file,params=variables,data=data,inits=inits,n.chains=length(inits),n.update=n.update,n.iter=n.iter, n.thin=n.thin)
  stopCluster(cl)
  elapsed.time = Sys.time() - start.time

  out=MCMCsummary(zc, params = c("p", "mu.b.all", "mu.e.all", "sigma.e.all", "sigma.b.all"))
  print(select_type)
  print(max(out$Rhat))
  hist(out$Rhat)
  print(min(out$n.eff))
  print(elapsed.time)
  print("ME next")
  print(MCMCsummary(zc,"mu.e.all")$mean)
  #chains is  the indices of a random sample of size n.sample.chains. Note that mu.b.all is used here simply to get the dimension of the zc object
  chains = sample(1:dim(MCMCchains(zc, params = c("mu.b.all")))[1], n.sample.chains)
  d[[i]]$Veg_type = unique(z$Veg_type_consumer_species)
  d[[i]]$body_mass = unique(z$body_mass)
  d[[i]]$mu.b.all = MCMCchains(zc, params = c("mu.b.all"))[chains,]
  d[[i]]$mu.e.all = MCMCchains(zc, params = c("mu.e.all"))[chains,]
  d[[i]]$sigma.b.all = MCMCchains(zc, params = c("sigma.b.all"))[chains,]
  d[[i]]$sigma.e.all = MCMCchains(zc, params = c("sigma.e.all"))[chains,]
  d[[i]]$p = MCMCchains(zc, params = "p")[chains,]
  d[[i]]$max_Rhat = max(out$Rhat, na.rm=TRUE)
  d[[i]]$MCMCsummary = out
  d[[i]]$out.high.Rhat = out[out$Rhat > 1.1, ]
  d[[i]]$n.high.Rhat = nrow(out[out$Rhat > 1.1, ])
  d[[i]]$season = unique(z$season)
  d[[i]]$Hobbs_switch = FALSE
  d[[i]]$body_mass = z$body_mass
  d[[i]]$n.species = length(unique(y$species.index))
} #end of Shipley fits
if(ME.switch) save(d, file = "Shipley_d_ME.Rdata") else save(d, file = "Shipley_d_N.Rdata")
rm(d)

###fit Spalinger distributions
Spalinger_types = unique( y.all %>% filter(source == "Spalinger") %>% pull(Veg_type))
d = vector("list", length(Spalinger_types))

for(i in 1 : length(Spalinger_types)){
  print(i)
  # #######
  select_type = Spalinger_types[i]
  # ##all data for selected type for presence - absence
  z = filter(y.all, Veg_type == select_type)
  y = z %>% filter(presence == 1)
  
  site.means_and_sd = y %>% group_by(site) %>% summarize(mean = mean(KgDM_Ha, na.rm = TRUE), sd = sd(KgDM_Ha, na.rm = TRUE), n = n()) %>% ungroup()
  
  species.site.means_and_sd = y %>% group_by(Species, site) %>% summarize(mean = mean(KgDM_Ha, na.rm = TRUE), sd = sd(KgDM_Ha, na.rm = TRUE), n = n()) %>%
    ungroup()
  
  n.sites.present = species.site.means_and_sd %>% group_by(Species) %>% summarize(n = n())
  
    species.site.means_and_sd = y %>% group_by(Species, site) %>% summarize(mean = mean(KgDM_Ha, na.rm = TRUE), sd = sd(KgDM_Ha, na.rm = TRUE), n = n()) %>%
    ungroup()
  
  n.sites.present = species.site.means_and_sd %>% group_by(Species) %>% summarize(n = n())
  
  
  
  y = y %>% mutate(site.index = as.integer(as.factor(site)), species.index = as.integer(as.factor(Species)))
  z = z%>% mutate(site.index = as.integer(as.factor(site)),  species.index = as.integer(as.factor(Species)))
  n.site = length(unique(z$site.index))
  
  
  
  #check to be sure that species indexes are aligned
  sort(unique(y$species.index))
  sort(unique(z$species.index))
  
  sort(unique(y$Species))
  sort(unique(z$Species))
  n.species = length(unique(y$species.index))
  print(n.species)
  n.site = length(unique(y$site.index))
  
  if(ME.switch) y.quality = y$ME else y.quality = y$N
  
  ##get values for inits
  if(ME.switch)  species.means_and_sd = y %>% group_by(Species) %>% summarize(mean.b = mean(KgDM_Ha, na.rm = TRUE), sd.b = sd(KgDM_Ha, na.rm = TRUE), mean.e = mean(ME, na.rm = TRUE), sd.e = sd(ME, na.rm = TRUE), n = n()) %>% ungroup() else species.means_and_sd = y %>% group_by(Species) %>% summarize(mean.b = mean(KgDM_Ha, na.rm = TRUE), sd.b = sd(KgDM_Ha, na.rm = TRUE), mean.e = mean(N, na.rm = TRUE), sd.e = sd(N, na.rm = TRUE), n = n()) %>% ungroup()
  
  
  data = list(y = y$KgDM_Ha,
              w = y.quality,
              z = n.sites.present$n,
              y.species.index = y$species.index,
              z.species.index = z$species.index,
              w.species.index = y$species.index,
              y.site.index = y$site.index,
              w.site.index = y$site.index,
              n.species = n.species,
              y.scale.b = species.means_and_sd$mean.b * 1.5,  # scale parameter for half Cauchy
              y.scale.e = species.means_and_sd$mean.e * 1.5, # scale parameter for half Cauchy
              y.n.site = n.site
              
  )
  base_init= list(mu.b.all = species.means_and_sd$mean.b[1:n.species],
                  sigma.b.all = species.means_and_sd$sd.b[1:n.species],  #assume sd  = mean for inits
                  mu.e.all = species.means_and_sd$mean.e[1:n.species],
                  sigma.e.all = abs(species.means_and_sd$mean.e[1:n.species] + .0001) *.3,
                  mu.n = matrix(.01, nrow = length(unique(y$Species)), ncol = n.site),
                  sigma.n = matrix(.001, nrow = length(unique(y$Species)), ncol = n.site),
                  mu.sigma.e.all = abs(species.means_and_sd$mean.e[1:n.species] *.5) + .0001,
                  p = rep(.5, n.species)
  )
  if(!ME.switch) base_init = lapply(base_init,"*",.003)
  
  make_init = function(multiplier, base = base_init){
    inits = list()
    inits = lapply(base, FUN = function(x) x * multiplier)
    inits$.RNG.name = "base::Mersenne-Twister"
    inits$.RNG.seed= round(runif(1,1,200))
    return(inits)
  }
  inits = list(
    make_init(1.2),
    make_init(1),
    make_init(.8)
  )

  if(ME.switch) jags_model_file = "SpalingerDistFit_JAGS.R" else jags_model_file = "SpalingerDistFit_N.R"
  variables = c("p", "mu.b.all", "mu.e.all", "sigma.b.all", "sigma.e.all")
  
  
  start.time = Sys.time()
  cl=makeCluster(length(inits))
  registerDoParallel(cl)
  zc=jags.parfit(cl,model=jags_model_file,params=variables,data=data,inits=inits,n.chains=length(inits),n.update=n.update,n.iter=n.iter, n.thin=n.thin)
  stopCluster(cl)
  elapsed.time = Sys.time() - start.time
  
  out=MCMCsummary(zc, params = c("p", "mu.b.all", "mu.e.all", "sigma.e.all", "sigma.b.all"))
  print(select_type)
  print(max(out$Rhat))
  print(min(out$n.eff))
  print(elapsed.time)
  print("mean quality next")
  print(MCMCsummary(zc,"mu.e.all")$mean)
  # out[out$Rhat > 1.1,]
  #chains is a random sample of the indices of chains to be output of length n.sample.chains
  chains = sample(1:dim(MCMCchains(zc, params = c("mu.b.all")))[1], n.sample.chains)
  d[[i]]$Veg_type = unique(z$Veg_type_consumer_species)
  d[[i]]$body_mass = unique(z$body_mass)
  d[[i]]$mu.b.all = MCMCchains(zc, params = c("mu.b.all"))[chains,]
  d[[i]]$mu.e.all = MCMCchains(zc, params = c("mu.e.all"))[chains,]
  d[[i]]$sigma.b.all = MCMCchains(zc, params = c("sigma.b.all"))[chains,]
  d[[i]]$sigma.e.all = MCMCchains(zc, params = c("sigma.e.all"))[chains,]
  d[[i]]$p = MCMCchains(zc, params = "p")[chains,]
  d[[i]]$max_Rhat = max(out$Rhat, na.rm=TRUE)
  d[[i]]$MCMCsummary = out
  d[[i]]$out.high.Rhat = out[out$Rhat > 1.1, ]
  d[[i]]$n.high.Rhat = nrow(out[out$Rhat > 1.1, ])
  d[[i]]$season = unique(z$season)
  d[[i]]$Hobbs_switch = FALSE
  d[[i]]$body_mass = z$body_mass
  d[[i]]$n.species = n.species
  if (ME.switch) save(d, file = "Spalinger_d_ME.Rdata") else save(d, file = "Spalinger_d_N.Rdata")
} #end of Spalinger_fits
rm(d)


# #
# #

} #end of function
#
# #test function
#
#fit_momnets(ME.switch = TRUE, y.all = y.all, n.adapt = 10, n.iter = 10, n.update = 10, n.sample.chains = 1, n.thin = 1)
#
