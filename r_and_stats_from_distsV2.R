#####
#get_r function
#Author: N. Thompson Hobbs
#Last revision date 1/26/2023
#Purpose: Returns parameters of population growth given the vector d and a proportion of the biomass consumed.  It also returns statistics describing the distribution.  Implements steps described in paragraph two, page 4 of Appendix S1.
#Dependencies: Output of fit_moments function, get_dist function

#Inputs########################
#[1] fx is the vector d estimated in the get_dist function and described on pages 2 - 4 of Appendix S1. 
#[2] total.biomass is the estimate of total biomass across species obtain from the get_dist function. 
#[3] n.species is the number of plant species (and plant parts) in a vegetation type.  Not used in this analysis.
#[4] Body mass is the mass in kg of a herbivore species.
#[5] phi is the proportion of biomass consumed. Eq. 6 in the manuscript
#[6] eta is the rate of daily intake of plant dry matter as a proportion of body mass derived from literature values described on lines 268 - 269 of the manuscript.
#[7] act_mul is the mutiplicative cost of activity over basal metabolism dervied from literature values describe on lines 262 - 263 of the manuscript
#[7] ME.swtich is a boolean value that detemines whether computations are based on metabolizable energy or digestible nitrogen.  It is set a TRUE because the all of the analyses in the manuscript are based on metabolizable energy.
library(tidyverse)
  get_r = function(fx, total.biomass, n.species = NULL, body_mass, phi, eta = .02145, act_mul = 2.08, ME.switch = TRUE){ 
    #eta is from weighted mean method from Meta_analysis_of_DMI.R
        N.days.km = (phi * total.biomass) / (eta * body_mass) * 100  #animal days per km^2. Not used in manuscript
      N = N.days.km / 365  #animals per km^2 for 365 days
      z = quantile(fx, 1 - phi) #Approximation of Eq. 8 in manuscript
      mu.eat = mean(fx[fx >= z]) #Approximation of Eq. 9 in manuscript
      if(ME.switch){
        e_intake = eta * body_mass * mu.eat  * 1000
        e_req = (act_mul * 70 * body_mass ^ 0.75) * 4.184
        q = e_intake / e_req #Eq. 10 in manuscript
      }
      else{
        n_intake = eta * body_mass * mu.eat  * 1000 #mu.eat is in gm N / gm dry matter
        endog_fecal_N = eta * body_mass * 5.6 #gms N per kg intake, Robbins and Moen
        endog_urinary_N = 0.16*body_mass^0.75  #Robbins citation,
        n_req = endog_fecal_N + endog_urinary_N
        q = n_intake / n_req
      
        }
      #from Sinclair paper
      rmax = 1.375 * body_mass^-(.315) #Sinclair, A. R., 2003. Mammal population regulation, keystone processes and ecosystem dynamics. Philos Trans R Soc Lond B Biol Sci 358:1729–40.
      if(ME.switch){
        
        # #solutions for k1 and k2 based on maximum me  Douhard et al. Frontiers in Zoology (2016) 13:32 and described on manuscript pages 14 and 15.
        
        k1 =  0.1375e-1 * body_mass ^ (-0.63e2 / 0.200e3) * exp(0.9589363630e18 * (-0.2207274913e1 + log(0.1e1 / (0.11e2 + 0.1396160558e1 * body_mass ^ (0.19e2 / 0.200e3)))) * exp(0.217147241e-1 * log(body_mass) ^ 2) / (-0.9589363630e18 * exp(0.217147241e-1 * log(body_mass) ^ 2) + 0.5000000000e18 * body_mass ^ (0.415000001e9 / 0.5000000000e10)))

      k2 = -(-0.2207274913e1 + log(0.1e1 / (0.11e2 + 0.1396160558e1 * body_mass ^ (0.19e2 / 0.200e3)))) * body_mass ^ (0.415000001e9 / 0.5000000000e10) / (-0.1917872726e1 * exp(0.217147241e-1 * log(body_mass) ^ 2) + body_mass ^ (0.415000001e9 / 0.5000000000e10))

      
      }
      else{
        k1 = 0.1375e-1 * body_mass ^ (-0.63e2 / 0.200e3) * exp(-0.31104e1 * (-0.2207274913e1 + log(0.1e1 / (0.11e2 + 0.1396160558e1 * body_mass ^ (0.19e2 / 0.200e3)))) * body_mass ^ (0.79e2 / 0.100e3) / (0.3110400000e1 * body_mass ^ (0.79e2 / 0.100e3) - 0.1400000000e0 * body_mass - 0.1600000000e0 * body_mass ^ (0.3e1 / 0.4e1)))
        
        
        k2 =  -(0.35e1 * body_mass + 0.4e1 * body_mass ^ (0.3e1 / 0.4e1)) * (-0.2207274913e1 + log(0.1e1 / (0.11e2 + 0.1396160558e1 * body_mass ^ (0.19e2 / 0.200e3)))) / (-0.7776e2 * body_mass ^ (0.79e2 / 0.100e3) + 0.35e1 * body_mass + 0.4e1 * body_mass ^ (0.3e1 / 0.4e1))
        
      }
      r <- rmax - k1 * exp(k2 * q) # primary result for manuscript
      dNdt <- r * N #not used in manucript
      variance.fx = var(fx)
      mean.fx = mean(fx)
      k <- moments::kurtosis(fx)
      gamma <- moments::skewness(fx)
      bimod = (gamma + 1)/k #Pearson, K (1916). "Mathematical contributions to the theory of evolution, XIX: Second supplement to a memoir on skew variation". Philosophical Transactions of the Royal Society A. 216 (538–548): 429–457
      r_stats = c(N, N.days.km, z, mu.eat, r, dNdt, mean.fx, variance.fx, k, gamma, bimod, n.species)
      names(r_stats) = NULL
    names(r_stats) = c("N/km^2*365days", "N.days/km^2", "z", "mu.eat", "r", "dNdt","mean.fx", "variance.fx", "k", "gamma", "bimod", "n.species")
    return(r_stats)
}# end of function
