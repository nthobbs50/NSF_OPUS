#This script executes all of the code needed to reproduce results in the manuscript "A general, resource-based explanation for density dependence in populations of large herbivores" by N.T. Hobbs


##Get functions needed for workflow being careful to comment out testing lines
source("fit_moments_function.R") #find biomass and quality moments from field data
source("get_dist_functionIII.R") ##approximate distributions from moments (also contains code for plotting)
source("r_and_stats_from_distsV2.R") # function to compute population growth parameters from nutrient distributions and summarize statistics for those distributions


##get the data
load("y_all.Rdata")  #data will be stored as .csv file in permanent, publicly accessible repository when the manuscript is accepted. 

######ME analysis#######################
########################################
##Provides plots of statistics of distributions and their effect on the strength of density dependence

##fit moments
fit_moments(ME.switch = TRUE, y.all = y.all, n.adapt = 1000, n.iter = 50000, n.update = 5000, n.sample.chains = 1000, n.thin = 1)

load("Spalinger_d_ME.Rdata")

#plot_dists(d = d, MEflag = TRUE, Hobbs_switch = FALSE)
#dev.off()
d.Spalinger = d
for(i in 1:length(d.Spalinger)) d.Spalinger[[i]]$Hobbs_switch = FALSE

#pdf(file = "ME_distribution_plots_Shipley.pdf")
load("Shipley_d_ME.Rdata")
#plot_dists(d = d, MEflag = TRUE, Hobbs_switch = FALSE)
#dev.off()
d.Shipley = d
for(i in 1:length(d.Shipley)) d.Shipley[[i]]$Hobbs_switch = FALSE

l.Spal = length(d.Spalinger)
l.Ship = length(d.Shipley)



##Combine lists
for(i in 1:sum(length(d.Spalinger), length(d.Shipley))){
  if(i <= l.Spal) d[[i]] = d.Spalinger[[i]]
  if(i > l.Spal & i <= (l.Spal + l.Ship)) d[[i]] = d.Shipley[[i - l.Spal]]
 }
length(d)
for(i in 1:length(d)) print(c(i, d[[i]]$n.species))
for(i in 1:length(d)) print(c(i, d[[i]]$Veg_type))
v = numeric()
for(i in 1:length(d)){
 v[i]= (d[[i]]$Veg_type)
 print(i)
 print(v[i])
}
v[duplicated(v)]
save(d, file = "All_types_d_ME.Rdata")
class(d[[1]])
names(d[[1]])
class(d[[1]]$mu.e.all)
head(d[[1]]$mu.e.all)
nrow(d[[1]]$mu.e.all)  ##Number of columns is number of species, number of rows is reps
#####Make matrix of r vs phi and distribution statistics
phi = seq(.01,1,.01)
n.types = length(d)
M = as.data.frame(matrix(nrow = n.types * 100 * length(phi), ncol = 17))
names_r_stats = c("N/km^2*365days", "N.days/km^2", "z", "mu.eat", "r", "dNdt","mean.fx", "variance.fx", "k", "gamma", "bimod", "n.species")
names(M)=c("Veg_type", "rep", "phi", "body_mass", "total.biomass", names_r_stats)
#d is a list of length n types that contains chains for means and variances of each plant species and identifying information. It should not be confused with the vector d described in Appnedix S1
ME.switch = TRUE
counter = 1


for(i in  1:n.types){
  pdf_file = paste0("dist_plots/",unique(d[[i]]$Veg_type),".pdf")
  pdf(file = pdf_file)
  print(i)
  start.time = Sys.time()
  out = get_dist(mu.b.all = d[[i]]$mu.b.all, sigma.b.all = d[[i]]$sigma.b.all, mu.e.all = d[[i]]$mu.e.all, sigma.e.all = d[[i]]$sigma.e.all, p = d[[i]]$p, tilde.K = 1000, MEflag = ME.switch, Hobbs_switch = d[[i]]$Hobbs_switch)
  dist = out$me.dist
  total.biomass = out$total.biomass
  n.species = out$n.species
  n.rep = 100
  if(is.null(n.species)) next
  for(j in 1 : n.rep){
    for(k in 1:length(phi)){
      r_stats = get_r(fx=dist[j,], phi = phi[k] , body_mass = unique(d[[i]]$body_mass), total.biomass = total.biomass[j], n.species = n.species, ME.switch = ME.switch) 
      
      M[counter,1] =  unique(d[[i]]$Veg_type)
      M[counter,2] = j  #reps 
      M[counter,3] = phi[k] #phi
      M[counter,4] = unique(d[[i]]$body_mass)
      M[counter,5] = out$total.biomass[j]
      M[counter, 6:17] = t(r_stats[1:length(r_stats)]) #t needed to create row of columns
     counter = counter + 1
    } #end of j
  }# end of k
  par(mfrow = c(1,2))
  plot(density(out$me.dist[1,]), xlab =  "Metabolizable energy (kjoule / gm", ylab = "Probability density", main =    d[[i]]$Veg_typ, lwd = .5, col = "grey", ylim = c(0,.25))
  
  for(j in 2:100) lines(density(out$me.dist[j,]),  col = adjustcolor("grey", alpha = 0.8), lwd = .5 )
 
 M_means = M %>% filter(Veg_type ==d[[i]]$Veg_type) %>% group_by(phi) %>% summarize(mean_r = mean(r), lower = quantile(r, .025, na.rm = TRUE), upper = quantile(r, .95, na.rm = TRUE)) 
  
  plot(M_means$phi,M_means$mean_r, type = "l", xlab = expression(paste("Proportion of bioimass consumed ", phi)), ylab = expression(paste("Per-capita rate of increase ", italic(r))), ylim = c(-.05, .3))
  lines(M_means$phi, M_means$upper, lty = "dashed")
  lines(M_means$phi, M_means$lower, lty = "dashed")
  print(Sys.time() - start.time)
  dev.off()
} #end of i

save(M, file = "M_matrix_r_estimate_from_ME.Rdata")

#Make plots of analytical results, figures 3 - 6 in manuscript
source("analytical_r_as_f(phi)_FINAL.R")

#Make plots of derivative of r vs phi for different values of  mean and variance of nutrient distribution, figures 9 and 10 in the manuscript. 
source("d_r_d_phi_plot.R")
   
#Make plots of approximations of nutrient distributions and r vs phi curve, figures 7 and 8 in manuscript
source("two_x_axis_plotv3.R")






