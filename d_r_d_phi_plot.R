#This script creates Figures 9 and 10 in the manuscript and associated statistics. The R data object M_matrix_r_estimate_from_ME.Rdata is created in the main workflow.

load("M_matrix_r_estimate_from_ME.Rdata")
library(rjags)
library(MCMCvis)
#======Model fit


temp.M = filter(M, n.species >= 5)
unique(temp.M$n.species)
dr_dphi = temp.M %>% group_by(Veg_type, phi)%>% 
  summarize(
    mean.r = mean(r),
    var.fx = mean(variance.fx),
    mean.fx = mean(mean.fx),
    mean.k = mean(k),
    mean.gamma = mean(gamma),
    mean.bimod = mean(bimod)
  ) %>% 
  mutate(dr.dphi = (lead(mean.r)-lag(mean.r)) / (lag(phi) + lead(phi))) %>% filter(round(phi,2) == .1 | phi == .25 | phi == .50) %>%
  mutate(phi.index = as.numeric(as.factor(as.character(phi)))) %>% ungroup
unique(dr_dphi$phi)
unique(dr_dphi$phi.index)

d_r_d_phi_plot = function(dr_dphi, predictor, x.label, xlim, ylim, new.x, color_var, leg_title){
  
 { #extra brackets only needed for R markdown
  sink("dr_dphi.R")
  cat("
  model{
  #priors
  for(i in 1:3){
    B0[i] ~ dnorm(0,.0001)
    B1[i] ~ dnorm(0, .0001)
    sigma[i] ~ dunif(0,2)
    tau[i] = 1/sigma[i]^2
  }
  #likelihood
  for(i in 1:length(y)){
    mu[i] = B0[phi.index[i]] + B1[phi.index[i]] * x[i]
    y[i] ~ dnorm(mu[i], tau[phi.index[i]])
  }
  #plotting
  for(j in 1:3){
    for(k in 1:length(new.x)){
      mu.hat[j,k] = B0[j] + B1[j] * new.x[k]
    }
  }
  
  } #end of model
      ",fill=TRUE)
  sink()

}

jags_x = pull(dr_dphi, {{predictor}})

j.data = list(y = dr_dphi$dr.dphi, x = jags_x, phi.index = dr_dphi$phi.index, new.x = new.x)


jm = jags.model("dr_dphi.R", n.chains = 3, data = j.data)
update(jm, n.iter = 10000)
zc = coda.samples(jm, n.iter = 50000, variable.names = c("B0","B1","sigma", "mu.hat"))

MCMCsummary(zc, excl = "mu.hat")

pred = MCMCpstr(zc, params = "mu.hat", func = function(x) quantile(x, probs = c(.025,.5,.975)))

phi.10 = as.data.frame(cbind(j.data$new.x, pred$mu.hat[1,,]))
names(phi.10)[1] = "new.x"
phi.10$phi = ".10"
temp = filter(dr_dphi, round(phi,2) == .1)
ggplot() + geom_point(data = temp, aes( x = var.fx, y = dr.dphi)) + geom_line(data = phi.10, aes( x = new.x, y = `50%`))


phi.25 = as.data.frame(cbind(j.data$new.x, pred$mu.hat[2,,]))
names(phi.25)[1] = "new.x"
phi.25$phi = ".25"
temp = filter(dr_dphi, phi == .25)
ggplot() + geom_point(data = temp, aes( x = {{predictor}}, y = dr.dphi)) + geom_line(data = phi.25, aes( x = new.x, y = `50%`))


phi.50 = as.data.frame(cbind(j.data$new.x, pred$mu.hat[3,,]))
names(phi.50)[1] = "new.x"
phi.50$phi = ".50"
temp = filter(dr_dphi, phi == .50)
ggplot() + geom_point(data = temp, aes( x = {{predictor}}, y = dr.dphi)) + geom_line(data = phi.50, aes( x = new.x, y = `50%`))



library(latex2exp)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
p1.points = filter(dr_dphi, round(phi,2) == .1 )
p1 = ggplot() + geom_point(data = p1.points, aes(x = {{predictor}}, y = dr.dphi, color = {{color_var}})) + labs(y = NULL, x = NULL, color = leg_title, title = expression(paste("A    ",phi, " = .1"))) + xlim(xlim)  + ylim(ylim) + geom_line(data = phi.10, aes(x = new.x, y = `50%`) ) + geom_line(data = phi.10, aes(x = new.x, y = `97.5%`), linetype = "dashed" ) + geom_line(data = phi.10, aes(x = new.x, y = `2.5%`), linetype = "dashed" )

p2.points = filter(dr_dphi, phi == .25 )
p2 = ggplot() + geom_point(data = p2.points, aes(x = {{predictor}}, y = dr.dphi, color = {{color_var}})) + labs(y = NULL, x = NULL, color = leg_title, title = expression(paste("B    ",phi, " = .25"))) + xlim(xlim)  + ylim(ylim) + geom_line(data = phi.25, aes(x = new.x, y = `50%`) ) + geom_line(data = phi.25, aes(x = new.x, y = `97.5%`), linetype = "dashed" ) + geom_line(data = phi.25, aes(x = new.x, y = `2.5%`), linetype = "dashed" ) + theme_classic()+ theme(legend.position="none") 

p3.points = filter(dr_dphi, phi == .50 )
p3 = ggplot() + geom_point(data = p3.points, aes(x = {{predictor}}, y = dr.dphi, color = mean.fx)) + labs(y = NULL, x = NULL, color = leg_title, title = expression(paste("C    ",phi, " = .50"))) + xlim(xlim)  + ylim(ylim) + geom_line(data = phi.50, aes(x = new.x, y = `50%`) ) + geom_line(data = phi.50, aes(x = new.x, y = `97.5%`), linetype = "dashed" ) + geom_line(data = phi.50, aes(x = new.x, y = `2.5%`), linetype = "dashed" ) + theme_classic()+ theme(legend.position="none") 

legend = get_legend(p1)
p1 <- p1 + theme_classic() + theme(legend.position="none") 

p_all = gridExtra::grid.arrange(p1,p2,p3, legend, nrow = 1, bottom = x.label, left = "Derivative of per capita growth rate dr/d\u03C6", widths = c(2.25, 2.25, 2.25, 1))

name_part = colnames(select(dr_dphi, {{predictor}}))
ggsave(p_all, filename = paste("../Manuscript/figures/Derivative_",name_part, "_collage.pdf"), device = cairo_pdf, width = 8, height = 4, units = "in")
}


d_r_d_phi_plot(dr_dphi = dr_dphi, predictor = var.fx, x.label = "Variance of metabolizable energy (ME) distribution", xlim = c(0,8), ylim = c(-.08,0), new.x = seq(0,8,.1), color_var = mean.fx, leg_title = "Mean ME")

d_r_d_phi_plot(dr_dphi = dr_dphi, predictor = mean.fx, x.label = "Mean of metabolizable energy (ME) distribution", xlim = c(7,12), ylim = c(-.08,0), new.x = seq(7,12,.1), color_var = var.fx, leg_title = "Variance ME")



