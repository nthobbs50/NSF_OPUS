#This script creates Figures 7 and 8 in the manuscript.  It requires list d created by the function get_dist and the matrix M created in the main workflow.  
dist_rate_panel_plot = function(d, M, d.indices, tilde.K, eta = .02145){
  library(cowplot)
  library(tidyverse)
  p1 = list()
  p2 = list()
tags = c("A", "B", "C")
for(j in 1:length(d.indices)){
  i = d.indices[j]
  nut_dist = get_dist(mu.b.all = d[[i]]$mu.b.all, sigma.b.all = d[[i]]$sigma.b.all, mu.e.all = d[[i]]$mu.e.all, sigma.e.all = d[[i]]$sigma.e.all, p = d[[i]]$p, tilde.K = tilde.K, MEflag = TRUE, Hobbs_switch = d[[i]]$Hobbs_switch)
print(unique(d[[i]]$Veg_type))
    ME = as_tibble(nut_dist$me.dist)
    ME$rep = rownames(ME)
    rep.sample = sample(as.integer(ME$rep),300)
    ME = ME[rep.sample, ]
    names(ME)[ncol(ME)] = "rep"
    ME = ME %>% pivot_longer(!rep)
    ME$color_line = "grey40"
    color_line = ME$color_line
    
    p1[[j]] = ggplot() + geom_density(data = ME, (aes(x=value,color = rep)), size = .05) + theme_classic() + theme(legend.position = "none") + scale_color_manual(values = color_line) + labs( x = "Metabolizable energy (kJ / g dm)", y = "Probability density") + labs(tag = tags[j]) + ylim(0, .75) + xlim(0,18) + theme(plot.margin=unit(c(.5,1,1,1), 'cm')) 
   
} #end of j    
M_means1 = M %>% filter(Veg_type == d[[d.indices[1]]]$Veg_type) %>% group_by(Veg_type, phi, body_mass) %>% summarize(mean_r = mean(r), lower_r = quantile(r, .025), upper_r = quantile(r, .975), mean_biomass = mean(total.biomass), mean_dNdt = mean(dNdt), lower_dNdt = quantile(dNdt, .025), upper_dNdt = quantile(dNdt, .975)) %>% ungroup()
mean_biomass[j] = unique(M_means[,,j]$mean_biomass)
body_mass = M_means1$body_mass
print(M_means1$mean_biomass * phi / (eta*body_mass)) 
top.label = latex2exp::TeX(r"(\textrm{Population density (animal days } ha$^{-1}$))")
p2[[1]] = M_means1 %>% ggplot(aes(x=phi, y = mean_r)) + geom_line()+ scale_x_continuous(sec.axis = sec_axis(~. * (unique(M_means1$mean_biomass) / (eta * body_mass)), name = top.label)) + geom_line(aes(x = phi, y = upper_r), linetype = "dashed") + geom_line(aes(x = phi, y = lower_r), linetype = "dashed") + labs(x = expression(paste("Proportion of biomass consumed ", phi)), y = expression(paste("Per-capita rate of increase ", italic(r)," (", yr^-1, ")"))) + ylim(-.10, .3) + theme_classic()  + theme(plot.margin=unit(c(.5,1,1,1), 'cm'))

M_means2 = M %>% filter(Veg_type == d[[d.indices[2]]]$Veg_type) %>% group_by(Veg_type, phi, body_mass) %>% summarize(mean_r = mean(r), lower_r = quantile(r, .025), upper_r = quantile(r, .975), mean_biomass = mean(total.biomass), mean_dNdt = mean(dNdt), lower_dNdt = quantile(dNdt, .025), upper_dNdt = quantile(dNdt, .975)) %>% ungroup()
mean_biomass[j] = unique(M_means[,,j]$mean_biomass)
body_mass = M_means2$body_mass
print(M_means2$mean_biomass * phi / (eta*body_mass)) 
top.label = latex2exp::TeX(r"(\textrm{Population density (animal days } ha$^{-1}$))")
p2[[2]] = M_means2 %>% ggplot(aes(x=phi, y = mean_r)) + geom_line()+ scale_x_continuous(sec.axis = sec_axis(~. * (unique(M_means2$mean_biomass) / (eta * body_mass)), name = top.label)) + geom_line(aes(x = phi, y = upper_r), linetype = "dashed") + geom_line(aes(x = phi, y = lower_r), linetype = "dashed") + labs(x = expression(paste("Proportion of biomass consumed ", phi)), y = expression(paste("Per-capita rate of increase ", italic(r)," (", yr^-1, ")"))) + ylim(-.10, .3) + theme_classic()  + theme(plot.margin=unit(c(.5,1,1,1), 'cm'))

M_means3 = M %>% filter(Veg_type == d[[d.indices[3]]]$Veg_type) %>% group_by(Veg_type, phi, body_mass) %>% summarize(mean_r = mean(r), lower_r = quantile(r, .025), upper_r = quantile(r, .975), mean_biomass = mean(total.biomass), mean_dNdt = mean(dNdt), lower_dNdt = quantile(dNdt, .025), upper_dNdt = quantile(dNdt, .975)) %>% ungroup()
mean_biomass[j] = unique(M_means[,,j]$mean_biomass)
body_mass = M_means3$body_mass
print(M_means3$mean_biomass * phi / (eta*body_mass)) 
top.label = latex2exp::TeX(r"(\textrm{Population density (animal days } ha$^{-1}$))")
p2[[3]] = M_means3 %>% ggplot(aes(x=phi, y = mean_r)) + geom_line()+ scale_x_continuous(sec.axis = sec_axis(~. * (unique(M_means3$mean_biomass) / (eta * body_mass)), name = top.label)) + geom_line(aes(x = phi, y = upper_r), linetype = "dashed") + geom_line(aes(x = phi, y = lower_r), linetype = "dashed") + labs(x = expression(paste("Proportion of biomass consumed ", phi)), y = expression(paste("Per-capita rate of increase ", italic(r)," (", yr^-1, ")"))) + ylim(-.10, .3) + theme_classic()  + theme(plot.margin=unit(c(.5,1,1,1), 'cm'))






p_grid = plot_grid(p1[[1]], p2[[1]], p1[[2]], p2[[2]], p1[[3]], p2[[3]], nrow = 3, scale = 1.05)
return(p_grid)
} #end of function
load("All_types_d_ME.Rdata")
load("M_matrix_r_estimate_from_ME.Rdata")
p_grid1 = dist_rate_panel_plot(d = d, M = M, d.indices = c(26, 21, 23), tilde.K = 1000)  #26 = spring grazed bluebunch wheatgrass, 21 = Doug fir hemlock, 23 = Ponderosa pine
pdf(file = "../Manuscript/figures/Distributions_and_r1.pdf", width = 7, height = 10)
p_grid1
dev.off()

p_grid2 = dist_rate_panel_plot(d = d, M = M, d.indices = c(1,8,16), tilde.K = 1000)
pdf(file = "../Manuscript/figures/Distributions_and_r2.pdf", width = 8 , height = 11)
p_grid2
dev.off()

