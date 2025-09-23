# Ship Application

library(tidyverse)
library(stats)
library(qqplotr)
library(evd)


##################################
#### Function for Gumbel JML  ####
##################################
# (sourced from Gumbel_jointEstimation.R)

fA <- function(u, v, r) exp(-u / r) + exp(-v / r)
dAdm <- function(u, r, sig1) exp(-u / r) / (r * sig1)

I_int <- function(u, v, r, sig1, case = "mm") {
  A <- fA(u, v, r)
  B <- dAdm(u, r, sig1)
  D <- B * u
  G <- -B / sig1 + D / (r * sig1)
  
  if (case == "mm") {
    return(-(r * (r - 1) * A^(r - 2) * B^2 + A^(r - 1) * B / sig1) *
             (1 - 1 / (A^r - (r - 1) / r)) +
             (r - 2) * (B / (r * sig1 * A) - B^2 / A^2) -
             (r * A^(r - 1) * B)^2 / (1 / r - 1 + A^r)^2)
  } else if (case == "ss") {
    return(-(r * (r - 1) * A^(r - 2) * D^2 + (u - 2 * r) * A^(r - 1) * D / sig1) *
             (1 - 1 / (A^r - (r - 1) / r)) +
             (r - 2) * (D * (u - 2 * r) / (r * sig1 * A) - D^2 / A^2) +
             (r - 2 * u) / (r * sig1^2) -
             (r * A^(r - 1) * D)^2 / (1 / r - 1 + A^r)^2)
  } else if (case == "ms") {
    return(-(r * (r - 1) * A^(r - 2) * D * B + r * A^(r - 1) * G) *
             (1 - 1 / (A^r - (r - 1) / r)) +
             (r - 2) * (G / A - B * D / A^2) -
             1 / (r * sig1^2) -
             (r * A^(r - 1))^2 * B * D / (1 / r - 1 + A^r)^2)
  }
}

integrate_inner <- function(v, r, sig1, case) {
  int_u <- function(u) {
    I_int(u, v, r, sig1, case) *
      dbvevd(c(u, v), dep = r, model = "log", mar1 = c(0, 1, 0), mar2 = c(0, 1, 0))
  }
  integrate(Vectorize(int_u), -10, 25)$value
}
integrate_inner <- Vectorize(integrate_inner)

integrate_outer <- function(r, sig1, case = "mm") {
  integrate(function(v) integrate_inner(v, r, sig1, case), -10, 25)$value
}

Var_joint <- Vectorize(function(r, sig1) {
  I_mm <- -integrate_outer(r, sig1, "mm")
  I_ms <- -integrate_outer(r, sig1, "ms")
  I_ss <- -integrate_outer(r, sig1, "ss")
  solve(matrix(c(I_mm, I_ms, I_ms, I_ss), 2))
})

#############################
#### Data Visualization  ####
#############################
setwd("./pMFMC/scripts")
# --- Read LAMP and SC values: record maxima ---
SCextreme = as_tibble(readRDS("../data/SCextreme.rds"))
LAMPextreme = as_tibble(readRDS("../data/LAMPextreme.rds"))

# --- plot for one record ---
seed = 25001
LAMPone =  read.table(sprintf("../data/L2_ONRFL_%05d.mot",seed),skip = 2, header = T)
LAMPone = LAMPone[-1,]
SCone =  read.table(sprintf("../data/flh_irreg4a-th-%07d.mot",seed),skip = 2, header = T)
colnames(SCone) = c("Time","Xcg","Ycg","Zcg","Roll","Pitch","Yaw")

p = data.frame(time = LAMPone$Time, LAMP = LAMPone$Zcg, SC = SCone$Zcg) %>%
  pivot_longer(c(LAMP,SC), values_to = "heave", names_to = "type") %>%
  filter(time >= 50, time <= 150, heave >= -5, heave <= 8) %>%
  ggplot() + geom_line(aes(time,heave, linetype=type)) + 
  xlim(50,150) +
  ylim(-5,8) + 
  ylab("Heave (m)") +
  xlab("time (s)") + 
  theme(legend.position = c(0.25, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17)) 
p
ggsave(filename = "../images/lamp_sc_heave.png", plot = p, width = 5, height = 4.5)



# ---- QQplot for one record ----
LAMPone =data.frame(heave = LAMPone$Zcg[seq(1,18000,20)])
SCone =data.frame(heave = SCone$Zcg[seq(1,17999,20)])

p=LAMPone %>%
  ggplot( mapping = aes(sample = heave)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  ggtitle("LAMP Heave")+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17))
p
ggsave(filename = "../images/qqplot_heave_lamp_one.png", plot = p, width = 5, height = 4.5)

p=SCone %>%
  ggplot( mapping = aes(sample = heave)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  ggtitle("SC Heave")+
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17))
p
ggsave(filename = "../images/qqplot_heave_sc_one.png", plot = p, width = 5, height = 4.5)


# ---- scatterplot for joint observation ----
SCjoint = SCextreme %>% filter(25000 < seed & seed <= 25100)
p=data.frame(seed = SCjoint$seed, LAMP = LAMPextreme$max, SC = SCjoint$max) %>%
  ggplot() + geom_point(aes(SC,LAMP)) + 
  ylab("LAMP heave record maxima (m)") +
  xlab("SC heave record maxima (m)") + 
  #ylim(5,8)+
  #xlim(5,8.5)+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", linewidth=0.5)+
  theme(legend.position = c(0.225, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17))
p
ggsave(filename = "../images/scatter_heave.png", plot = p, width = 5, height = 4.5)


#### Fit parametric distribution (Gumbel) ####

dgumbel <- function(x, loc, scale) {
  z <- (x - loc) / scale
  lp = - log(scale) -(z + exp(-z))
  return(exp(lp))
}

pgumbel <- function(q, loc, scale) {
  lp = -exp(-(q - loc) / scale)
  exp(lp)
}

qgumbel <- function(p, loc, scale) {
  loc - scale * log(-log(p))
}
gamma = -digamma(1)
fit.lofi = fgev(SCextreme$max, shape=0)


p=SCextreme %>% ggplot(aes(x=max))+
  geom_histogram(aes(y = after_stat(density)),
                 colour = "darkgray", fill = "gray") +
  stat_function(fun= function(x){dgumbel(x, loc=fit.lofi$estimate[1], scale = fit.lofi$estimate[2])}) +
  ggtitle("Histogram and PDF of SC heave record maxima") +
  labs(x = "record maxima (m)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17))
p
ggsave(filename = "../images/density_heave_sc.png", plot = p, width = 5, height = 4.5)


n <- nrow(SCextreme)
probs <- (1:n) / (n + 1)
theoretical_quantiles <- qgumbel(probs, loc=fit.lofi$estimate[1], scale = fit.lofi$estimate[2] )

# Create a data frame with sample and theoretical quantiles
qq_data <- data.frame(
  sample = sort(SCextreme$max),
  theoretical = theoretical_quantiles
)

di = "gumbel"
dp = list(loc = fit.lofi$estimate[1], scale = fit.lofi$estimate[2] )
gg <- ggplot(data = qq_data, mapping = aes(sample = sample)) +
    stat_qq_band(distribution = di, dparams = dp, bandType="boot") +
    stat_qq_line(distribution = di, dparams = dp) +
    stat_qq_point(distribution = di, dparams = dp) +
    ggtitle("SC Heave record maxima") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme(plot.title = element_text(hjust = 0.5, size = 17),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17))
gg
ggsave(filename = "../images/qqplot_gumbel_heave_sc.png", plot = gg, width = 5, height = 4.5)


##########################
#### MFMC estimation  ####
##########################

# moment formulation
dgdy_mom <- function(fitted){
  m1 = fitted[1]
  m2 = fitted[2]
  gamma = -digamma(1)
  g11 = 1 + gamma * sqrt(6/pi^2)* m1/sqrt(m2-m1^2)
  g12 = -gamma * sqrt(6/pi^2) / (2 *sqrt(m2-m1^2))
  g21 = -sqrt(6/pi^2) * m1/sqrt(m2-m1^2)
  g22 = sqrt(6/pi^2)/ (2 *sqrt(m2-m1^2))
  return(matrix(c(g11, g12, g21, g22),nrow=2, byrow=T))
}
g_mom <- function(y){
  y1= y[1]
  y2 = y[2]
  gamma = -digamma(1)
  return(c(y1 - gamma * sqrt(6/pi^2 *(y2-y1^2)), sqrt(6/pi^2 *(y2-y1^2))))
}

data.joint = data.frame(seed = SCjoint$seed, LAMP = LAMPextreme$max, SC = SCjoint$max)

#### Baselines ####
# bl-ml
fit.hifi = fgev(data.joint$LAMP, shape=0)
fit.hifi$var.cov
gamma = -digamma(1)
bl.ml.est = fit.hifi$estimate
bl.ml.var= fit.hifi$estimate[2]^2 * 6/pi^2 * matrix(c((gamma-1)^2+pi^2/6, 1-gamma, 1-gamma, 1),nrow=2)
nrow(data.joint) *diag(fit.hifi$var.cov)
pMFMC = data.frame(est = c(fit.hifi$estimate,NA), var=as.vector(bl.ml.var)[c(1,4,2)], type=rep("bl-ml",3), pType= c("mu","sigma","cov"))

# bl-mom
muY.est = c(mean(data.joint$LAMP), mean(data.joint$LAMP^2))
G = dgdy_mom(muY.est)
bl.mom.var = G %*% cov(data.frame(y1=data.joint$LAMP, y2= data.joint$LAMP^2)) %*% t(G)
pMFMC = rbind(pMFMC, data.frame(est = c(g_mom(muY.est),NA), var=as.vector(bl.mom.var)[c(1,4,2)], type=rep("bl-mom",3), pType= c("mu","sigma","cov")))

pMFMC

#### MML ####
fit.lofi.all = fgev(SCextreme$max, shape=0)
fit.lofi.joint = fgev(data.joint$SC, shape=0)
est.lofi.all = fit.lofi.all$estimate
est.lofi.joint = fit.lofi.joint$estimate


a1 = 6/pi^2*(gamma-1)^2+1; b1= 6/pi^2*(1-gamma)
a2 = 6/pi^2*(1-gamma); b2 =  6/pi^2
varP <- function(a,b){
  var_int <- function(u, a,b) {
        (a-b-a*exp(-u)+b*u-b*u*exp(-u))^2* dgev(u)
  }
  var_int = Vectorize(var_int)
  integrate(function(x) var_int(x,a,b), lower = -10, upper = 25)$value
}

h_Y <- function(x,mu,sig,a,b){
  z = (x-mu)/sig
  sig * (a-b-a*exp(-z) + b*z-b*z*exp(-z))
}
Y1.hifi=sapply(data.joint$LAMP, function(x) h_Y(x, bl.ml.est[1],bl.ml.est[2], a1,b1))
Y2.hifi=sapply(data.joint$LAMP, function(x) h_Y(x, bl.ml.est[1],bl.ml.est[2], a2,b2))
Y1.lofi=sapply(data.joint$SC, function(x) h_Y(x, est.lofi.joint[1],est.lofi.joint[2], a1,b1))
Y2.lofi=sapply(data.joint$SC, function(x) h_Y(x, est.lofi.joint[1],est.lofi.joint[2], a2,b2))

beta_opt= c(cov(Y1.hifi, Y1.lofi)/var(Y1.lofi), cov(Y2.hifi, Y2.lofi)/var(Y2.lofi))
mf.mml.est = bl.ml.est + beta_opt * (est.lofi.all - est.lofi.joint)
mf.mml.var = c(var(Y1.hifi) * (1- cor(Y1.hifi,Y1.lofi)^2), var(Y2.hifi) * (1- cor(Y2.hifi,Y2.lofi)^2))
mf.mml.var = c(mf.mml.var, cov(Y1.hifi, Y2.hifi) - (beta_opt[1]* cov(Y2.hifi, Y1.lofi) + beta_opt[2]* cov(Y1.hifi, Y2.lofi) - beta_opt[1]* beta_opt[2] * cov(Y1.lofi, Y2.lofi)))
pMFMC = rbind(pMFMC, data.frame(est = c(mf.mml.est,NA), var=mf.mml.var, type=rep("mf-mml",3), pType= c("mu","sigma","cov")))
pMFMC

#### MoM ####
Yh1 = data.joint$LAMP; Yl1 = data.joint$SC
Yh2 = data.joint$LAMP^2; Yl2 = data.joint$SC^2
Vh1 = var(Yh1); Vl1 = var(Yl1); Vh2 = var(Yh2); Vl2 = var(Yl2);
Chl11 = cov(Yh1,Yl1); Chl21 = cov(Yh2, Yl1); Chh12 = cov(Yh1, Yh2)
Chl22 = cov(Yh2, Yl2); Chl12 = cov(Yh1, Yl2); Cll12 = cov(Yl1, Yl2); 

# --- Appendix B proposition ---
alpha_opt <-function(g1,g2){
  return(c((Vl2*(g1*Chl11 + g2*Chl21)-Cll12*(g2*Chl22 + g1*Chl12))/(g1*(Vl1*Vl2 - Cll12^2)),
         (Vl1*(g2*Chl22 + g1*Chl12) - Cll12*(g1*Chl11 + g2*Chl21))/(g2*(Vl1*Vl2 - Cll12^2)))
  )
}
muY.est = c(mean(data.joint$LAMP), mean(data.joint$LAMP^2))
G = dgdy_mom(muY.est)
alpha_opt_1 = alpha_opt(G[1,1], G[1,2])
mu_Y1 = c( mean(Yh1) + alpha_opt_1[1]*(mean(SCextreme$max) - mean(Yl1)), 
           mean(Yh2) + alpha_opt_1[2]*(mean(SCextreme$max^2) - mean(Yl2)) )
  
mu.update = g_mom(mu_Y1)[1]


muY.est = c(mean(data.joint$LAMP), mean(data.joint$LAMP^2))
G = dgdy_mom(muY.est)
alpha_opt_2 = alpha_opt(G[2,1], G[2,2])
mu_Y2 = c( mean(Yh1) + alpha_opt_2[1]*(mean(SCextreme$max) - mean(Yl1)), 
           mean(Yh2) + alpha_opt_2[2]*(mean(SCextreme$max^2) - mean(Yl2)) )
sig.update = g_mom(mu_Y2)[2]

mf.mom.est = c(mu.update, sig.update)

fit.Sig.opt <- function(par1,par2){
  VmuY1 = Vh1 - 2 * par1[1] * Chl11 + par1[1]^2 * Vl1
  CmuY12 = Chh12 - par1[1] * Chl21 - par1[2] * Chl12 + par1[1] * par1[2] * Cll12
  VmuY2 = Vh2 - 2 * par1[2] * Chl22 + par1[2]^2 * Vl2
  var.muY.1 = matrix(c(VmuY1,CmuY12,CmuY12, VmuY2),nrow=2)
  
  VmuY1 = Vh1 - 2 * par2[1] * Chl11 + par2[1]^2 * Vl1
  CmuY12 = Chh12 - par2[1] * Chl21 - par2[2] * Chl12 + par2[1] * par2[2] * Cll12
  VmuY2 = Vh2 - 2 * par2[2] * Chl22 + par2[2]^2 * Vl2
  var.muY.2 = matrix(c(VmuY1,CmuY12,CmuY12, VmuY2),nrow=2)
  
  mu_Y1 = c( mean(Yh1) + par1[1]*(mean(SCextreme$max) - mean(Yl1)), 
             mean(Yh2) + par1[2]*(mean(SCextreme$max^2) - mean(Yl2)) )
  mu_Y2 = c( mean(Yh1) + par2[1]*(mean(SCextreme$max) - mean(Yl1)), 
               mean(Yh2) + par2[2]*(mean(SCextreme$max^2) - mean(Yl2)) )
  
  G = dgdy_mom(muY.est)

  covY1 = Vh1 - (par1[1]+par2[1]) * Chl11 + par1[1]*par2[1] * Vl1
  covY12 = Chh12 - par1[1] * Chl21 - par2[2] * Chl12 + par1[1] * par2[2] * Cll12
  covY21 = Chh12 - par2[1] * Chl21 - par1[2] * Chl12 + par2[1] * par1[2] * Cll12
  covY2 = Vh2 - (par1[2]+par2[2]) * Chl22 + par1[2]*par2[2] * Vl2
  cov.muY = matrix(c(covY1,covY12,covY21, covY2),nrow=2, byrow=T)
  
  c(G[1,] %*% var.muY.1 %*% (G[1,]),
    G[2,] %*% var.muY.2 %*% (G[2,]),
    G[1,] %*% cov.muY %*% (G[2,])
    )
}


pMFMC = rbind(pMFMC, data.frame(est = c(mf.mom.est,NA), var=fit.Sig.opt(alpha_opt_1, alpha_opt_2), type=rep("mf-mom-constant",3), pType= c("mu","sigma","cov")))
pMFMC


#### JML ####

bv_gumbel_mf <- function(pars, x,y){
  # negative of ll
  # x is (2, n) matrix with n > 1 number of samples
  # pars = (scale1, scale2, alpha, loc1, loc2), where alpha = 1 is indep, alpha -> inf is de>
  x[1, ] <- (x[1, ] - pars[4])/pars[1]
  x[2, ] <- (x[2, ] - pars[5])/pars[2]
  y <- ( y - pars[5])/pars[2]
  
  A <- exp(-pars[3]* x[1,]) + exp(-pars[3]* x[2,])

  ll <-  - pars[3]*(x[1, ] + x[2, ]) + (1/pars[3]-2) * log(A) - exp(log(A)/pars[3]) + log(exp(log(A)/pars[3])  - 1 + pars[3]) - log(pars[1]) - log(pars[2])
  
  lmarg <- -log(pars[2]) - y - exp(-y)

  -sum(ll) - sum(lmarg)
}


m <- optim(c("scale1"=bl.ml.est[2], "scale2"=est.lofi.all[2], "1/r"=1, "loc1"=bl.ml.est[1], "loc2"=est.lofi.all[1]),
           bv_gumbel_mf, x =t(data.joint[-1]), y = SCextreme$max[SCextreme$seed>100],
           method = "L-BFGS-B",
           lower = c(0, 0, 1, -20, -20),
           upper = c(50, 50, 10000, 20, 20))

m$par[c(4,1)]
joint.ml.var= Var_joint(1/m$par[3], m$par[1])

pMFMC = rbind(pMFMC, data.frame(est = c(m$par[c(4,1)],NA), var=joint.ml.var[c(1,4,2)], type=rep("joint-ml",3), pType= c("mu","sigma","cov")))
pMFMC = pMFMC %>% mutate(var = var/nrow(data.joint))

# --- aggregate ---
pMFMC = pMFMC %>% filter(type!="mf-mom-opt") %>% mutate(type = case_when(grepl("mf-mom",type)~"mf-mom",.default = type))
pMFMC=pMFMC %>% mutate(type = case_when(type=="bl-ml" ~ "BL-ML", type=="bl-mom"~"BL-MoM", type=="mf-mml" ~ "MML",
                                  type=="mf-mom" ~ "MoM", type == "joint-ml" ~ "JML"))


###############
#### Plot  ####
###############

alpha=0.05
p=pMFMC %>% filter(pType=="sigma" ) %>%
  mutate(plot="estimator", 
         upper=(est+qnorm(1-alpha/2)*sqrt(var)),
         lower=(est-qnorm(1-alpha/2)*sqrt(var))) %>%
  ggplot(aes(type,(est))) +
  ylab("scale parameter") + xlab("Estimation Methods")+
  geom_point()+
  geom_errorbar(aes(ymin=(lower), ymax=upper), width=.1) +
    theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16), 
        axis.title.x = element_text(margin = margin(t = 15)))
p
ggsave(filename = "../images/var_plot_sig_wide.png", plot = p, width = 8, height = 3.5)


p = pMFMC %>% filter(pType=="mu") %>%
  mutate(plot="estimator", 
         upper=(est+qnorm(1-alpha/2)*sqrt(var)),
         lower=(est-qnorm(1-alpha/2)*sqrt(var))) %>%
  ggplot(aes(type,(est))) +
  ylab("location parameter") + xlab(" ")+
  geom_point()+
  geom_errorbar(aes(ymin=(lower), ymax=upper), width=.1)+
    theme(
        axis.text=element_text(size=13),
        axis.title=element_text(size=16))
p
ggsave(filename = "../images/var_plot_mu_wide.png", plot = p, width = 8, height = 3.5)


#######################
#### QoI analysis  ####
#######################

q_theta <- function(mu,sig,a, QoI="prob"){
  if(QoI=="prob"){
    z = (a-mu)/sig
    return(log10(1-exp(-exp(-z))))
  }
  else if(QoI=="quantile"){
    return(mu-sig*log(-log(a)))
  }
}
a_prob = 12
a_quan = 1-1e-2

dqdt <- function(mu, sig,a, QoI="prob"){
  if(QoI=="prob"){
    z = (a-mu)/sig
    return(c(exp( -exp(-z) -z -log(sig)), exp(-exp(-z) -z +log(z) -log(sig)))/(1-exp(-exp(-z)))/log(10))
  }
  else if(QoI=="quantile"){
    return(c(1, -log(-log(a))))
  }
}

QoImf =  pMFMC %>% rowwise() %>%
  mutate(upper = if (pType != "cov") est + qnorm(1 - alpha/2) * sqrt(var) else NA_real_,
         lower = if (pType != "cov") est - qnorm(1 - alpha/2) * sqrt(var) else NA_real_) %>% ungroup() %>%
  pivot_wider(values_from=c(est, var,upper,lower), names_from=pType) %>% 
  mutate(mu=est_mu, sig=est_sigma) %>% dplyr::select(-c(est_mu,est_sigma, est_cov, upper_cov, lower_cov)) %>%
  group_by(type) %>%
  mutate(prob = pmap_dbl(list(mu, sig), ~ q_theta(..1, ..2, a=a_prob, QoI = "prob")),
         prob_var = pmap_dbl(list(mu, sig,var_mu, var_sigma, var_cov), ~ dqdt(..1, ..2, a=a_prob, QoI = "prob") %*% matrix(c(..3,..5,..5,..4),nrow=2) %*% dqdt(..1, ..2, a=a_prob, QoI = "prob")),
         quan = pmap_dbl(list(mu, sig), ~ q_theta(..1, ..2, a_quan, QoI = "quantile")),
         quan_var = pmap_dbl(list(mu, sig,var_mu, var_sigma, var_cov), ~ dqdt(..1, ..2, a=a_quan, QoI = "quantile") %*% matrix(c(..3,..5,..5,..4),nrow=2) %*% dqdt(..1, ..2, a=a_quan, QoI = "quantile"))
         ) %>%
  ungroup()
p=QoImf %>% pivot_longer(c(prob,quan), names_to = "QoI", values_to = "estimate") %>%
  mutate(var=case_when(QoI=="prob"~prob_var, QoI=="quan"~quan_var)) %>%
  dplyr::select(-c(prob_var,quan_var)) %>%
  mutate(upper=estimate +qnorm(1-alpha/2)*sqrt(var), 
         lower= estimate-qnorm(1-alpha/2)*sqrt(var)) %>%
  filter(QoI=="prob" ) %>% #filter( type !="mf-mom"& type !="bl-mom") %>%
  ggplot(aes(type,(estimate))) +
  ylab("Log10 of Exceedance Probability") + xlab(" ")+
  geom_point()+
  geom_errorbar(aes(ymin=(lower), ymax=(upper)), width=.1)+
    theme(
        axis.text=element_text(size=13),
        axis.title=element_text(size=15))

p
ggsave(filename = "../images/exceedP_ship_12.png", plot = p, width = 8, height = 3.5)


p = QoImf %>% pivot_longer(c(prob,quan), names_to = "QoI", values_to = "estimate") %>%
  mutate(var=case_when(QoI=="prob"~prob_var, QoI=="quan"~quan_var)) %>%
  dplyr::select(-c(prob_var,quan_var)) %>%
  mutate(upper=estimate +qnorm(1-alpha/2)*sqrt(var), 
         lower= estimate-qnorm(1-alpha/2)*sqrt(var)) %>%
  filter(QoI=="quan") %>% #filter( type !="mf-mom"&type !="bl-mom") %>%
  ggplot(aes(type,(estimate))) +
  ylab("Estimate of 99% Quantile") + xlab("Estimation Methods")+
  geom_point()+
  #ylim(6.4, 7.1)+
  geom_errorbar(aes(ymin=(lower), ymax=upper), width=.1)+
    theme(
        axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        axis.title.x = element_text(margin = margin(t = 15)))

p
ggsave(filename = "../images/quantile_ship.png", plot = p, width = 8, height = 3.5)



