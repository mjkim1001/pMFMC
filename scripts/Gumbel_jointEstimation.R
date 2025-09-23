# Bivariate Gumbel Joint Estimation Analysis
# Estimating both mu and sigma

library(tidyverse)
library(evd)

# Fixed parameters
sig1 <- 4
sig2 <- 1
mu1 <- 2
mu2 <- 2
gamma <- -digamma(1)
r_seq <- seq(0.1, 1, by = 0.05)


#### Baseline MLE ####

asymp.var <- solve(matrix(c(1, gamma - 1,
                            gamma - 1, (gamma - 1)^2 + pi^2 / 6), 2)) * sig1^2

#### MML ####
a1 <- 6 / pi^2 * (gamma - 1)^2 + 1
b1 <- 6 / pi^2 * (1 - gamma)
a2 <- 6 / pi^2 * (1 - gamma)
b2 <- 6 / pi^2


varP <- function(a, b){
  var_int <- function(u, a, b) {
    (a - b - a * exp(-u) + b * u - b * u * exp(-u))^2 * dgev(u)
  }
  var_int = Vectorize(var_int)
  integrate(function(x) var_int(x,a,b), lower = -10, upper = 25)$value
}

integrate_cov_mle_inner <- function(v, r, a, b) {
  cov_fun <- function(u) {
    h1 <- a - b - a * exp(-u) + b * u - b * u * exp(-u)
    h2 <- a - b - a * exp(-v) + b * v - b * v * exp(-v)
    h1 * h2 * dbvevd(c(u, v), dep = r, model = "log",
                     mar1 = c(0, 1, 0), mar2 = c(0, 1, 0))
  }
  integrate(Vectorize(cov_fun), -10, 25)$value
}
integrate_cov_mle_inner=Vectorize(integrate_cov_mle_inner)
integrate_cov_mle <- function(r,a,b) {
  integrate(function(v) {
    integrate_cov_mle_inner(v,r,a,b) 
  }, lower = -10, upper = 25)$value
}
integrate_cov_mle = Vectorize(integrate_cov_mle)



#### JML ####
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

#### MoM ####
g11 <- 1 + gamma * (6 / pi^2) * (gamma + mu1 / sig1)
g12 <- -gamma * 3 / pi^2 / sig1
g21 <- -6 / pi^2 * (gamma + mu1 / sig1)
g22 <- 3 / pi^2 / sig1

mom_marg <- function(mu,sig,d){
  mom <- function(x) {
    x^d * dgev(x, loc = mu, scale = sig)
  }
  mom = Vectorize(mom)
  integrate(mom, lower = -10 * sig + mu, upper = 25 * sig + mu)$value
}

V1 <- pi^2 / 6 * sig1^2
V2 <- pi^2 / 6 * sig2^2
W1 <- mom_marg(mu1, sig1, 4) - mom_marg(mu1, sig1, 2)^2
W2 <- mom_marg(mu2, sig2, 4) - mom_marg(mu2, sig2, 2)^2
Cvw1 <- mom_marg(mu1, sig1, 3) - mom_marg(mu1, sig1, 2) * (mu1 + sig1 * gamma)
Cvw2 <- mom_marg(mu2, sig2, 3) - mom_marg(mu2, sig2, 2) * (mu2 + sig2 * gamma)

mom_joint_inner <- function(x2, mu1,r, d1,d2) {
  cov_int <- function(x1) {
    input_vector <- c(x1, x2)
    x1^d1*x2^d2 * dbvevd(input_vector, dep = r, model = "log", mar1 = c(mu1, sig1, 0), mar2 = c(mu2, sig2, 0))
  }
  integrate(Vectorize(cov_int), lower = -10*sig1+mu1, upper = 25*sig1+mu1)$value
}
mom_joint_inner=Vectorize(mom_joint_inner)

mom_joint <- function(mu1,r, d1, d2) {
  integrate(function(x2) {
    mom_joint_inner(x2, mu1,r, d1,d2) 
  }, lower = -10*sig2+mu2, upper = 25*sig2+mu2)$value
}
mom_joint = Vectorize(mom_joint)

# computed for each r ...
mom_var <- function(r){
  Cv12 <- mom_joint(mu1, r, 1, 1) - (mu1 + sig1 * gamma) * (mu2 + sig2 * gamma)
  Cw12 <- mom_joint(mu1, r, 2, 2) - mom_marg(mu1, sig1, 2) * mom_marg(mu2, sig2, 2)
  Cvw12 <- mom_joint(mu1, r, 1, 2) - mom_marg(mu2, sig2, 2) * (mu1 + sig1 * gamma)
  Cvw21 <- mom_joint(mu1, r, 2, 1) - mom_marg(mu1, sig1, 2) * (mu2 + sig2 * gamma)
  
  mom_g <- function(g1, g2) {
    alpha1 <- (W2 * (g1 * Cv12 + g2 * Cvw21) - Cvw2 * (g2 * Cw12 + g1 * Cvw12)) /
      (g1 * (V2 * W2 - Cvw2^2))
    alpha2 <- (V2 * (g2 * Cw12 + g1 * Cvw12) - Cvw2 * (g1 * Cv12 + g2 * Cvw21)) /
      (g2 * (V2 * W2 - Cvw2^2))
    
    Vmu1 <- V1 - 2 * alpha1 * Cv12 + alpha1^2 * V2
    Vmu2 <- W1 - 2 * alpha2 * Cw12 + alpha2^2 * W2
    Cmu12 <- Cvw1 - (alpha1 * Cvw21 + alpha2 * Cvw12 - alpha1 * alpha2 * Cvw2)
  
    return(g1^2 * Vmu1 + 2 * g1 * g2 * Cmu12 + g2^2 * Vmu2)
  }
  V11 <- mom_g(g11, g12)
  V22 <- mom_g(g21, g22)
  return(c(V11,V22))
}

# -------------------
# Generate and save results
# -------------------
Var_all <- data.frame(dep = numeric(), Vmm = numeric(), Vms = numeric(), Vss = numeric())
# JML
for (rr in r_seq) {
  temp <- as.vector(Var_joint(rr, sig1))
  Var_all <- rbind(Var_all, data.frame(dep = rr, Vmm = temp[1], Vms = temp[2], Vss = temp[4]))
}
# MML
Var_all <- Var_all %>%
  mutate(
    Vmle1 = sig1^2 * varP(a1, b1) * (1 - (integrate_cov_mle(dep, a1, b1) / varP(a1, b1))^2),
    Vmle2 = sig1^2 * varP(a2, b2) * (1 - (integrate_cov_mle(dep, a2, b2) / varP(a2, b2))^2)
  ) 
# MoM
Var_all = Var_all %>%
  rowwise() %>%
  mutate(var = list(mom_var(dep)),
         Vmom1 = var[1],
         Vmom2 = var[2]) %>%
  dplyr::select(-var)

# -------------------
# Plot
# -------------------

p <- Var_all %>%
  pivot_longer(c("Vmm", "Vss", "Vmle1", "Vmle2", "Vmom1", "Vmom2"),
               names_to = "Vtype", values_to = "variance") %>%
  mutate(
    Ptype = case_when(
      Vtype %in% c("Vmm", "Vmle1", "Vmom1") ~ "mu",
      Vtype %in% c("Vss", "Vmle2", "Vmom2") ~ "sigma"
    ),
    Vtype = case_when(
      Vtype %in% c("Vmm", "Vss") ~ "JML",
      Vtype %in% c("Vmle1", "Vmle2") ~ "MML",
      Vtype %in% c("Vmom1", "Vmom2") ~ "MoM"
    )
  ) %>%
  ggplot(aes(x = dep, y = variance, color = Vtype)) +
  geom_point() + geom_line() +
  geom_hline(data = data.frame(
    Vtype = rep(c("JML", "MML", "MoM"), each = 2),
    Ptype = rep(c("mu", "sigma"), times = 3),
    variance = c(asymp.var[1, 1], asymp.var[2, 2],
                 sig1^2 * varP(a1, b1), sig1^2 * varP(a2, b2),
                 g11^2 * V1 + 2 * g11 * g12 * Cvw1 + g12^2 * W1,
                 g21^2 * V1 + 2 * g21 * g22 * Cvw1 + g22^2 * W1)
  ),
  aes(yintercept = variance, color = Vtype),
  linetype = "dashed") +
  xlab("r") + ylab("Asympt Variance") +
  facet_grid(~Ptype) +
  labs(color = "Method") +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 13)
  )

ggsave("../images/Gumbel_jointParam.png", plot = p, width = 8, height = 3.7)
print(p)
