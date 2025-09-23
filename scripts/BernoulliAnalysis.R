# Bernoulli Analysis

library(tidyverse)
library(mvtnorm)

#################################
#### Gumbel–Hougaard Copula  ####
#################################

# Copula function C(p1, theta)
p2 <- 0.5
th <- 5 

Cop <- function(p1, th) {
  exp(-((-log(p1))^th + (-log(p2))^th)^(1/th))
}

# Asymptotic Variance for Baseline Estimator
AVar_Base <- function(p1, th) {
  p1 * (1 - p1)
}

# Asymptotic Variance for Multi-Fidelity Estimator
AVar_MF <- function(p1, th) {
  p1 * (1 - p1) * (1 - (Cop(p1, th) - p1 * p2)^2 / (p1 * (1 - p1) * p2 * (1 - p2)))
}
# Plot for Multiple Dependence Parameters (Gumbel–Hougaard)
p11 <- seq(0.001, 0.999, 0.001)
th <- c(0.1, 0.5, 1)  # r values
label_r <- function(value) {
  paste0("r = ", value)
}

p <- data.frame(p1 = rep(p11, times = length(th)), dep = rep(th, each = length(p11))) %>%
  mutate(Baseline = AVar_Base(p1, 1 / dep), MF = AVar_MF(p1, 1 / dep)) %>%
  pivot_longer(c("Baseline", "MF"), values_to = "variance", names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("MF", "Baseline"))) %>%
  ggplot(aes(p1, variance, linetype = Method)) +
  geom_line() +
  facet_grid(~dep, labeller = labeller(dep = label_r)) +
  ylab("Asympt Variance") +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 13)
  )

ggsave(filename = "../images/Bernoulli-G.png", plot = p, width = 7, height = 3.2)
print(p)


##############################
###### Gaussian Copula  ######
##############################

p2 <- 0.5
rh <- 0.95

Cop <- function(p1, rho, p2val = p2) {
  pmvnorm(
    lower = c(-Inf, -Inf),
    upper = c(qnorm(p1), qnorm(p2val)),
    sigma = matrix(c(1, rho, rho, 1), 2, 2),
    abseps = 1e-12,
    maxpts = c(250, 250, 250, 250)
  )[1]
}
Cop <- Vectorize(Cop, c("p1", "rho", "p2val"))

# Plot for Multiple Dependence Parameters (Gaussian Copula)
th <- c(0, 0.75, 0.95)  # rho values

label_rho <- function(value) {
  lapply(value, function(v) bquote(rho == .(v)))
}

p <- data.frame(p1 = rep(p11, times = length(th)), dep = rep(th, each = length(p11))) %>%
  mutate(Baseline = AVar_Base(p1, dep), MF = AVar_MF(p1, dep)) %>%
  pivot_longer(c("Baseline", "MF"), values_to = "variance", names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("MF", "Baseline"))) %>%
  ggplot(aes(p1, variance, linetype = Method)) +
  geom_line() +
  facet_grid(~dep, labeller = labeller(dep = label_rho, .default = label_parsed)) +
  ylab("Asympt Variance") +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 13)
  )

ggsave(filename = "../images/Bernoulli-N.png", plot = p, width = 7, height = 3.2)
print(p)


##############################
######   Special Case   ######
##############################

# Asymptotic Variance: Method of Moments
AVar_Base <- function(p1) {
  2 * p1 * (1 - p1 / 2)
}

# Asymptotic Variance: Joint MLE
AVar_MLE <- function(p1) {
  (1 / (2 * p1) + 1 / (6 * (3 - p1)) + 2 / (3 * (3 - 2 * p1)))^(-1)
}

p11 <- seq(0.001, 0.999, 0.001)
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols <- gg_color_hue(3)

p <- data.frame(p = p11) %>%
  mutate(MoM = AVar_Base(p), JML = AVar_MLE(p)) %>%
  pivot_longer(c("MoM", "JML"), values_to = "variance", names_to = "Method") %>%
  mutate(Method = factor(Method, levels = c("JML", "MoM"))) %>%
  ggplot(aes(p, variance, color = Method, linetype = Method)) +
  geom_line(linewidth = 1) +
  ylab("Asympt Variance") +
  scale_color_manual(values = c("JML" = cols[1], "MoM" = cols[3])) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 13)
  )

ggsave(file = "../images/Ber-special.png", plot = p, width = 6, height = 4)
print(p)
