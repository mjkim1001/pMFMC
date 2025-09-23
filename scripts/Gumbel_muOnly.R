# Bivariate Gumbel Analysis
# Estimating mu only

library(MASS)
library(tidyverse)
library(ggplot2)
library(evd)

# Parameters
sig1 <- 4
sig2 <- 1
mu1 <- 2
mu2 <- 2
r <- 0.75
r_seq <- seq(0.1, 1, by = 0.05)
Var_dep <- data.frame(dep = r_seq)

# Simulated data plot
set.seed(1)
simdata.joint <- rbvevd(100, dep = r, model = "log",
                        mar1 = c(mu1, sig1, 0),
                        mar2 = c(mu2, sig2, 0))
colnames(simdata.joint) <- c("hifi", "lofi")

ggplot(as.data.frame(simdata.joint)) +
  geom_point(aes(x = lofi, y = hifi))

#### JML ####
fA <- function(u, v, r) exp(-u / r) + exp(-v / r)
dAdm <- function(u, r) exp(-u / r) / (r * sig1)

I_int <- function(u, v, r) {
  A <- fA(u, v, r)
  B <- dAdm(u, r)
  term1 <- - (r * (r - 1) * A^(r - 2) * B^2 + A^(r - 1) * B / sig1) *
    (1 - 1 / (A^r - (r - 1) / r))
  term2 <- (r - 2) * (B / (r * sig1 * A) - B^2 / A^2)
  term3 <- - (r * A^(r - 1) * B)^2 / (1 / r - 1 + A^r)^2
  return(term1 + term2 + term3)
}

integrate_inner <- function(v, r) {
  int_u <- function(u) {
    I_int(u, v, r) * dbvevd(c(u, v), dep = r, model = "log",
                            mar1 = c(0, 1, 0), mar2 = c(0, 1, 0))
  }
  integrate(Vectorize(int_u), -12, 25)$value
}
integrate_inner <- Vectorize(integrate_inner)

integrate_outer <- function(r) {
  integrate(function(v) integrate_inner(v, r), -12, 25)$value
}
Var_est <- Vectorize(function(r) -1 / integrate_outer(r))

Var_dep <- Var_dep %>% mutate(Vjoint = Var_est(dep))

#### MoM ####
# To compute Cov(Y1, Y2)
integrate_cov_inner <- function(y2, mu1, r) {
  cov_int <- function(y1) {
    y1 * y2 * dbvevd(c(y1, y2), dep = r, model = "log",
                     mar1 = c(mu1, sig1, 0), mar2 = c(mu2, sig2, 0))
  }
  integrate(Vectorize(cov_int), -12 * sig1 + mu1, 25 * sig1 + mu1)$value
}
integrate_cov_inner <- Vectorize(integrate_cov_inner)

integrate_cov <- function(mu1, r) {
  integrate(function(y2) integrate_cov_inner(y2, mu1, r),
            -12 * sig2 + mu2, mu2 + 25 * sig2)$value
}
integrate_cov <- Vectorize(integrate_cov)

gamma <- -digamma(1)
Var_dep <- Var_dep %>%
  mutate(correlation = 6 * (integrate_cov(mu1, dep) -
                              (mu1 + sig1 * gamma) * (mu2 + sig2 * gamma)) / (pi^2 * sig1 * sig2),
         Vmme = pi^2 * sig1^2 / 6 * (1 - correlation^2))

ggplot(Var_dep, aes(dep, correlation)) +
  geom_point() + geom_line()

#### MML ####
integrate_cov_mle_inner <- function(v, r) {
  cov_int <- function(u) {
    (1 - exp(-u)) * (1 - exp(-v)) *
      dbvevd(c(u, v), dep = r, model = "log",
             mar1 = c(0, 1, 0), mar2 = c(0, 1, 0))
  }
  integrate(Vectorize(cov_int), -12, 25)$value
}
integrate_cov_mle_inner <- Vectorize(integrate_cov_mle_inner)

integrate_cov_mle <- Vectorize(function(r)
  integrate(function(v) integrate_cov_mle_inner(v, r), -12, 25)$value)

Var_dep <- Var_dep %>%
  mutate(Vmle = sig1^2 * (1 - integrate_cov_mle(dep)^2))

# --- Final plot for mu1 only ---
p <- Var_dep %>%
  pivot_longer(c("Vjoint", "Vmme", "Vmle"),
               names_to = "Vtype", values_to = "variance") %>%
  mutate(Vtype = recode(Vtype, Vjoint = "JML", Vmme = "MoM", Vmle = "MML"),
         Vtype = factor(Vtype, levels = c("JML", "MoM", "MML")),
         Ptype = "mu") %>%
  ggplot(aes(x = dep, y = variance, color = Vtype)) +
  geom_point() + geom_line() +
  geom_hline(data = data.frame(
    Vtype = c("MoM", "MML"),
    variance = c(pi^2 * sig1^2 / 6, sig1^2)),
    aes(yintercept = variance, color = Vtype),
    linetype = "dashed") +
  xlab("r") + ylab("Asympt Variance") +
  labs(color = "Method") +
  facet_grid(~Ptype) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 13)
  )

ggsave(filename = "../images/Gumbel_muOnly.png", plot = p, width = 5.5, height = 4)
print(p)
