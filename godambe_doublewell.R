
# DOUBLE WELL
compute_expression <- function(x_old1,h, par, center) {
  sigma <- par[3]
  
  Omega_wrap <- function(p) {
    Omega_scalar(h, p, center)
  }
  
  Omega <- Omega_wrap(par)
  OmegaInv <- 1/Omega
  
  # first derivative
  grad_Omega <- grad(Omega_wrap, par)
  dO1 <- grad_Omega[3]
  dOinv1 <- -dO1/Omega^2
  
  # second derivative
  hess_Omega <- hessian(Omega_wrap, par)
  d2O11 <- hess_Omega[3,3]
  d2Oinv11 <- 2*dO1^2/Omega^3-d2O11/Omega^2
  
  #EG2 <- 2*(1/sigma^2)*(2*sigma)*(1/sigma^2)*2*sigma
  
  #monte carlo
  N <- 600
  em_h <- 0.0001
  
  times <- seq(0, N * em_h, length.out = N + 1)
  F <- function(x){
    (x[1]-x[1]^3-par[2])/par[1]
  }
  G <- function(x) par[3]
  x_new1 <- EM_step(F, G, x_old1, times, n = 1000)
  fh_new1 <- apply(
    X = x_new1,
    MARGIN = 1,  
    FUN = fh_scalar,
    h = -h/2,
    par = par,
    b = center,
    const_term_in_N = FALSE
  )
  
  fh_old1 <- fh_scalar(x_old1,h/2,par,center,FALSE)
  mu1 <- mu_scalar(fh_old1,h,par,center,FALSE)
  z <- fh_new1 - mu1
  EG2_1 <- mean((dO1*OmegaInv+z^2*dOinv1)^2)
  S <- mean(dOinv1*dO1+OmegaInv*d2O11+d2Oinv11*z^2)
  return(S^2/EG2_1)
}
S_wrapper <- function(b) {
  sapply(b, function(bi) {
    compute_expression(
      x= -0.9,
      h = 0.1,
      par = c(0.2,-0.05,1.2),
      bi
    )
  })
}

curve(S_wrapper, -1.2, 1.2,n=40,xname = "b")
optimize(S_wrapper, interval = c(-5,5), maximum = TRUE)