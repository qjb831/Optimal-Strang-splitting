

y <- 0
sigma <- 0.5
epsilon <- 0.1


B <- function(b,x){
  c0 <- -sigma^2*x^3/epsilon^2 - sigma^2*x/(2*epsilon^2) - x^3/(3*epsilon^3)
  c1 <-  (3*x^2)/(4*epsilon^3) + x^4/(8*epsilon^3) + sigma^2/(2*epsilon^2)
  c2 <- (3*x^3)/(2*epsilon^3) + (9*sigma^2*x)/(4*epsilon^2) - x/(2*epsilon^3)
  c3 <- -(15*x^2)/(4*epsilon^3) + 1/(12*epsilon^3) - (5*sigma^2)/(4*epsilon^2) - (3*x^4)/(8*epsilon^3)
  c4 <- -(3*x^3)/(2*epsilon^3) + (5*x)/(2*epsilon^3)
  c5 <- (9*x^2)/(2*epsilon^3) - 3/(8*epsilon^3)
  c6 <- -3*x/epsilon^3
  c7 <- 3/(8*epsilon^3)
  return(c0+c1*b+c2*b^2+c3*b^3+c4*b^4+c5*b^5+c6*b^6+c7*b^7)
}



inv <- function(x){exp(2*(0.5*x^2-0.25*x^4-y*x)/(sigma^2*epsilon))}
Z <- integrate(inv,-Inf,Inf)$value
p <- function(x){ inv(x)/Z}

avg_bias <- function(b){
  integrand <- function(x){B(b,x)*p(x)}
  integrate(integrand,-Inf,Inf)$value
}
absbias <- function(b,x){abs(B(b,x))}
avg_absbias <- function(b){
  integrand <- function(x){absbias(b,x)*p(x)}
  integrate(integrand,-Inf,Inf)$value
}

#avg_bias(1)

bseq <- seq(-1,1,0.01)
Bias <- rep(NA,length(bseq))
for(i in 1:length(bseq)){
  Bias[i] <- avg_bias(bseq[i])}
plot(bseq,Bias)