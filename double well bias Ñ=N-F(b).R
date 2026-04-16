

y <- 0
sigma <- 0.7
epsilon <- 0.1


B <- function(b,x){
  
  denom <- 1/(24*epsilon^3)
  top <- (-6*b^3+3*y)*x^4 +
    x^3*(-36*b^4-24*sigma^2*epsilon+36*b^2-8) +
    x^2*(54*(b^3-y)*(b^2-1/3)) +
    x*(-24*b^6+54*b^2*epsilon*sigma^2+36*b^3*y-12*sigma^2*epsilon-12*y^2) +
    6*b^5-9*b^4*y+b^3*(-18*sigma^2*epsilon-4) +
    6*b^2*y+12*epsilon*sigma^2*y
  return(top*denom)
}
absB <- function(b,x) abs(B(b,x))

x_seq = seq(-1.4,1.4,0.001)
b_seq = seq(-1.4,1.4,0.001)
z = outer(b_seq, x_seq, absB)
colnames(z) = b_seq
rownames(z) = x_seq

dat = as.data.frame(z) %>% 
  rownames_to_column(var="x_seq") %>% 
  gather(b_seq, value, -x_seq) %>% 
  mutate(b_seq=as.numeric(b_seq), 
         x_seq=as.numeric(x_seq)) %>% 
  mutate(value = pmax(value, 1e-4))

ggplot(dat, aes(b_seq, x_seq, fill=value)) + 
  geom_raster() +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()



library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ---------------- PARAMETERS ----------------
sigma <- 0.7
epsilon <- 0.1

B <- function(b,x,y,epsilon,sigma){
  denom <- 1/(24*epsilon^3)
  top <- (-6*b^3+3*y)*x^4 +
    x^3*(-36*b^4-24*sigma^2*epsilon+36*b^2-8) +
    x^2*(54*(b^3-y)*(b^2-1/3)) +
    x*(-24*b^6+54*b^2*epsilon*sigma^2+36*b^3*y-12*sigma^2*epsilon-12*y^2) +
    6*b^5-9*b^4*y+b^3*(-18*sigma^2*epsilon-4) +
    6*b^2*y+12*epsilon*sigma^2*y
  return(top*denom)
}
get_roots <- function(y){
  r <- polyroot(c(-y, 1, 0, -1))  # x^3 - x + y = 0
  r <- Re(r[abs(Im(r)) < 1e-8])  # keep real roots
  sort(r)
}

inv <- function(x,y,sigma,epsilon){exp(2*(0.5*x^2-0.25*x^4-y*x)/(sigma^2*epsilon))}


avg_absbias <- function(b,y,epsilon,sigma){
  uns_root <- median(get_roots(y))
  Z <- integrate(inv,-Inf,Inf,y=y,sigma=sigma,epsilon=epsilon)$value
  p_inv <- 1/integrate(function(x,y,sigma,epsilon) {inv(x,y,sigma,epsilon) / Z} , -Inf,uns_root,y,sigma,epsilon)$value
  integrand <- function(x){
    abs(B(b,x,y,epsilon,sigma) * inv(x,y,sigma,epsilon) / Z)
  }
  
  
  abs(integrate(integrand, -Inf, uns_root)$value*p_inv)
}

b_seq <- seq(-1.2,-0.25, by = 0.001)
y_seq <- seq(-0.2, 0.2, by = 0.005)

z <- outer(b_seq, y_seq,
           Vectorize(function(b,y) avg_absbias(b,y,epsilon,sigma)))

rownames(z) <- b_seq
colnames(z) <- y_seq

# ---------------- DATAFRAME ----------------
dat <- as.data.frame(z) %>% 
  rownames_to_column("b_seq") %>% 
  pivot_longer(-b_seq, names_to="y_seq", values_to="value") %>% 
  mutate(
    b_seq = as.numeric(b_seq),
    y_seq = as.numeric(y_seq),
    value = pmax(value, 1e-6)
  )

# ---------------- ROOT CURVES ----------------
root_curves <- lapply(y_seq, function(y){
  r <- get_roots(y)
  if(length(r) == 3){
    data.frame(
      y = y,
      left = r[1],
      middle = r[2],
      right = r[3]
    )
  }
}) %>% bind_rows()

# ---------------- PLOT ----------------
ggplot(dat, aes(b_seq, y_seq, fill=value)) + 
  geom_raster() +
  
  # left root
  geom_line(data = root_curves,
            aes(x = left, y = y),
            color = "black",
            linewidth = 0.3,
            inherit.aes = FALSE) +
  
  # right root
  geom_line(data = root_curves,
            aes(x = right, y = y),
            color = "black",
            linewidth = 0.3,
            inherit.aes = FALSE) +
  
  scale_fill_gradient(low="yellow", high="red", trans="log10") +
  theme_classic() +
  labs(x="b", y="y", fill="avg |bias|")


