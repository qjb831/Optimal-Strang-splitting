# ---------------- PARAMETERS ----------------
sigma_vec <- c(0.5, 1.5)
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

b_seq <- seq(-1.5,1.5, by = 0.01)
y_seq <- seq(-0.37,0.37, by = 0.01)

data <- NULL
for(sigma in sigma_vec){
    df <- cross_join(data.frame(b_seq),data.frame(y_seq))
    colnames(df) <- c("b_seq","y_seq")
    df$epsilon <- epsilon
    df$sigma <- sigma
    func <- Vectorize(function(b,y) avg_absbias(b,y,epsilon,sigma))
    df$z <-  func(df$b_seq,df$y_seq) 
    data <- rbind(data,df)
}

root_curves <- lapply(y_seq, function(y){
  r <- get_roots(y)
  if(length(r) == 3){
    data.frame(
      y = y,
      left = r[1],
      middle = r[2],
      right = r[3],
      test = -r[1]/2
    )
  }
}) %>% bind_rows()

append_y <- function(string) 
  TeX(paste("$y = $", string)) 
append_sigma <- function(string) 
  TeX(paste("$\\sigma = $", string)) 


ggplot(data, aes(b_seq, y_seq, fill=z)) + 
  geom_raster() +
  scale_y_continuous(breaks=seq(-0.3,0.3,0.15)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  geom_line(data = root_curves,
            aes(x = left, y = y),
            color = "black",
            linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_line(data = root_curves,
            aes(x = right, y = y),
            color = "black",
            linewidth = 0.6,
            inherit.aes = FALSE) +
  geom_line(data = root_curves,
            aes(x =middle, y = y),
            color = "black",
            linewidth = 0.6,
            linetype=2,
            inherit.aes = FALSE) +
  theme_bw(base_size=16)+
  theme(panel.grid.major = element_blank(), # remove the ugly lines
        panel.grid.minor = element_blank())+
  theme(legend.position="right",legend.key.size= unit(1.2, "cm"),)+
  theme(strip.text.x = element_text(size = 16, color = "black", face = "bold"),
        strip.text.y = element_text(size = 16, color = "black", face = "bold"),
        strip.text = element_text(face="bold"),
        strip.background = element_rect(fill="grey")) +
  facet_wrap(~sigma ,    labeller = as_labeller(append_sigma, 
                                               default = label_parsed)            )+
  labs(x="c",y="y",fill=TeX('Mean $|B_\\theta(x,c)|$'))

