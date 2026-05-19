library(tidyverse)
library(ggplot2)
library(latex2exp)
library(cowplot)

y_vec <- c( 0, 0.1)
sigma_vec <- c( 0.75, 1.5)
epsilon <- 0.1


B <- function(b,x,y,sigma){
  denom <- 1/(24*epsilon^3)
  top <- (-6*b^3+3*y)*x^4 +
    x^3*(-36*b^4-24*sigma^2*epsilon+36*b^2-8) +
    x^2*(54*(b^3-y)*(b^2-1/3)) +
    x*(-24*b^6+54*b^2*epsilon*sigma^2+36*b^3*y-12*sigma^2*epsilon-12*y^2) +
    6*b^5-9*b^4*y+b^3*(-18*sigma^2*epsilon-4) +
    6*b^2*y+12*epsilon*sigma^2*y
  return(top*denom)
}
absB <- function(b,x,y,sigma) abs(B(b,x,y,sigma))


x_seq = data.frame(seq(-1.4,1.4,0.01))
b_seq = data.frame(seq(-1.4,1.4,0.01))


data <- NULL
for(i in y_vec){
  for(sigma in sigma_vec){
    df <- cross_join(x_seq,b_seq)
    colnames(df) <- c("x_seq","b_seq")
    df$y <- i
    df$sigma <- sigma
    df$z <- pmax( absB(df$b_seq, df$x_seq, df$y, df$sigma), 1e-4)
    roots <- sort(Re(polyroot(c(-i,1,0,-1))))
    df$LR <- roots[1]
    df$MR <- roots[2]
    df$RR <- roots[3]
    data <- rbind(data,df)
  }
}

rootfunction <- function(y){ 
  roots <- Re(polyroot(c(-y,1,0,-1))) 
  return(c(min(roots),median(roots),max(roots)))}


vlines_df1 <- data %>%
  distinct(y, sigma, rootfunction(y)[1])
vlines_df2 <- data %>%
  distinct(y, sigma, rootfunction(y)[2])
vlines_df3 <- data %>%
  distinct(y, sigma, rootfunction(y)[3])

head(vlines_df1)

vlines_df2 <- heat_df %>%
  distinct(y, sigma, x_star2, y_lab, sigma_lab) %>%
  mutate(
    y_lab = factor(y_lab, levels = y_labels_ordered),
    sigma_lab = factor(sigma_lab, levels = sigma_labels_ordered)
  )
vlines_df3 <- heat_df %>%
  distinct(y, sigma, unstable_fix, y_lab, sigma_lab) %>%
  mutate(
    y_lab = factor(y_lab, levels = y_labels_ordered),
    sigma_lab = factor(sigma_lab, levels = sigma_labels_ordered)
  )





summary(data)

ggplot(data, aes(b_seq, x_seq, fill=z)) + 
  geom_raster() +
  scale_y_continuous(breaks=seq(-1,1,0.5)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()+
  facet_grid(sigma~y,labeller=label_both)+
  geom_vline(xintercept=data$MR,linetype=2,inherit.aes=FALSE)+
  geom_vline(xintercept=data$LR,linetype=1,inherit.aes=FALSE)+
  geom_vline(xintercept=data$RR,linetype=1,inherit.aes=FALSE)+
  labs(x="c",y="x",fill="|B(x,c)|")





sort(Re(polyroot(c(-0.1,1,0,-1))))

df <- cross_join(x_seq,b_seq)
colnames(df) <- c("x_seq","b_seq")
data11 <- df
data12 <- df
data21 <- df
data22 <- df
data11$y <- 0
data11$sigma <- 0.75
data11$z <- pmax( absB(data11$b_seq, data11$x_seq, data11$y, data11$sigma), 1e-4)
data12$y <- 0
data12$sigma <- 1.5
data12$z <- pmax( absB(data12$b_seq, data12$x_seq, data12$y, data12$sigma), 1e-4)
data21$y <- 0.1
data21$sigma <- 0.75
data21$z <- pmax( absB(data21$b_seq, data21$x_seq, data21$y, data21$sigma), 1e-4)
data22$y <- 0.1
data22$sigma <- 1.5
data22$z <- pmax( absB(data22$b_seq, data22$x_seq, data22$y, data22$sigma), 1e-4)
p11 <- ggplot(data11, aes(b_seq, x_seq, fill=z)) + 
  geom_raster() +
  geom_vline(xintercept=0,linetype=2,col="darkslategrey")+
  geom_vline(xintercept=-1,linetype=1,col="darkslategrey")+
  geom_vline(xintercept=1,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=0,linetype=2,col="darkslategrey")+
  geom_hline(yintercept=-1,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=1,linetype=1,col="darkslategrey")+
  scale_y_continuous(breaks=seq(-1,1,0.5)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()+
  labs(x="c",y="x",fill="|B(x,c)|",title=TeX('$y=0,$ $\\sigma=0.75$'))
p12 <- ggplot(data12, aes(b_seq, x_seq, fill=z)) + 
  geom_raster() +
  geom_vline(xintercept=0,linetype=2,col="darkslategrey")+
  geom_vline(xintercept=-1,linetype=1,col="darkslategrey")+
  geom_vline(xintercept=1,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=0,linetype=2,col="darkslategrey")+
  geom_hline(yintercept=-1,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=1,linetype=1,col="darkslategrey")+
  scale_y_continuous(breaks=seq(-1,1,0.5)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()+
  labs(x="c",y="x",fill="|B(x,c)|",title=TeX('$y=0,$ $\\sigma=1.5$'))
p21 <- ggplot(data21, aes(b_seq, x_seq, fill=z)) + 
  geom_raster() +
  geom_vline(xintercept=0.1010313  ,linetype=2,col="darkslategrey")+
  geom_vline(xintercept=-1.0466805  ,linetype=1,col="darkslategrey")+
  geom_vline(xintercept=0.9456493,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=0.1010313  ,linetype=2,col="darkslategrey")+
  geom_hline(yintercept=-1.0466805  ,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=0.9456493,linetype=1,col="darkslategrey")+
  scale_y_continuous(breaks=seq(-1,1,0.5)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()+
  labs(x="c",y="x",fill="|B(x,c)|",title=TeX('$y=0.1,$ $\\sigma=0.75$'))
p22 <- ggplot(data22, aes(b_seq, x_seq, fill=z)) + 
  geom_raster() +
  geom_vline(xintercept=0.1010313  ,linetype=2,col="darkslategrey")+
  geom_vline(xintercept=-1.0466805  ,linetype=1,col="darkslategrey")+
  geom_vline(xintercept=0.9456493,linetype=2,col="darkslategrey")+
  geom_hline(yintercept=0.1010313  ,linetype=2,col="darkslategrey")+
  geom_hline(yintercept=-1.0466805  ,linetype=1,col="darkslategrey")+
  geom_hline(yintercept=0.9456493,linetype=2,col="darkslategrey")+
  scale_y_continuous(breaks=seq(-1,1,0.5)) +
  scale_x_continuous(breaks=seq(-1,1,0.5)) +
  scale_fill_gradient(low="yellow", high="red",trans="log10") +
  theme_classic()+
  labs(x="c",y="x",fill="|B(x,c)|",title=TeX('$y=0.1,$ $\\sigma=1.5$'))



plot_grid(p11, p12, p21, p22, 
          labels = NULL,
          ncol = 2, nrow = 2)
