#################################################################################### 
# Script for Natural disasters paper
####################################################################################
# !diagnostics off

##### 
# Clean working space, packages and functions
##### 

rm(list=ls(all=TRUE))
library(Rsolnp)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(ggrepel)

##### 
# Load required functions to perform the analysis
##### 

source('corefx.R')

#####
# Specify the parameters of the model
#####

# Biological parameters
r <- .3 # Growth rate
cc <- 100 # Carrying capacity

# Economic parameters (TT,x0,p,w,q,c,r,cc,delta,beta,i_a,i_b)
p <- 5 # Price
w <- 1 # Cost
q <- .005 # Catchability coefficient
c <- 10 # Cost of capital 
delta <- 0.1 # Depreciation
beta <- 0.95 # Discount factor

##### 
# Simulations
##### 

# Initial conditions
TT <- 80
k0 <- 0
nk0 <- 0
X0 <- 100
d <- 0 # No subsidies

# Number of firms to evaluate the model 
NN = c(seq(1,50,1))

# Maximum number of iterations to calculate the closed loop
maxit = 100

# Empty data frames to store steady state conditions
A <- data.frame(N = NN, SR = NN) 
C <- data.frame(N = NN, SR = NN) 

# Run simulations for different number of firms 
for(mm in 1:length(NN)){
  N = NN[mm]
  
  i <- matrix(0,nrow=TT,ncol=N)
  
  print(c("Currently evaluating equilibrium for ",N,"firms"))
  for (nn in 1:maxit){ # Max number of iterations for equilibrium
    for(jj in 1:ncol(i)){ # Index over firm jj
      # Initiate calculations, assume other players have zero investment 
      if(N==1){
        ni <- rep(0,TT) # Sole owner has no competiion
      }else if(N==2){
        ni <- i[,-jj] # N=2 It's just the opposite column
      }else{
        ni <- rowSums(i[,-jj]) # Add all other columns 
      }
    
      # Estimate the best response of firm jj
      OR <- solnp(rep(0.1/N,TT), # Starting values
                  fun = nash.handle, # Function to minimize
                  LB = rep(0,TT), # Lower bound for decision variables
                  UB = rep(1e3,TT), # Upper bound for decision variables
                  control = list(trace = 0)) # Omit output
      
      # Store the results
      i[,jj] <- OR$par 
    }
    # Visulaize convergence
    matplot(i,type="l",xlab="Time",ylab="Optimal investment",
            main=c("Iteration ",nn))
    
    # Stopping criteria
    if(N==1){break}else{
      if(sum(i[,1]-i[,N])^2<1e-6){break}
    }
    if(nn==maxit){print("EQUILIIBRUM CONDITIONS ARE NOT SATISFIED!")} 
  }
  # Store results
  A$SR[mm] <- list(fishery.simulation.nash(TT,N,X0,k0,nk0,p,w,q,c,r,cc,delta,beta,rowSums(cbind(i[,-1],rep(0,TT))),i[,1],d))
  
  # Cooperative counterfactual
  OR <- solnp(rep(0.1/N,TT), # Starting values
              fun = coop.handle, # Function to minimize
              LB = rep(0,TT), # Lower bound for decision variables
              UB = rep(1e3,TT), # Upper bound for decision variables
              control = list(trace = 0)) # Omit output
  i <- OR$par
  # Store results
  C$SR[mm] <- list(fishery.simulation.coop(TT,N,X0,k0,p,w,q,c,r,cc,delta,beta,i))
}

# Store full set of competitive results as data frame
A <- A %>% rowwise() %>% 
  mutate(N.X = SR$Total.stock[TT/2],
         N.i = SR$Investment[TT/2],
         N.I = N*N.i,
         N.k = SR$Capital[TT/2],
         N.K = SR$Total.capital[TT/2],
         N.Pi = SR$Profits[TT/2],
         N.TPI = N*N.Pi,
         N.NPV = SR$NPV)
A$SR <- NULL

# Store full set of cooperative results as data frame
C <- C %>% rowwise() %>% 
  mutate(C.X = SR$Total.stock[TT/2],
         C.i = SR$Investment[TT/2],
         C.I = C.i,
         C.k = SR$Capital[TT/2],
         C.K = SR$Total.capital[TT/2],
         C.Pi = SR$Profits[TT/2],
         C.TPI = N*C.Pi,
         C.NPV = SR$NPV)
C$SR <- NULL

# Create table with results
R <- merge(A,C,by='N')

# Calculate path back to steady state after a natural disaster destroys 90% of the capital stock

di <- .9 # Percentage of capital destroyed

# Empty data frames to store results
AD <- data.frame(N = NN, SR = NN) 
CD <- data.frame(N = NN, SR = NN) 

# Evaluate recovery path for each N scenario
for(mm in 1:length(NN)){
  N = NN[mm]
  k0 <- R$N.k[mm] * di # Initial capital is reduced 
  nk0 <- k0 * (N-1)
  X0 <- R$N.X[mm]
  
  i <- matrix(0,nrow=TT,ncol=N)
  
  print(c("Currently evaluating equilibrium for ",N,"firms"))
  for (nn in 1:maxit){ # Max number of iteration iterations for equilibrium
    for(jj in 1:ncol(i)){ # Index over firm jj
      # Initiate calculations, assume other players have zero investment 
      if(N==1){
        ni <- rep(0,TT) # Sole owner has no competiion
      }else if(N==2){
        ni <- i[,-jj] # N=2 It's just the opposite column
      }else{
        ni <- rowSums(i[,-jj]) # Add all other columns 
      }
      
      # Estimate the best response of firm jj
      OR <- solnp(rep(0.1/N,TT), # Starting values
                  fun = nash.handle, # Function to minimize
                  LB = rep(0,TT), # Lower bound for decision variables
                  UB = rep(1e3,TT), # Upper bound for decision variables
                  control = list(trace = 0)) # Omit output
      
      # Store the results
      i[,jj] <- OR$par 
    }
    # Visulaize convergence
    matplot(i,type="l",xlab="Time",ylab="Optimal investment",
            main=c("Iteration ",nn))
    
    # Stopping criteria
    if(N==1){break}else{
      if(sum(i[,1]-i[,N])^2<1e-6){break}
    }
    if(nn==maxit){print("EQUILIIBRUM CONDITIONS ARE NOT SATISFIED!")} 
  }
  
  # Store results
  AD$SR[mm] <- list(fishery.simulation.nash(TT,N,X0,k0,nk0,p,w,q,c,r,cc,delta,beta,rowSums(cbind(i[,-1],rep(0,TT))),i[,1],d))
  
  # Cooperative counterfactual
  
  k0 <- R$C.k[mm] * .1
  X0 <- R$C.X[mm]
  
  OR <- solnp(rep(0.1/N,TT), # Starting values
              fun = coop.handle, # Function to minimize
              LB = rep(0,TT), # Lower bound for decision variables
              UB = rep(1e3,TT), # Upper bound for decision variables
              control = list(trace = 0)) # Omit output
  
  # Visualize convergence
  i <- OR$par
  CD$SR[mm] <- list(fishery.simulation.coop(TT,N,X0,k0,p,w,q,c,r,cc,delta,beta,i))
}

# Store competitive results as data frame
AD <- AD %>% rowwise() %>% 
  mutate(N.X = SR$Total.stock[TT/2],
         N.i = SR$Investment[TT/2],
         N.I = N*N.i,
         N.k = SR$Capital[TT/2],
         N.K = SR$Total.capital[TT/2],
         N.Pi = SR$Profits[TT/2],
         N.TPI = N*N.Pi,
         N.NPV = SR$NPV)

# Store cooperative results as data frame
CD <- CD %>% rowwise() %>% 
  mutate(C.X = SR$Total.stock[TT/2],
         C.i = SR$Investment[TT/2],
         C.I = C.i,
         C.k = SR$Capital[TT/2],
         C.K = SR$Total.capital[TT/2],
         C.Pi = SR$Profits[TT/2],
         C.TPI = N*C.Pi,
         C.NPV = SR$NPV)

RD <- merge(AD,CD,by='N')

AD$SR <- NULL
CD$SR <- NULL

# Create table for comparison
RR <- merge(AD,CD,by='N')

##### 
# Plot simulation results
##### 

# Grab steady state values
ND <- data.frame(time = rep((1:(TT*2)), times = length(NN)),
                 N = rep(NN, times = 1 , each = (TT*2)),
                 id = rep((1:length(NN)), times = 1 , each = (TT*2)))

ND <- ND %>% rowwise() %>%
  mutate(X = R$N.X[id],
         K = R$N.K[id],
         I = R$N.I[id],
         Pi = R$N.Pi[id])

for (nn in 1:length(NN)){
  ND[ND$time>TT & ND$id==nn,4] = RD$SR.x[[nn]][['Total.stock']]
  ND[ND$time>TT & ND$id==nn,5] = RD$SR.x[[nn]][['Total.capital']]
  ND[ND$time>TT & ND$id==nn,6] = RD$SR.x[[nn]][['Total.investment']]
  ND[ND$time>TT & ND$id==nn,7] = RD$SR.x[[nn]][['Profits']]
}

ND <- ND %>% group_by(N) %>% 
  mutate(i = I/N,
         si = i/i[1],
         k = K/N,
         rr = i/k,
         tt = ifelse(time<=TT,0,time-TT),
         counter = ifelse(time<=TT,0,(c/2*i[1]^2)*beta^(tt-1)),
         actual = ifelse(time<=TT,0,(c/2*i^2)*beta^(tt-1)))

# PLot a
R %>% ggplot() +
  geom_point(aes(N.X,N.K),size = 7) +
  geom_label_repel(aes(N.X,N.K,label = N),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   direction = "y",
                   hjust=0.5,
                   size=10) +
  xlab(TeX('Steady State Resource ($X^*$)')) +
  ylab(TeX('Steady State Capital ($K^*$)')) +
  ggtitle(TeX('Resource - Capital (a)')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title=element_text(size=30,face="bold"),
        axis.text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
  ggsave(file="f1.eps",width=12, height=9)

# Plot b
ND %>% filter(time>(TT-10) & time<TT+30) %>%
  ggplot(aes(group=N)) +
  geom_line(aes(time,si,linetype=as.factor(N)),size=1.5) +
  xlab(TeX('Time')) +
  ylab(TeX('Relative Investment ($i/i^*$)')) +
  ggtitle(TeX('Response in Investment (b)')) +
  labs(linetype=TeX('Firms ($N$)'))+
  annotate("text", label="ND", x=TT-5, y=1.12, size=10) +
  annotate("segment", x=TT-5, y=1.1, xend=TT-2, yend=1.04, size=2, 
           arrow=arrow(length=unit(.5, "cm")))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title = element_text(size=30,face="bold"),
        axis.text = element_text(size=20),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = c(.9, .5),
        legend.background = element_rect(color = "black", size = 1, linetype = "solid"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        legend.key.height=unit(2,"line"))+ 
  ggsave(file="f2.eps",width=12, height=9)

# Plot c
hh <- ND %>% filter(time==TT+1) %>% select(N,rr)
spline_int <- as.data.frame(spline(hh$N,hh$rr))
ND %>% filter(time==TT+1) %>%
  ggplot() +
  geom_line(data = spline_int, aes(x = x, y = y), size = 3)+
  geom_point(aes(N,rr), size = 7) +
  xlab(TeX('Firms ($N$)')) +
  ylab(TeX('Rate of investment ($i/k$)')) +
  ggtitle(TeX('Investment Pulse (c)')) +
  labs(linetype=TeX('$N$'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title=element_text(size=30,face="bold"),
        axis.text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
  ggsave(file="f3.eps",width=12, height=9)

# Plot d
hh <- ND %>% filter(time<TT*1.5) %>% 
  group_by(N) %>% 
  summarize(actual = sum(actual),counter = sum(counter)) %>%
  mutate(loss =(actual - counter)/counter)

spline_int <- as.data.frame(spline(hh$N,hh$loss))
hh %>%  ggplot() +
  geom_line(data = spline_int, aes(x = x, y = y), size = 3)+
  geom_point(aes(N,loss), size = 7) +
  xlab(TeX('Firms ($N$)')) +
  ylab(TeX('Relative cost of rebuilding')) +
  ggtitle(TeX('Direct Capital Impact (d)')) +
  labs(linetype=TeX('$N$'))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=30),
        axis.title=element_text(size=30,face="bold"),
        axis.text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
  ggsave(file="f4.eps",width=12, height=9)
