####################################################################################### 
# Functions for the numerical simulations of capital investment
#######################################################################################

# General functions
#####
# Stock growth
stock.growth<-function(x,r,cc){
  g=r*x*(1-x/cc)+x
  g
}

# Profit function
profit.function<-function(y,i,x,p,w,q,c,N){
  Pi = (p-w/(q*x))*y - c/2*(i)^2
  Pi
}

# Competitive extraction simulation
fishery.simulation.nash <- function(TT,N,X0,k0,nk0,p,w,q,c,r,cc,delta,beta,ni,i,d){
  
  x <- rep(0,TT) # Biomass vector
  
  y <- rep(0,TT) # Extraction
  ny <- rep(0,TT) # Others Extraction
  
  Y <- rep(0,TT) # Total Extraction
  
  k <- rep(0,TT) # Firm's capital
  nk <- rep(0,TT) # Others' capital
  
  K <- rep(0,TT) # Total capital
  
  I <- rep(0,TT) # Total investment
  
  Pi <- rep(0,TT) # Individual profits
  
  DPi <- rep(0,TT) # Discounted profit vector for the firm
  
  for (t in 1:TT){
    if(t == 1){
      x[t] = X0 # Initial biomass
      k[t] = k0
      nk[t] = nk0
    }else{
      x[t] = stock.growth(x[t-1],r,cc)-Y[t-1] # Growth dynamics
      k[t] = (1-delta)*k[t-1]+i[t-1]
      nk[t] = (1-delta)*nk[t-1]+ni[t-1]
    }
    
    K[t] <- k[t] + nk[t] # Total capital
    I[t] <- i[t] + ni[t] # Total investment
    
    y[t] <- x[t]*q*k[t] # Individual harvest
    ny[t] <- x[t]*q*nk[t] # Others' harvest
    
    Y[t] <- y[t] + ny[t] # Total harvest
    
    if(d==1 & t<=5){ # Handle that alllows for 70% subsidy during first 5 periods
      kc = c*.3
    }else{
      kc = c
    }
    
    Pi[t] = profit.function(y[t],i[t],x[t],p,w,q,kc,N) # Individual profits
    
    DPi[t] = beta^(t-1)*Pi[t] # Discounted profits
  }
  
  NPV = sum(DPi) # Net present value for the firm
  
  list(Total.stock = x, Harvest = y, Total.harvest  = Y, Investment = i, Total.investment = I,
       Capital = k, Total.capital = K, Profits = Pi, NPV = NPV
  )
} 

# Competitive extraction handle for optimization
nash.handle <- function(i){
  OH <- fishery.simulation.nash(TT,N,X0,k0,nk0,p,w,q,c,r,cc,delta,beta,ni,i,d)
  z  <- -OH$NPV
  z
}

# Cooperative extraction simulation
fishery.simulation.coop <- function(TT,N,X0,k0,p,w,q,c,r,cc,delta,beta,i){
  
  x <- rep(0,TT) # Biomass vector
  
  y <- rep(0,TT) # Extraction
  
  Y <- rep(0,TT) # Total Extraction
  
  k <- rep(0,TT) # Firms' capital
  
  K <- rep(0,TT) # Total capital
  
  I <- rep(0,TT) # Total investment
  
  Pi <- rep(0,TT) # Individual profits
  
  DPi <- rep(0,TT) # Discounted profit vector for the firms
  
  for (t in 1:TT){
    if(t == 1){
      x[t] = X0 # Initial biomass
      k[t] = k0
    }else{
      x[t] = stock.growth(x[t-1],r,cc)-Y[t-1] # Growth dynamics
      k[t] = (1-delta)*k[t-1]+i[t-1]
    }
    
    K[t] <- N*k[t] # Total capital
    I[t] <- N*i[t] # Total investment
    
    y[t] <- N*x[t]*q*k[t] # Individual harvest
    
    Y[t] <- N*y[t] # Total harvest
    
    Pi[t] = profit.function(y[t],i[t],x[t],p,w,q,c,N) # Individual profits
    
    DPi[t] = N*beta^(t-1)*Pi[t] # Discounted TOTAL profits
  }
  
  NPV = sum(DPi) # Net present value for the FISHERY
  
  list(Total.stock = x, Harvest = y, Total.harvest  = Y, Investment = i, Total.investment = I,
       Capital = k, Total.capital = K, Profits = Pi, NPV = NPV
  )
} 

# Competitive extraction handle for optimization
coop.handle <- function(i){
  OH <- fishery.simulation.coop(TT,N,X0,k0,p,w,q,c,r,cc,delta,beta,i)
  z  <- -OH$NPV
  z
}

