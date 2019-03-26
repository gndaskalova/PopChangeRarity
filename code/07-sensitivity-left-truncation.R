# Calculating population change using state-space models
# Based on Humbert et al. 2009 and Dennis et al. 2006

# Libraries ----
library(dplyr)
library(readr)
library(tidyr)

# Optional clean of the working environment
rm(list = ls())

# Load data ----
vertebrates <- read.csv("data/input/LPIdata_Feb2016.csv", na.strings = c("NA", "", " "))
colnames(vertebrates)[26:70] <- parse_number(colnames(vertebrates)[26:70])
# 16 624 rows in the wide version


# Remove populations which have completely identical abundance throughout the whole timeseries
vertebrates <- vertebrates %>% gather("year", "pop", 26:70) %>% drop_na(pop)
vertebrates$pop <- as.numeric(as.character(vertebrates$pop))

vertebrates <- vertebrates %>% group_by(id) %>% 
  filter(length(unique(year)) > 11)

vertebrates <- distinct(vertebrates)

vertebrates <- vertebrates %>% group_by(id) %>%
  dplyr::slice(-c(1:6))

vertebrates <- vertebrates  %>%  mutate(scalepop = (pop-min(pop))/(max(pop)-min(pop))) %>%
  dplyr::select (-pop) %>% drop_na(scalepop) %>%
  ungroup()

vertebrates$scalepop <- vertebrates$scalepop + 0.01
vertebrates$year <- as.factor(as.character(vertebrates$year))


vertebrates <- vertebrates %>% group_by(id) %>%
  spread(year, scalepop)

# Count years in each time-series
#vertebrates$Years <- NA

for(i in 1:nrow(vertebrates)){
  vertebrates$Years[i] <- sum(is.na(vertebrates[i, 26:64]) == FALSE)
}

# Remove populations with less than 5 years of data
#vertebrates_trim <- filter(vertebrates, Years > 5)

vertebrates_trim <- vertebrates

# Remove abundance in each year, just to simplify the dataframe
# Pop change values (mu) will be added to the simplified dataframe later
Vertebrate_data <- vertebrates_trim[,-c(26:65)]

# Compile all time-series into a list to estimate state-space population trends##
Vert <- list()
for(i in 1:nrow(vertebrates_trim)){
  data <- vertebrates_trim[i, 26:64]
  Year <- colnames(data)[is.na(data) == FALSE]
  data <- data[Year]
  N <- as.vector(t(data))
  Vert[[i]] <- data.frame(vertebrates_trim$id[i], Year, N)
}

# Define functions ----

# Functions from Humbert et al. 2009
# Define REML log-likelihoods
#  REML objective function "negloglike.reml" is negative of log-likelihood
#  for second differences of the log-scale observations.  The REML objective
#  function uses equations A18-A25 from Humbert et al. (2009).  The three
#  function arguments are:  theta, vector of parameters (transformed to the
#  real line), yt, vector of time series observations (log scale), and
#  tt, vector of observation times.  Function performs the differencing.
negloglike.reml = function(theta, yt, tt)
{
  sigsq = exp(theta[1]);         #  Constrains ssq > 0.
  tausq = exp(theta[2]);         #  Constrains tsq > 0.
  q = length(yt) - 1;
  qp1 = q + 1;
  vx = matrix(0, qp1, qp1);
  for (ti in 1:q)
  {
    vx[(ti + 1):qp1,(ti + 1):qp1] = matrix(1, 1, (qp1 - ti))*tt[ti + 1];
  }
  Sigma.mat = sigsq*vx;
  Itausq = matrix(rep(0, (qp1*qp1)), nrow = q + 1, ncol = q + 1);
  diag(Itausq) = rep(tausq, q + 1);
  V = Sigma.mat + Itausq;
  ss = tt[2:qp1] - tt[1:q];
  D1mat = cbind(-diag(1/ss), matrix(0, q, 1)) + cbind(matrix(0, q, 1), diag(1/ss));
  D2mat = cbind(-diag(1, q - 1), matrix(0,q - 1,1)) +
    cbind(matrix(0, q - 1, 1), diag(1,q - 1));
  Phi.mat = D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat);
  wt = (yt[2:qp1] - yt[1:q])/ss;
  ut = wt[2:q] - wt[1:q - 1];
  ofn = (q/2)*log(2*pi) + (0.5*log(det(Phi.mat)))+
    (0.5*(ut%*%ginv(Phi.mat)%*%ut));
  return(ofn);
}

# Loop through every time-series, estimating mu from state-space (REML) models ----
for(i in 1:nrow(vertebrates_trim)){
  Observed.t <- as.numeric(Vert[[i]]$N)
  Time.t <- as.numeric(Vert[[i]]$Year)
  print(i)
  
  library(MASS);            #  Loads miscellaneous functions (ginv, etc.)
  T.t = Time.t - Time.t[1];     #  Time starts at zero.
  
  for(z in 1:length(Observed.t)){
    if(Observed.t[z] == 0){
      if(max(Observed.t %% 1) == 0){
        Observed.t <- Observed.t + 1}
      if(max(Observed.t %% 1) > 0){
        Observed.t <- Observed.t + max(Observed.t %% 1)
      }
    }
  }
  Y.t = Observed.t;        #  Log-transform the observations.
  q = length(Y.t) - 1;          #  Number of time series transitions, q.
  qp1 = q + 1;                  #  q+1 gets used a lot, too.
  S.t = T.t[2:qp1] - T.t[1:q];  #  Time intervals.
  m = rep(1, qp1);              #  Will contain Kalman means for Kalman calculations.
  v = rep(1, qp1);              #  Will contain variances for Kalman calculations.
  
  # Calculating EGOE AND EGPN estimates for use as initial values ----
  
  # The EGOE estimates
  Ybar = mean(Y.t);
  Tbar = mean(T.t);
  mu.egoe = sum((T.t - Tbar)*(Y.t - Ybar))/sum((T.t - Tbar)*(T.t - Tbar));
  x0.egoe = Ybar - mu.egoe*Tbar;
  ssq.egoe = 0;
  Yhat.egoe = x0.egoe + mu.egoe*T.t;
  tsq.egoe = sum((Y.t - Yhat.egoe)*(Y.t - Yhat.egoe))/(q - 1);
  
  # The EGPN estimates
  Ttr = sqrt(S.t);
  Ytr = (Y.t[2:qp1] - Y.t[1:q])/Ttr;
  mu.egpn = sum(Ttr*Ytr)/sum(Ttr*Ttr);
  Ytrhat = mu.egpn*Ttr;
  ssq.egpn = sum((Ytr - Ytrhat)*(Ytr - Ytrhat))/(q - 1);
  tsq.egpn = 0;
  x0.egpn = Y.t[1];
  
  # Initial values for EGSS are averages of EGOE and EGPN values 
  ssq0 = ssq.egpn/2;            #  For ML and REML
  tsq0 = tsq.egoe/2;            #  For ML and REML
  
  # To set different initial values for iterations, enter manually a value
  #   after the equal sign of the concern parameter instead of the
  #   automatically generated value. Then run again the line and the program
  #   section 5 below.
  #   Initial values near the EGOE and EGPN models are good for exploring
  #   possible alternative local maxima. The values which produce the largest
  #   log-likelihood should be used. To see the log-likelihood for the REML
  #   estimates type:
  #   EGSSreml$value[1];
  #   See Dennis et al. 2006 for more details.
  
  # Calculate ML & REML parameter estimates ----
  
  # The REML estimates.
  EGSSreml = optim(par = c(ssq0, tsq0),
                   negloglike.reml, NULL, method = "Nelder-Mead", yt = Y.t, tt = T.t);
  params.reml = c(exp(EGSSreml$par[1]), exp(EGSSreml$par[2]))
  ssq.reml = params.reml[1];   
  #  These are the REML estimates.
  tsq.reml = params.reml[2];   
  vx = matrix(0, qp1, qp1);
  for (ti in 1:q)
  {
    vx[(ti + 1):qp1, (ti + 1):qp1] = matrix(1, 1, (qp1 - ti))*T.t[ti + 1];
  }
  Sigma.mat = ssq.reml*vx;
  Itausq = matrix(rep(0,(qp1*qp1)), nrow = q + 1, ncol = q + 1);
  diag(Itausq) = rep(tsq.reml, q + 1);
  V = Sigma.mat + Itausq;
  D1mat = cbind(-diag(1/S.t), matrix(0, q, 1)) + cbind(matrix(0, q, 1), diag(1/S.t));
  V1mat = D1mat%*%V%*%t(D1mat);
  W.t = (Y.t[2:qp1] - Y.t[1:q])/S.t;
  j1 = matrix(1, q, 1);
  V1inv = ginv(V1mat);
  mu.reml = (t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1);
  j = matrix(1, qp1, 1);
  Vinv = ginv(V);
  x0.reml = (t(j)%*%Vinv%*%(Y.t-as.vector(mu.reml)*T.t))/(t(j)%*%Vinv%*%j);  # Gergana added as.vector()
  Var_mu.reml = 1/(t(j1)%*%V1inv%*%j1);         #  Variance of mu
  mu_hi.reml = mu.reml + 1.96*sqrt(Var_mu.reml);  #  95% CI for mu
  mu_lo.reml = mu.reml - 1.96*sqrt(Var_mu.reml);
  
  #  Calculate estimated population sizes for EGSS model with Kalman filter, for plotting.
  #  Choose REML estimates here for calculating model values
  mu = mu.reml;  ssq = ssq.reml;  tsq = tsq.reml;  x0 = x0.reml;
  m[1] = x0;       
  #  Initial mean of Y(t).
  v[1] = tsq;      
  #  Initial variance of Y(t).
  for (ti in 1:q)   #  Loop to generate estimated population abundances
  {                 #    using Kalman filter (see equations 6 & 7, Dennis et al. (2006)).
    m[ti + 1] = mu + (m[ti] + ((v[ti] - tsq)/v[ti])*(Y.t[ti] - m[ti]));
    v[ti + 1] = tsq*((v[ti] - tsq)/v[ti]) + ssq + tsq;
  }
  Predict.t = exp(m + ((v - tsq)/v)*(Y.t - m))
  Vert[[i]]$Pred.N <- Predict.t
  
  #  The following statement calculates exp{E[X(t) | Y(t), Y(t-1),...,Y(0)]};
  #  Print the parameter estimates
  parms.egoe = c(mu.egoe, ssq.egoe, tsq.egoe, x0.egoe); #  Collect for printing
  parms.egpn = c(mu.egpn, ssq.egpn, tsq.egpn, x0.egpn);
  parms.reml = c(mu.reml, ssq.reml, tsq.reml, x0.reml); 
  names = c("mu", "ssq", "tsq", "x0");             
  types = c("EGOE","EGPN","EGSS-ML","EGSS-REML");    
  
  # Add to dataframe
  Vertebrate_data[i,26] <- parms.reml[1]
  Vertebrate_data[i,27] <- mu_lo.reml[1,1]
  Vertebrate_data[i,28] <- mu_hi.reml[1,1]
  Vertebrate_data[i,29] <- parms.reml[2]
  Vertebrate_data[i,30] <- parms.reml[3]
  Vertebrate_data[i,31] <- parms.reml[4]
  Vertebrate_data[i,32] <- min(Time.t)
  Vertebrate_data[i,33] <- max(Time.t)
}
colnames(Vertebrate_data)[26:33]<- c("mu", "lCI", "uCI", "sigma^2",
                                     "tau^2", "x0", "Start", "End")

# Save data frame ----
Vertebrate_data_left_trunc <- Vertebrate_data
write.csv(Vertebrate_data_left_trunc, "data/output/global_mus_left_trunc.csv")