library(MASS)
library(moments)

#### POPULATION
set.seed(123)
mu <- c(1, 2)
sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
n_pop <- 1000

#### Generate correlated normal variables
rawvars <- mvrnorm(n_pop, mu, sigma)
pvars <- pnorm(rawvars)
data <- qgamma(pvars, 1:2)
data1 <- data.frame(Y = data[,1], X = data[,2])

#######POPULATION PARAMETERS
Xb <- mean(data1$X)
Yb <- mean(data1$Y)
sdx <- sd(data1$X)
sdy <- sd(data1$Y)
bx1 <- skewness(data1$X)
bx2 <- kurtosis(data1$X)
rho <- cor(data1$X, data1$Y)
cx <- sdx/Xb
cy <- sdy/Yb

#### STRATIFICATION
H <- 3  # Number of strata
data1$stratum <- cut(data1$X, 
                     breaks = quantile(data1$X, seq(0, 1, length.out = H+1)),
                     include.lowest = TRUE)
strata <- split(data1, data1$stratum)
Nh <- sapply(strata, nrow)
Wh <- Nh/n_pop

## Summary statistics
m <- 3
; r <- 3; n <- m*r
f <- n/n_pop
lambda <- (1-f)/n
P <- 0.5
u <- rho
v <- bx1

ki <- (u*Xb)/(u*Xb + v)
J1 <- 0.5*ki + (2*P-1)
J2 <- P + 0.5*(2*P-1)*ki + 0.375*ki^2

####### Summary Statistic
Wy <- 0; Wx <- 0; Wyx <- 0
Hy <- (1/(m*r))*cy - Wy^2
Gx <- (1/(m*r))*cx - Wx^2
cyx <- rho*cy*cx
Jyx <- (1/(m*r))*cyx - Wyx

##Optimization constants
T1 <- Yb*(14*Gx^6 + 18*Gx^4*Hy^2 + 15*Hy^4*Jyx - 30*Gx^4*Hy^2 + 16*Jyx^2) /
  (4*(4*Gx^4 + 5*Gx^4*Jyx - 8*Gx^2*Hy^2 - 4*Gx^2*Jyx^2 + 4*Jyx^2))
T2 <- (Gx^4*(2*Gx^2 + Jyx)) /
  (4*(4*Gx^6 + 5*Gx^4*Hy^2 + 4*Gx^4*Jyx - 4*Gx^2*Hy^2 - 4*Gx^2*Jyx^2 + 4*Jyx^2))
T3 <- -Yb*(6*Gx^6 + 7*Gx^4*Hy^2 + 5*Gx^4*Jyx - 14*Gx^4*Hy^2 - 6*Gx^2*Jyx^2 + 8*Jyx^2) /
  (4*Xb*(4*Gx^6 + 5*Gx^4*Hy^2 + 4*Gx^4*Jyx - 4*Gx^2*Hy^2 - 4*Gx^2*Jyx^2 + 4*Jyx^2))
T4 <- (10*Gx^4*Hy^2*Yb + 9*Gx^4*Jyx*Yb + 26*Gx^2*Jyx^2*Yb + 8*Gx^2*Hy^2 - 12*Gx^2*Jyx - 8*Jyx^2) /
  (Gx^2*(4*Gx^2*Hy^2 + 5*Jyx^2))
T5 <- (-3*Jyx*(3*Gx^2*Yb - 4)) /
  (4*Yb*(4*Gx^2*Hy^2 + 5*Jyx^2))
A1 <- (cx^2*Xb)/(cx^2*Xb + n_pop)
T6 <- -(4*Gx^4*Hy^2*Yb + 7*Gx^2*Jyx^2*Yb + 8*Gx^2*Hy^2 + 4*Jyx^2) /
  (4*Xb*Gx^2*(4*Gx^2*Hy^2 + 5*Jyx^2))
T7_numer <- Gx^4*Yb*(A1^2 + 2*A1 + 1) - A1^2*Gx^2*Jyx*(2*A1 + 1) - 
  Gx^2*Jyx*Yb*(A1^2 + A1) + (16*A1^2*Yb + A1*Gx^2 + 16*A1*Yb - 8*A1 + 4*Yb + 4)*Jyx + 
  (2*A1 - 1)*Jyx^2 - 16*Hy^2*Yb
T7_denom <- (A1^4 + 2*A1^3 + A1^2)*Gx^4 - (2*A1^2 + 2*A1)*Gx^2*Jyx + 
  16*(A1^2 + A1)*Jyx - 16*Hy^2 + Jyx^2
T7 <- T7_numer / T7_denom
T8 <- 4*Jyx*(2*A1 - Yb - 1)/Yb * T7_numer
T9 <- -T7_numer / (2*Xb*T7_numer)  # Simplified form
### SIMULATION SETUP
sim_iter <- 10000
U <- 1.2
# Initialize result vectors
estimators <- list(
  SRS = numeric(sim_iter),
  Ratio = numeric(sim_iter),
  Product = numeric(sim_iter),
  BT = numeric(sim_iter),
  ET = numeric(sim_iter),
  P1 = numeric(sim_iter),
  P2 = numeric(sim_iter),
  P3 = numeric(sim_iter)
)

#### STRATIFIED RSS IMPLEMENTATION
for(i in 1:sim_iter) {
  stratum_ests <- lapply(1:H, function(h) {
    # Per-stratum operations
    data_h <- strata[[h]]
    
    # Ranked set sampling
    idx <- sample(nrow(data_h), m*r, replace = TRUE)
    X_vals <- data_h$X[idx]
    Y_vals <- data_h$Y[idx]
    
    # Matrix formation and ranking
    X_matrix <- matrix(X_vals, nrow = r, ncol = m)
    Y_matrix <- matrix(Y_vals, nrow = r, ncol = m)
    ranks <- apply(X_matrix, 2, rank)
    
    # Selection mechanism
    selected <- apply(ranks, 2, which.min)
    y_measured <- Y_matrix[cbind(selected, 1:m)]
    x_measured <- X_matrix[cbind(selected, 1:m)]
    
    list(
      yb = mean(y_measured),
      xb = mean(x_measured)
    )
  })
  
  # Combine stratum estimates
  yb <- sum(Wh * sapply(stratum_ests, function(x) x$yb))
  xb <- sum(Wh * sapply(stratum_ests, function(x) x$xb))
  
  #===== [8] ESTIMATOR CALCULATIONS =====#
  # Basic estimators
  estimators$SRS[i] <- yb
  estimators$Ratio[i] <- yb * (Xb/xb)
  estimators$Product[i] <- yb * (xb/Xb)
  estimators$BT[i] <- yb * exp((Xb - xb)/(Xb + xb))
  A <- lambda*(cy^2 + (J1^2 + 2*J2)*cx^2 - 4*J1*rho*cy*cx)
  B <- 1 + lambda*(J1^2 + 2*J2)*cx^2
  C <- lambda*(J2*cx^2 - J1*rho*cy*cx)
  D <- 1 + lambda*J2*cx^2
  E <- 1 + lambda*((J1^2 + 2*J2)*cx^2 - 2*J1*rho*cy*cx)
  t1 <- (B*C - D*E + B)/(A*B - E^2 + B)
  t2 <- Yb*(A*D - C*E + D - E)/(A*B - E^2 + B)
  estimators$ET[i] <- (t1*yb + t2) * 
    (P*(xb/Xb) + (1-P)*(Xb/xb)) * 
    exp(u*(Xb - xb)/(u*(Xb + xb) - 2*v))
  
  # Proposed estimators P1 and P3(P2 in the article)
  estimators$P1[i] <- T1 + T2*yb + T3*(Xb - xb)*exp((Xb - xb)/(Xb + xb))
 # estimators$P2[i] <- T4 + T5*yb + T6*(Xb - xb)*exp((xb - Xb)/(xb + Xb))
  estimators$P3[i] <- (T7 + T8*yb*U + T9*(Xb - xb)) * exp((Xb - xb)/(Xb + xb))
}

############ RESULTS ANALYSIS 
results <- data.frame(
  Estimator = names(estimators),
  Variance = sapply(estimators, var)
)

print(results)