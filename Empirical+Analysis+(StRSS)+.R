data <- read.csv(file.choose())# Brows your CSV data file
attach(data)
N <- length(X)  # Population size
X <- data[,2]  # column consist of Auxiliary variable
Y <- data[,1]  # target your Study variable acordingly
data <- data.frame(X, Y)

# Stratify using K-means on X and Y
kmeans_result <- kmeans(data, centers = 4)
data$stratum <- as.factor(kmeans_result$cluster)
strata <- split(data, data$stratum)
L <- 4  # Number of strata
P_h <- 1/L  # Equal stratum weights

# Function to perform RSS within a stratum
rss_sample <- function(stratum_data, m_h = 5, r = 4) {
  ranked_data <- list()
  for (j in 1:r) {
    pool <- stratum_data[sample(nrow(stratum_data), m_h^2), ]
    pool <- pool[order(pool$X), ]  # Rank by X
    # Select the i-th ranked unit in each set (using lapply)
    selected <- lapply(1:m_h, function(i) pool[(i-1)*m_h + i, ])
    ranked_data[[j]] <- do.call(rbind, selected)
  }
  return(do.call(rbind, ranked_data))  # Returns a data frame
}

# Apply RSS to each stratum
rss_samples <- lapply(strata, rss_sample)
compute_stratum_params <- function(stratum_data, rss_sample_df) {
  Y_h <- mean(stratum_data$Y)
  X_h <- mean(stratum_data$X)
  S_yh2 <- var(stratum_data$Y)
  S_xh2 <- var(stratum_data$X)
  S_yhxh <- cov(stratum_data$Y, stratum_data$X)
  
  m_h <- 5; r <- 4
  gamma_h <- 1 / (m_h * r)
  
  # Extract Y and X from the RSS sample (now a data frame)
  mu_yhi <- sapply(1:m_h, function(i) mean(rss_sample_df$Y[seq(i, nrow(rss_sample_df), by = m_h)]))
  mu_xhi <- sapply(1:m_h, function(i) mean(rss_sample_df$X[seq(i, nrow(rss_sample_df), by = m_h)]))
  
  tau_yhi <- mu_yhi - Y_h
  tau_xhi <- mu_xhi - X_h
  tau_yhxh <- tau_yhi * tau_xhi
  
  W_yh2 <- sum(tau_yhi^2) / (m_h^2 * r * Y_h^2)
  W_xh2 <- sum(tau_xhi^2) / (m_h^2 * r * X_h^2)
  W_yhxh <- sum(tau_yhxh) / (m_h^2 * r * X_h * Y_h)
  
  C_yh <- sqrt(S_yh2) / Y_h
  C_xh <- sqrt(S_xh2) / X_h
  
  Q_yh2 <- gamma_h * C_yh^2 - W_yh2
  R_xh2 <- gamma_h * C_xh^2 - W_xh2
  S_yhxh_term <- gamma_h * (S_yhxh / (Y_h * X_h)) - W_yhxh
  
  return(list(
    Y_h = Y_h, X_h = X_h, Q_yh2 = Q_yh2, R_xh2 = R_xh2, S_yhxh = S_yhxh_term,
    W_yh2 = W_yh2, W_xh2 = W_xh2, W_yhxh = W_yhxh
  ))
}

# Compute parameters for each stratum using Map
stratum_params <- Map(compute_stratum_params, strata, rss_samples)
names(stratum_params) <- paste0("Stratum_", 1:L)

# Calculate MSEs and PREs
var_classical <- sum(sapply(stratum_params, \(p) P_h^2 * p$Y_h^2 * p$Q_yh2))

alpha1_opt <- sapply(stratum_params, \(p) 1 / (1 + p$Q_yh2 - p$S_yhxh / p$R_xh2))
mse_Ty1 <- sum(sapply(1:L, \(h) P_h^2 * stratum_params[[h]]$Y_h^2 * (1 - alpha1_opt[h]) / 20))

P2 <- sapply(stratum_params, \(p) 1 + p$S_yhxh/2 - (p$S_yhxh^2)/(2 * p$R_xh2))
Q2 <- sapply(stratum_params, \(p) 1 + p$Q_yh2 + p$S_yhxh - 2 * (p$S_yhxh^2)/p$R_xh2)
mse_Ty2 <- sum(sapply(1:L, \(h) P_h^2 * stratum_params[[h]]$Y_h^2 * (1 - (P2[h]^2 / Q2[h])) / 20))

mse_Ty3 <- sum(sapply(1:L, \(h) P_h^2 * stratum_params[[h]]$Y_h^2 * (1 - alpha1_opt[h]) / 20))

# Proposed Estimator yStP1
mse_StP1 <- sum(sapply(stratum_params, \(p) {
  numerator <- p$R_xh2^2 * (2 * p$R_xh2 + p$S_yhxh)
  denominator <- 4 * (4 * p$R_xh2^3 + (5 * p$Q_yh2 + 4 * p$S_yhxh) * p$R_xh2^2 - 
                        4 * (p$Q_yh2 + p$S_yhxh^2) * p$R_xh2 + 4 * p$S_yhxh^2)
  P_h^2 * numerator / denominator
}))

# Compute PREs
pre <- \(mse) (var_classical / mse) * 100

results <- data.frame(
  Estimator = c("Classical", "T_{y1}^c", "T_{y2}^c", "T_{y3}^c", "yStP1"),
  #MSE = c(var_classical, mse_Ty1, mse_Ty2, mse_Ty3, mse_StP1)
  PRE = pre(c(var_classical, mse_Ty1, mse_Ty2, mse_Ty3, mse_StP1))
)

print(results)
