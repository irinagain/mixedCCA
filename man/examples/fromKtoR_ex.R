set.seed(1)

n <- 100 # sample size
r <- 0.8 # true correlation

# Generate one dataset X: truncated continuous
Z <- mvrnorm(n, mu = c(0, 1, 0), Sigma = matrix(c(1, r, r, r, 1, r, r, r, 1), nrow = 3))
X <- ifelse(Z > 0, Z, 0)
zratio <- colMeans(X==0)
K <- Kendall_matrix(X)
# Estimate latent correlation matrix based on K.
fromKtoR(K, zratio = zratio, type = "trunc")

# Generate two dataset with different types.
# X1: truncated continuous (n by p1), X2: continuous (n by p2).
p1 <- 2; p2 <- 3; p <- p1 + p2
JSigma <- matrix(r, nrow=p, ncol=p); diag(JSigma) <- 1
Z <- mvrnorm(n, mu = c(0, 1, 0, 1, 0), Sigma = JSigma)
X1 <- Z[, 1:2]
X1[Z[, 1:2] < 0] <- 0
X2 <- Z[, 3:5]
zratio1 <- colMeans(X1==0)
bridge <- bridge_select(type1 = "trunc", type2 = "continuous")
K12 <- Kendall_matrix(X1, X2)
# Estimate latent correlation based on K12.
fromKtoR_mixed(K12, zratio1 = zratio1, type1 = "trunc", type2 = "continuous")

# X1: truncated continuous, X2: truncated continuous
X2[Z[, 3:5] < 0] <- 0
zratio2 <- colMeans(X2==0)
bridge <- bridge_select(type1 = "trunc", type2 = "trunc")
K12 <- Kendall_matrix(X1, X2)
# Estimate latent correlation based on K12.
fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc")

# X1: truncated continuous, X2: binary
X2[X2 != 0] <- 1
bridge <- bridge_select(type1 = "trunc", type2 = "binary")
K12 <- Kendall_matrix(X1, X2)
# Estimate latent correlation based on K12.
fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary")

