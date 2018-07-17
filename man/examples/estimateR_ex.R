
# Generate Data
set.seed(1)
n <- 100; p1 <- 5; p2 <- 7
maxcancor <- 0.9 # true canonical correlation
Sigma1 <- autocor(p1, 0.7)
groupind <- c(rep(1, 2), rep(2, p2-2))
Sigma2 <- blockcor(0.7, groupind)
copula1 <- "exp"; copula2 <- "cube" # transformation for Gaussian copula model
mu <- rep(0, p1+p2)
type1 <- type2 <- "trunc" # X1: truncated continuous, X2: truncated continuous
c1 <- rep(0, p1); c2 <- rep(0, p2) # threshold for truncation of underlying continuous variable
trueidx1 <- c(0, 0, 1, 1, 1) # true variable indices for dataset X1
trueidx2 <- c(1, 0, 1, 0, 0, 1, 1) # true variable indices for dataset X2
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = copula1, copula2 = copula2,
                        muZ = mu,
                        type1 = type1, type2 = type2, c1 = c1, c2 = c2
)
X1 <- simdata$X1
X2 <- simdata$X2

# Estimate latent correlation matrix
estimateR(X1, type = "trunc")
estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc")
