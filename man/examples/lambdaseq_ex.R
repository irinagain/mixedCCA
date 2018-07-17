# Generate Data
n <- 100; p1 <- 5; p2 <- 7
maxcancor <- 0.9
Sigma1 <- autocor(p1, 0.7)
groupind <- c(rep(1, 2), rep(2, p2-2))
Sigma2 <- blockcor(0.7, groupind)
copula1 <- "exp"; copula2 <- "cube"
mu <- rep(0, p1+p2)
type1 <- type2 <- "trunc"
c1 <- rep(0, p1); c2 <- rep(0, p2)
trueidx1 <- c(0, 0, 1, 1, 1)
trueidx2 <- c(1, 0, 1, 0, 0, 1, 1)
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = copula1, copula2 = copula2,
                        muZ = mu,
                        type1 = type1, type2 = type2, c1 = c1, c2 = c2
)

X1 <- simdata$X1
X2 <- simdata$X2

R <- estimateR_mixed(X1, X2, type1 = type1, type2 = type2)
R1 <- as.matrix(R$R1)
R2 <- as.matrix(R$R2)
R12 <- as.matrix(R$R12)

# Tuning parameter grid
lambda_seq <- lambdaseq_generate(nlambda = 10, Sigma1 = R$R1, Sigma2 = R$R2, Sigma12 = R$R12)
