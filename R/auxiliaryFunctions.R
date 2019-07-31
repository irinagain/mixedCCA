
# Calculate norm of canonical vector w. It is calculated as t(w) %*% Sigma %*% w.
# If w is all zero vector, then 0 will be returned.
w_norm <- function(W, Sigma){
  W <- as.vector(W)
  # Check first if the inputs are valid.
  if(!isSymmetric(Sigma)){
    stop("Sigma should be symmetric.")
  }
  if(length(W) != nrow(Sigma)){
    stop("non-conformable arguments. Check W and Sigma.")
  }
  wnorm <- as.numeric(crossprod(W, Sigma %*% W))
  return(wnorm)
}



# Function to normalize the canonical vector W1. The returned vector should satisfy t(W) %*% Sigma %*% W = 1.
normalizedW <- function(W, Sigma){
  W <- as.vector(W)
  # Check first if the inputs are valid.
  if(!isSymmetric(Sigma)){
    stop("Sigma should be symmetric.")
  }
  if(length(W) != nrow(Sigma)){
    stop("non-conformable arguments. Check W and Sigma.")
  }
  normW <- w_norm(W, Sigma)
  if(normW == 0){
    normalizedW <- 0
  } else {
    normalizedW <- W/sqrt(normW)
  }
  return(normalizedW)
}


# Estimate canonical correlation based on two estimates w1 and w2
# Sigma1, Sigma2, Sigma12 are all true latent correlation matrices or from out-of-sample correlation matrices.
# w1 is a vector of length p1. w2 is a vector of length p2.
# Sigma1 is a matrix of size p1 by p1. Sigma2 is a matrix of size p2 by p2. Sigma12 is a matrix of size p1 by p2.
cancorhat <- function(w1, w2, Sigma1, Sigma2, Sigma12){
  w1 <- as.vector(w1)
  w2 <- as.vector(w2)

  # Check first if the inputs are valid.
  if(length(w1) != ncol(Sigma1)){
    stop("non-conformable arguments: xcoef and Sigma1")
  }
  if(length(w2) != ncol(Sigma2)){
    stop("non-conformable arguments: ycoef and Sigma2")
  }
  if(ncol(Sigma1) != nrow(Sigma12)){
    stop("non-conformable arguments: Sigma1 and Sigma12")
  }
  if(ncol(Sigma2) != ncol(Sigma12)){
    stop("non-conformable arguments: Sigma2 and Sigma12")
  }
  # when one of the estimates w1 or w2 is all zero vector, the output will be zero.
  norm1 <- w_norm(w1, Sigma1)
  norm2 <- w_norm(w2, Sigma2)
  if (norm1 == 0 | norm2 == 0){
    output <- 0
  } else {
    output <- as.numeric(crossprod(w1, Sigma12 %*% w2))/sqrt(norm1)/sqrt(norm2)
  }
  return(output)
}
