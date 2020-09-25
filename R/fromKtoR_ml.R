##### This is going to be approximate version of fromKtoR which is much faster using multilinear interpolation. (07-23-2020)
##### multilinear interpolation method is from ipol in chebpol package.

# K: Kendall's tau matrix.
# zratio: a column vector of zero proportion values.
fromKtoR_ml <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  # If this is just 1 variable, then correlation is automatically 1
  if (d1 == 1){return(as.matrix(1))}

  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else { # if the type is either "trunc" or "binary"
    upperR <- c(upper.tri(K)) # length p^2 of true/false with true corresponding to upper.tri
    hatRupper <- rep(NA, sum(upperR)) # length p(p-1)/2
    Kupper <- c(K[upperR]) # upper triangle of K matrix

    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type, type2 = type)
    cutoff <- cutoff_select(type1 = type, type2 = type)

    # check if there is any element that is outside of the safe boundary for interpolation.
    zratio1vec = rep(zratio, d1)[upperR] # length p(p-1)/2
    zratio2vec = rep(zratio, each = d1)[upperR] # length p(p-1)/2
    ind_cutoff <- which(abs(Kupper) > cutoff(zratio1vec, zratio2vec))

    if (length(ind_cutoff) == 0){
      # multi-linear interpolation part using saved ipol function.
      hatRupper <- bridgeInv(Kupper, zratio1 = zratio1vec, zratio2 = zratio2vec)
    }else{
      # Interpolate only those elements that are inside
      hatRupper[-ind_cutoff] <- bridgeInv(Kupper[-ind_cutoff], zratio1 = zratio1vec[-ind_cutoff], zratio2 = zratio2vec[-ind_cutoff])

      # Apply original method to the elements outside
      bridge <- bridge_select(type1 = type, type2 = type)

      for(ind in ind_cutoff){
        f1 <- function(r)(bridge(r, zratio1 = zratio1vec[ind], zratio2 = zratio2vec[ind]) - Kupper[ind])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          hatRupper[ind] <- NA
        } else {
          hatRupper[ind] <- unlist(op)
        }
      }
    }

    # Get upperR into hatR
    hatR <- matrix(0, d1, d1)
    hatR[upperR] <- hatRupper
    hatR <- hatR + t(hatR)
    diag(hatR) <- 1
  }
  return(hatR)
}

# K12: Kendall's tau matrix.
# zratio1: a vector of zero proportion values for row variables. The length should match with nrow of K12.
# zratio2: a vector of zero proportion values for column variables. The length should match with ncol of K12.
fromKtoR_ml_mixed <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", tol = 1e-3) {

  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)

  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # if the case is either CT, TC, TT, BC, BB or TB.
    hatR <- matrix(NA, d1, d2)

    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type1, type2 = type2)
    cutoff <- cutoff_select(type1 = type1, type2 = type2)

    # check if there is any element that is outside of the safe boundary for interpolation.
    zratio1vec = rep(zratio1, d2)
    zratio2vec = rep(zratio2, each = d1)
    ind_cutoff <- which(abs(K12) > cutoff(zratio1vec, zratio2vec))

    # much faster multi-linear interpolation part using saved ipol function.
    if (length(ind_cutoff) == 0){
      # Interpolate all the elements
      hatR <- matrix(bridgeInv(c(K12), zratio1 = zratio1vec, zratio2 = zratio2vec), d1, d2)
    }else{
      # Interpolate only those elements that are inside
      hatR[-ind_cutoff] <- bridgeInv(c(K12[-ind_cutoff]), zratio1 = zratio1vec[-ind_cutoff], zratio2 = zratio2vec[-ind_cutoff])

      # Apply original method to the elements outside
      bridge <- bridge_select(type1 = type1, type2 = type2)

      for(ind in ind_cutoff){
        f1 <- function(r)(bridge(r, zratio1 = zratio1vec[ind], zratio2 = zratio2vec[ind]) - K12[ind])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          hatRupper[ind] <- NA
        } else {
          hatR[ind] <- unlist(op)
        }
      }
    }
   # done with for the pairs that are outside of the safe boundary for multi-linear interpolation.
  }

  return(hatR)
}
