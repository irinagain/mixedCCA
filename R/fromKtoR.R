fromKtoR <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  # If this is just 1 variable, then correlation is automatically 1
  if (d1 == 1){return(as.matrix(1))}

  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else {
    # Select bridge function based on the type of variables
    bridge <- bridge_select(type1 = type, type2 = type)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d1)

    for(i in 1:(d1 - 1)) {
      for(j in (i + 1):d1){
        # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
        f1 <- function(r)(bridge(r, zratio1 = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
        if(op == 100){
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          hatR[i, j] <- hatR[j, i] <- NA
        }else {
          hatR[i, j] <- hatR[j, i] <- unlist(op)
        }
      }
    }
  }
  return(hatR)
}

fromKtoR_mixed <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", tol = 1e-3) {

  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)

  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # Select bridge function based on the type of variables
    bridge <- bridge_select(type1 = type1, type2 = type2)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d2)

    for(i in 1:d1) {
      for(j in 1:d2){
        # Optimize with given bridge
        f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
        op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          hatR[i,j] <- NA
        } else {
          hatR[i,j] <- unlist(op)
        }
      }
    }
  }

  return(hatR)
}

