fromKtoR <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  de1 <- NULL
  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else {
    if (is.null(zratio)){ stop ("The input for zratio is required for \"trunc\" and \"binary\" types.") }
    if (length(zratio)!=d1){ stop ("The length of zratio must match with the number of columns in K.") }
    bridge <- bridge_select(type1 = type, type2 = type)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d1)

    if (d1 > 1){ # for matrix K
      for(i in 1:(d1-1)) {
        for(j in (i+1):d1){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i, j] <- hatR[j, i] <-0
          }else {
            hatR[i, j] <- hatR[j, i] <- unlist(op)
          }
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
  ###################################################################
  de1 <- de2 <- NULL

  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    if (type1 != "continuous"){
      if (is.null(zratio1)){ stop ("The input for zratio1 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio1)!=d1){ stop ("The length of zratio1 must match with the number of rows in K12.") }
    }
    if (type2 != "continuous"){
      if (is.null(zratio2)){ stop ("The input for zratio2 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio2)!=d2){ stop ("The length of zratio2 must match with the number of columns in K12.") }
    }

    bridge <- bridge_select(type1 = type1, type2 = type2)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d2)

    if ( d1 == 1 & d2 == 1 ){ # for scalar K: one Kendall tau value
      i <- j <- 1
      f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        hatR[i, j] <- hatR[j, i] <- 0
      } else {
        hatR[i, j] <- hatR[j, i] <- unlist(op)
      }
    } else {
      for(i in 1:d1) {
        for(j in 1:d2){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i,j] <- 0
          } else {
            hatR[i,j] <- unlist(op)
          }
        }
      }
    }
  }
  return(hatR)
}

