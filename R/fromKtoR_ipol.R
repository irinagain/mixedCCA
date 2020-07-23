##### This is going to be approximate version of fromKtoR which is much faster using multilinear interpolation. (07-23-2020)
##### multilinear interpolation method is from ipol package.

# K: Kendall's tau matrix.
# zratio: a column vector of zero proportion values.
fromKtoR_ml <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  # K is d1 by d1 square matrix. zratio should be a vector of length d1. If the type is "continuous", no input is necessary for zratio.

  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else { # if the type is either "trunc" or "binary"

    if (is.null(zratio)){ stop ("The input for zratio is required for \"trunc\" and \"binary\" types.") }
    if (length(zratio)!=d1){ stop ("The length of zratio must match with the number of columns in K.") }

          # based on the data type, select bridgeInv and cutoff functions.
          bridgeInv <- bridgeInv_select(type1 = type, type2 = type)
          cutoff <- cutoff_select(type1 = type, type2 = type)

          # much faster multi-linear interpolation part using saved ipol function.
          hatR <- bridgeInv(K, zratio1 = zratio, zratio2 = zratio)

          # make sure the diagonal elements are all one and they are symmetric.
          diag(hatR) <- 1
          hatR <- (hatR + t(hatR))/2 # even without this step, hatR is very close to symmetric but not exactly. (symmetric within the error 1e-5)

          # check if there is any element that is outside of the safe boundary for interpolation.
          ##### only considering off-diagonal matrix.
          cutoff_matrix <- matrix(cutoff(zratio1 = rep(zratio, d1), zratio2 = rep(zratio, each = d1)), nrow = d1)
          ind_cutoff <- which(abs(K[upper.tri(K)]) > cutoff_matrix[upper.tri(cutoff_matrix)], arr.ind = TRUE)

          if (length(ind_cutoff) > 0){
            # if there is any element that is outside of the safe boundary for interpolation,
            # then calculate using original bridge inverse optimization elementwise.
            ind_cutoff_mat <- cbind(floor(sqrt(ind_cutoff*2)), ceiling(sqrt(ind_cutoff*2)))
            bridge <- bridge_select(type1 = type, type2 = type)

            for(rind in 1:length(ind_cutoff)){
              i <- ind_cutoff_mat[rind, 1]
              j <- ind_cutoff_mat[rind, 2]
              if(type == "trunc" & sum(zratio[c(i, j)]) < 1e-6){
                hatR[i, j] <- hatR[j, i] <- sin(pi/2 * K[i, j])
              } else {
                f1 <- function(r)(bridge(r, zratio1 = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
                op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
                if(op == 100) {
                  warning("Check pairs bewteen variable ", i, " and variable ", j, "\n")
                } else {
                  hatR[i, j] <- hatR[j, i] <- unlist(op)
                }
              }
            }
          } # done with for the pairs that are outside of the safe boundary for multi-linear interpolation.
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
  # K12 is d1 by d2 matrix.
  # zratio1 should be a vector of length d1. zratio2 should be a vector of length d2.
  # If the both types are "continuous", no input is necessary for zratio1 and zratio2.

  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # if the case is either CT, TC, TT, BC, BB or TB.
    if (type1 != "continuous"){
      if (is.null(zratio1)){ stop ("The input for zratio1 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio1)!=d1){ stop ("The length of zratio1 must match with the number of rows in K12.") }
    }
    if (type2 != "continuous"){
      if (is.null(zratio2)){ stop ("The input for zratio2 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio2)!=d2){ stop ("The length of zratio2 must match with the number of columns in K12.") }
    }

    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type1, type2 = type2)
    cutoff <- cutoff_select(type1 = type1, type2 = type2)

    # much faster multi-linear interpolation part using saved ipol function.
    hatR <- bridgeInv(K12, zratio1 = zratio1, zratio2 = zratio2)
    # Here, there is no step to make sure about the diagonal elements and symmetrize.

    # check if there is any element that is outside of the safe boundary for interpolation.
    cutoff_matrix <- matrix(cutoff(zratio1 = rep(zratio1, d2), zratio2 = rep(zratio2, each = d1)), nrow = d1)
    ind_cutoff <- which(abs(K12) > cutoff_matrix, arr.ind = TRUE)

    if (nrow(ind_cutoff) > 0){
      bridge <- bridge_select(type1 = type1, type2 = type2)

      for(rind in 1:nrow(ind_cutoff)){
        # if there is any element that is outside of the safe boundary for interpolation,
        # then calculate using original bridge inverse optimization elementwise.
        i <- ind_cutoff[rind, 1]
        j <- ind_cutoff[rind, 2]
        if(sum(c(type1, type2) %in% c("trunc", "continuous")) == 2 & sum(zratio1[i], zratio2[j]) < 1e-6){
          hatR[i, j] <- sin(pi/2 * K12[i, j])
        } else {
          f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            warning("Check pairs bewteen variable ", i, " and variable ", j, "\n")
          } else {
            hatR[i, j] <- unlist(op)
          }
        }
      }
    } # done with for the pairs that are outside of the safe boundary for multi-linear interpolation.
    }
  return(hatR)
}
