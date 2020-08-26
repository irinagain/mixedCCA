

############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################

#' @importFrom chebpol ipol
NULL

############################################################################################
# Cutoff criteria based on the combination of variable types
############################################################################################

cutoff_bb <- function(zratio1, zratio2){ 0.9*2*pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2)) }
cutoff_tt <- function(zratio1, zratio2){ 0.9*(1-pmax(zratio1, zratio2)^2) }
cutoff_tc <- function(zratio1, zratio2 = NULL){ 0.9*(1 - zratio1^2) }
cutoff_ct <- function(zratio1 = NULL, zratio2){ 0.9*(1 - zratio2^2) }
cutoff_bc <- function(zratio1, zratio2 = NULL){ 0.9*2*zratio1*(1 - zratio1) }
cutoff_cb <- function(zratio1 = NULL, zratio2){ 0.9*2*zratio2*(1 - zratio2) }
cutoff_tb <- function(zratio1, zratio2){ 0.9*pmin(1-zratio1^2, 2*zratio2*(1 - zratio2)) }
cutoff_bt <- function(zratio1, zratio2){ 0.9*pmin(2*zratio1*(1 - zratio1), 1-zratio2^2) }


cutoff_select <- function(type1, type2){
  if (type1 == "binary" & type2 == "binary") {
    cutoff_select <- cutoff_bb
  } else if (type1 == "trunc" & type2 == "trunc") {
    cutoff_select <- cutoff_tt
  } else if (type1 == "trunc" & type2 == "continuous") {
    cutoff_select <- cutoff_tc
  } else if (type1 == "continuous" & type2 == "trunc") {
    cutoff_select <- cutoff_ct
  } else if (type1 == "binary" & type2 == "continuous") {
    cutoff_select <- cutoff_bc
  } else if (type1 == "continuous" & type2 == "binary") {
    cutoff_select <- cutoff_cb
  } else if (type1 == "trunc" & type2 == "binary") {
    cutoff_select <- cutoff_tb
  } else if (type1 == "binary" & type2 == "trunc") {
    cutoff_select <- cutoff_bt
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}



############################################################################################
# Select which bridge inverse function based on the combinatino of variable types
############################################################################################

bridgeInv_select <- function(type1, type2) {
  if (type1 == "binary" & type2 == "binary") { bridgeInv_select <- bridgeInv_bb
  } else if (type1 == "trunc" & type2 == "trunc") { bridgeInv_select <- bridgeInv_tt
  } else if (type1 == "trunc" & type2 == "continuous") { bridgeInv_select <- bridgeInv_tc
  } else if (type1 == "continuous" & type2 == "trunc") { bridgeInv_select <- bridgeInv_ct
  } else if (type1 == "binary" & type2 == "continuous") { bridgeInv_select <- bridgeInv_bc
  } else if (type1 == "continuous" & type2 == "binary") { bridgeInv_select <- bridgeInv_cb
  } else if (type1 == "trunc" & type2 == "binary") { bridgeInv_select <- bridgeInv_tb
  } else if (type1 == "binary" & type2 == "trunc") { bridgeInv_select <- bridgeInv_bt
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}


# wrapper function to make matrix and vector input available and allow the different order of the combination too.
bridgeInv_tc <- function(tau, zratio1, zratio2 = NULL){
  tau <- as.matrix(tau)
  zratio1 <- as.matrix(zratio1, nrow = 1)
  out <- TCipol(rbind(c(tau), rep(zratio1, ncol(tau))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(out)
}

bridgeInv_ct <- function(tau, zratio1 = NULL, zratio2){
  zratio1 <- zratio2
  tau <- t(as.matrix(tau))
  # above is the only difference with "bridgeInv_tc" to switch the order of the variable, transposed the matrix "tau".
  # and the final returned output should be also transposed back again.

  zratio1 <- as.matrix(zratio1, nrow = 1)
  out <- TCipol(rbind(c(tau), rep(zratio1, ncol(tau))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(t(out))
}


# wrapper function to use matrix and vectors as inputs.
bridgeInv_tt <- function(tau, zratio1, zratio2){
  tau <- as.matrix(tau)
  zratio1 <- as.matrix(zratio1, nrow = 1)
  zratio2 <- as.matrix(zratio2, nrow = 1)
  out <- TTipol(rbind(c(tau), rep(zratio1, length(zratio2)), rep(zratio2, each = length(zratio1))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(out)
}


# wrapper function to use matrix and vectors as inputs.
bridgeInv_tb <- function(tau, zratio1, zratio2){
  tau <- as.matrix(tau)
  zratio1 <- as.matrix(zratio1, nrow = 1)
  zratio2 <- as.matrix(zratio2, nrow = 1)
  out <- TBipol(rbind(c(tau), rep(zratio1, length(zratio2)), rep(zratio2, each = length(zratio1))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(out)
}

bridgeInv_bt <- function(tau, zratio1, zratio2){
  tau <- t(as.matrix(tau)) # transpose the matrix to switch b and t
  zratio1 <- as.matrix(zratio2, nrow = 1)
  zratio2 <- as.matrix(zratio1, nrow = 1)
  # above are the only difference with bridgeInv_tb to switch the variable. and the final returned output should be also transposed back again.

  out <- TBipol(rbind(c(tau), rep(zratio1, length(zratio2)), rep(zratio2, each = length(zratio1))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(t(out))
}


# wrapper function to make matrix and vector input available.
bridgeInv_bc <- function(tau, zratio1, zratio2 = NULL){
  tau <- as.matrix(tau)
  zratio1 <- as.matrix(zratio1, nrow = 1)
  out <- BCipol(rbind(c(tau), rep(zratio1, ncol(tau))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(out)
}

bridgeInv_cb <- function(tau, zratio1 = NULL, zratio2){
  tau <- t(as.matrix(tau))
  zratio1 <- as.matrix(zratio2, nrow = 1)
  # above are the only difference with bridgeInv_cb to switch the variable.
  # and the final returned output should be also transposed back again.

  out <- BCipol(rbind(c(tau), rep(zratio1, ncol(tau))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(t(out))
}

# wrapper function to use matrix and vectors as inputs.
bridgeInv_bb <- function(tau, zratio1, zratio2){
  tau <- as.matrix(tau)
  zratio1 <- as.matrix(zratio1, nrow = 1)
  zratio2 <- as.matrix(zratio2, nrow = 1)
  out <- BBipol(rbind(c(tau), rep(zratio1, length(zratio2)), rep(zratio2, each = length(zratio1))))
  out <- matrix(out, nrow = nrow(tau), ncol = ncol(tau))
  return(out)
}
