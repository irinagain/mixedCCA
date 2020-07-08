## code to prepare `DATASET` dataset goes here

# usethis::use_data("DATASET")


############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################


############################################################################################
# For TC case
############################################################################################

load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tc_0405.Rda")

# grid values that used to create precomputed values.
# d1 <- log10(seq(1, 10^0.99, length = 50))
# tau <- seq(-0.99, 0.99, by = 0.005)

# create computed values (in matrix) and grid (in list) for ipol function.
value <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)
grid <- list(tau, d1) # the length of list should be the same as the kinds of inputs.

interp_multilin <- chebpol::ipol(value, grid = grid, method = "multilin")

# create input values for ipol
TCvalue <- matrix(unlist(gridTCinv), ncol = length(d1), byrow = FALSE)

# create grid input for ipol
TCipolgrid <- list(tau, d1)

# interpolation
TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")

############################################################################################
# For TT case
############################################################################################

load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tt_0405.Rda")

TTvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      TTvalue[k, i, j] <- gridTTinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# grid values that used to create precomputed values.
# d2 <- log10(seq(1, 10^0.99, length = 50))

# create grid input for ipol
TTipolgrid <- list(tau, d1, d2)

# interpolation.
TTipol <- chebpol::ipol(TTvalue, grid = TTipolgrid, method = "multilin")

############################################################################################
# For TB case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/tb_0405.Rda")

TBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      TBvalue[k, i, j] <- gridTBinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# create grid input for ipol
TBipolgrid <- list(tau, d1, d2)

# interpolation.
TBipol <- chebpol::ipol(TBvalue, grid = TBipolgrid, method = "multilin")

############################################################################################
# For BC case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bc_0405.Rda")

# create input values for ipol
BCvalue <- matrix(unlist(gridBCinv), ncol = length(d1), byrow = FALSE)

# create grid input for ipol
BCipolgrid <- list(tau, d1)

# interpolation
BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")

############################################################################################
# For BB case
############################################################################################
load("~/Dropbox/TAMU/Irina/mixedCCAfast/R1_PrecomputedResults/bb_0405.Rda")

BBvalue <- array(NA, dim = c(length(tau), length(d1), length(d2)))
for (i in 1:length(d1)){
  for ( j in 1:length(d2)){
    for ( k in 1:length(tau)){
      BBvalue[k, i, j] <- gridBBinv[[length(d2)*(i - 1) + j]][k]
    }
  }
}

# create grid input for ipol
BBipolgrid <- list(tau, d1, d2)

# interpolation.
BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")

usethis::use_data(TCipol, TTipol, TBipol, BCipol, BBipol, internal = TRUE, overwrite = TRUE, compress = "xz")

