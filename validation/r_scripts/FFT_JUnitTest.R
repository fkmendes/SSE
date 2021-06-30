# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# PropagatesQuaSSE
# (1) testPropagateTimeOneChQuaSSETest
# (2) testMakeNormalKernInPlaceAndFFT
# (3) testPropagateXOneChQuaSSETest


library(diversitree)

## First, diversitree's functions

fftR.propagate.t <- function(vars, lambda, mu, dt, ndat) {
  i <- seq_len(ndat) # creates a seq that starts with 1 and ends in ndat, with increments of 1
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) /
    (mu + z*(e0 - 1)*lambda - lambda*e0)
  dd <- (z * r * r)/(z * lambda - mu + (1-z)*lambda*e0)^2
  ##print(paste("numerator", (z * r * r)))
  ##print(paste("denominator", (z * lambda - mu + (1-z)*lambda*e0)))
  ##print(paste("numerator/denominator", dd))
  ##vars[i,-1] <- dd * d0
  vars
}

normalise <- function(x) x / sum(x)

fftR.make.kern <- function(mean, sd, nx, dx, nkl, nkr) {
  kern <- rep(0, nx)
  xkern <- (-nkl:nkr)*dx
  ikern <- c((nx - nkl + 1):nx, 1:(nkr + 1))
  ##print("xkern=")
  ##print(xkern)
  ##print(paste0("mean=", mean, " sd=", sd))
  ##print("dnorm=")
  ##print(dnorm(xkern, mean, sd))
  kern[ikern] <- normalise(dnorm(xkern, mean, sd))
  kern
}

ifft <- function(x) fft(x, inverse=TRUE)

fftR.propagate.x <- function(vars, nx, fy, nkl, nkr) {
  vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft))/nx
  ndat <- nx - (nkl + 1 + nkr)
  i.prev.l <- 1:nkl
  i.prev.r <- (ndat-nkr+1):ndat
  i.zero <- (ndat+1):nx
  vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
  vars.out[i.zero,] <- 0
  vars.out[vars.out < 0] <- 0
  vars.out
}


## Now some input values
vars <- cbind(rep(0.0001,48),
              c(8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688))
# vars <- data.frame(0.0, 1.0) # E and D
lambda <- 1.0
mu <- 0.5
dt <- 0.01
ndat <- 40 # just one useful quant ch bin

### (1) testPropagateTimeOneChQuaSSETest

#### eAndD

res <- fftR.propagate.t(vars, lambda, mu, dt, ndat)

#                 E         D
#  [1,] 0.005061285  8.383508
#  [2,] 0.005061285  9.953882
#  [3,] 0.005061285  8.954895
#  [4,] 0.005061285 11.671936
#  [5,] 0.005061285  8.386246
#  [6,] 0.005061285  9.769374
#  [7,] 0.005061285  8.817404
#  [8,] 0.005061285  9.168019
#  [9,] 0.005061285 10.913191
# [10,] 0.005061285 10.604198

### (2) testMakeNormalKernInPlaceAndFFT

drift <- 0.0
diffusion <- 0.001
nkl <- nkr <- 4
nx <- 48
dx <- 0.001 # (dx * nx) gives the range of x values we're considering

# '-' in -dt because we're going back in time
kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion), nx, dx, nkl, nkr)

plot(kern)

## x: the substitution rate
## y: gives you the probability of changing from x to (x + delta.x) over time dt,
## with delta.x being (dx * nkr) on the right-hand side of x, and (dx * nkl) on
## the left-hand side of x
##
## (later, y will be fourier-transformed into fy)

## worked example
##
## min x = 0.0
## max x = 0.1 = 1/10
## start x at rate = 0.01 = 1/100
## dx = 0.001
## nx = 100
## nkl = nkr = 5
##
## x is at 10th bin, so we ignore anything beyond nkl and nkr bins away from 10th, meaning
## x will influence from 5th - 15th bins
##
## the first bin gives the probability that x stays where it is
## the 2nd - 6th bin give the probability x changes to 11th - 15th bins
## the 95th - 100th bins give the probability x changes to the 9th - 5th

fy <- fft(kern)
ify <- ifft(fy) # note that this is NOT what we do in the likelihood calculation

Re(fy)

# [1]  1.000000000  0.956873917  0.835109753  0.655887828  0.449367562
# [5]  0.248270676  0.081152021 -0.033081907 -0.088189035 -0.090559850

Re(ify)

# [1]  7.149417e+00  6.800735e+00  5.853447e+00  4.558669e+00  3.212440e+00
# [6]  5.957485e-16  4.549779e-16  1.366926e-15  2.766926e-16  3.584917e-17



### (3) testPropagateXOneChQuaSSETest

vars <- cbind(rep(0.0001,48),
              c(8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688))

Re(apply(apply(vars, 2, fft) * fy, 2, ifft)) / nx

##            E         D
##   [1,] 1e-04  9.734176
##   [2,] 1e-04  9.639650
##   [3,] 1e-04  9.580337
##   [4,] 1e-04  9.613342
##   [5,] 1e-04  9.705725
##   [6,] 1e-04  9.842746
##   [7,] 1e-04  9.931957
##   [8,] 1e-04 10.041652
##   [9,] 1e-04 10.144672
##  [10,] 1e-04 10.437184



# (4)
fftR.propagate.x(vars, nx, fy, nkl, nkr)
