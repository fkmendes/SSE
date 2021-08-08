# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# PropagatesQuaSSE
# (1) testPropagateTimeOneChQuaSSETest
# (2) testMakeNormalKernInPlaceAndFFT
# (3) testConvolve
# (4) testPropagateChOneChQuaSSETest
# (5) testProcessBranch
# (6) testLogistic
# (7) testDimensions

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
  ndat <- nx - (nkl + 1 + nkr) # this plus 1 here is confusing, need to ask Xia
  i.prev.l <- 1:nkl
  i.prev.r <- (ndat-nkr+1):ndat
  i.zero <- (ndat+1):nx
  print(c(i.prev.l, i.prev.r))
  vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
  vars.out[i.zero,] <- 0
  vars.out[vars.out < 0] <- 0
  vars.out
}

quasse.extent <- function (control, drift, diffusion)
{
    nx <- control$nx
    dx <- control$dx
    dt <- control$dt.max
    xmid <- control$xmid
    r <- control$r
    w <- control$w
    if (control$method == "mol") {
        ndat <- nx * c(r, 1)
        padding <- NULL
    }
    else {
        mean <- drift * dt
        sd <- sqrt(diffusion * dt)
        nkl <- max(ceiling(-(mean - w * sd)/dx)) * c(r, 1)
        nkr <- max(ceiling((mean + w * sd)/dx)) * c(r, 1)
        ndat <- nx * c(r, 1) - (nkl + 1 + nkr)
        padding <- cbind(nkl, nkr)
        storage.mode(padding) <- "integer"
    }
    x0.2 <- xmid - dx * ceiling((ndat[2] - 1)/2)
    x0.1 <- x0.2 - dx * (1 - 1/r)
    x <- list(seq(x0.1, length.out = ndat[1], by = dx/r), seq(x0.2,
        length.out = ndat[2], by = dx))
    tr <- seq(r, length.out = ndat[2], by = r)
    list(x = x, padding = padding, ndat = ndat, tr = tr, nx = c(nx *
        r, nx))
}

expand.pars.quasse <- function (lambda, mu, args, ext, pars)
{
    pars.use <- vector("list", 2)
    for (i in c(1, 2)) {
        x <- list()
        pars.use[[i]] <- list(x = ext$x[[i]], lambda = do.call(lambda,
            c(ext$x[i], pars[args$lambda])), mu = do.call(mu,
            c(ext$x[i], pars[args$mu])), drift = pars[args$drift],
            diffusion = pars[args$diffusion], padding = ext$padding[i,
                ], ndat = ext$ndat[i], nx = ext$nx[i])
    }
    names(pars.use) <- c("hi", "lo")
    pars.use$tr <- ext$tr
    pars.use
}

sigmoid.x <- function (x, y0, y1, xmid, r) {
    to.add = (y1 - y0) / (1 + exp(r * (xmid - x)))
    y0 + to.add
}

## Imports that are generally hidden.
make.pde.quasse.fftC <- diversitree:::make.pde.quasse.fftC
make.pde.quasse.fftR <- diversitree:::make.pde.quasse.fftR
make.pde.quasse.mol <- diversitree:::make.pde.quasse.mol

make.branches.quasse.fftC <- diversitree:::make.branches.quasse.fftC
make.branches.quasse.fftR <- diversitree:::make.branches.quasse.fftR
make.branches.quasse.mol <- diversitree:::make.branches.quasse.mol

make.cache.quasse <- diversitree:::make.cache.quasse



### (1) testPropagateTimeOneChQuaSSETest

vars <- cbind(rep(0.0001,48),
              c(8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688))
# vars <- data.frame(0.0, 1.0) # E and D
lambda <- 1.0
mu <- 0.5
dt <- 0.01
ndat <- 40 # just one useful quant ch bin

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



### (3) testConvolve

vars <- cbind(rep(0.0001,48),
              c(8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264, 11.1384877145132, 9.40550089155345, 9.97880581995152, 11.4504630996011, 10.3369599590198, 10.3149165707367, 10.3046840297378, 8.32290274946024, 9.46368095367558, 8.81487516662079, 9.83971439912364, 11.886850507066, 11.6196319895886, 10.7171936473579, 9.00153746682918, 9.44772548737688))

vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft)) / nx
vars.out

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



# (4) testPropagateChOneChQuaSSETest
res <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
res

##           E         D
##  [1,] 1e-04  8.509768
##  [2,] 1e-04 10.103793
##  [3,] 1e-04  9.089761
##  [4,] 1e-04 11.847721
##  [5,] 1e-04  9.705725
##  [6,] 1e-04  9.842746
##  [7,] 1e-04  9.931957
##  [8,] 1e-04 10.041652
##  [9,] 1e-04 10.144672
## [10,] 1e-04 10.437184
## [11,] 1e-04 10.480127
## [12,] 1e-04 10.377378
## [13,] 1e-04 10.365478
## [14,] 1e-04 10.170118
## [15,] 1e-04 10.001616
## [16,] 1e-04  9.810126
## [17,] 1e-04  9.679028
## [18,] 1e-04  9.558105
## [19,] 1e-04  9.526896
## [20,] 1e-04  9.559773
## [21,] 1e-04  9.767578
## [22,] 1e-04  9.794544
## [23,] 1e-04  9.922514
## [24,] 1e-04 10.012601
## [25,] 1e-04 10.017203
## [26,] 1e-04 10.106680
## [27,] 1e-04  9.906855
## [28,] 1e-04  9.678018
## [29,] 1e-04  9.663567
## [30,] 1e-04  9.480333
## [31,] 1e-04  9.404464
## [32,] 1e-04  9.516894
## [33,] 1e-04  9.679203
## [34,] 1e-04  9.919722
## [35,] 1e-04 10.024894
## [36,] 1e-04 11.450463
## [37,] 1e-04 10.336960
## [38,] 1e-04 10.314917
## [39,] 1e-04 10.304684
## [40,] 0e+00  0.000000
## [41,] 0e+00  0.000000
## [42,] 0e+00  0.000000
## [43,] 0e+00  0.000000
## [44,] 0e+00  0.000000
## [45,] 0e+00  0.000000
## [46,] 0e+00  0.000000
## [47,] 0e+00  0.000000
## [48,] 0e+00  0.000000



# (5) testProcessBranchNoLowRes

## Basic control list.
control.fft <- list(tc=1.3, # time point at which we go from high -> low resolution of X
                    dt.max=1/20, # dt
                    nx=1024, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=5, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.01
sd <- 1/20
len <- 1 # Integrate down a branch length of 1, doesn't get to low res (@tc=1.3)

args <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft <- quasse.extent(control.fft, drift, diffusion) # prepares X axis stuff

# ext.fft
# x[[1]]: high-res X bins
# x[[2]]: low-res X bins
# padding: [1,] left and right padding for high-res
#          [2,] left and right padding for low-res
# ndat: number of X bins actually updated by propagate X
# tr: indices of high-res bins that become each i-th low-res bin
# nx: total number of bins in high-res and low-res, respectively

ndat <- ext.fft$ndat[2] # 999

# control.mol <- modifyList(control.fft, list(nx=ndat, method="mol")) # he didn't finish implementing mol stuff
# ext.mol <- quasse.extent(control.mol, drift, diffusion)

vars.fft <- matrix(0, control.fft$nx, 2)

# vars.fft
# [,1]: E, [,2]: D
# D comes from a normal distn, with states=dnorm(ext.fft$x[[2]], 0, sd), where sd=1/20
vars.fft[seq_len(ndat),2] <- dnorm(ext.fft$x[[2]], 0, sd) # states are the qu character values observed at the tips

pars.fft <- expand.pars.quasse(lambda, mu, args, ext.fft, pars) # adds lambda and mu vectors to ext.fft

# pars.fft
# $hi: high-res stuff
# $hi$x: high-res X bins
# $hi$lambda: high-res lambda values, calculated from the 'lambda' variable (sigmoid.x)
# $hi$mu: high-res mu values, calculated from the 'mu' variable (constant.x)
# $hi$drift, $hi$diffusion, $hi$padding, $hi$ndat, $hi$nx are self-explanatory
# $lo: low-res stuff (counterparts of high-res)

pde.fftC <- with(control.fft, make.pde.quasse.fftC(nx, dx, dt.max, 2L, flags)) # partial differential equation
pde.fftR <- with(control.fft, make.pde.quasse.fftR(nx, dx, dt.max, 2L))
ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0) # calculates answer with C
ans.fftR <- pde.fftR(vars.fft, len, pars.fft$lo, 0) # calculates answer with R



# (6) testLogistic

## See test (5) to get ext.fft done first

## The following code is inside quasse.extent
ndat.lo <- ext.fft$ndat[2]
ndat.hi <- ext.fft$ndat[1]
dx <- 0.01
y0 <- 0.1; y1 <- 0.2; xmid <- 0; r <- 2.5
hi.lo.ratio <- 4
xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

ls.hi <- sigmoid.x(x.hi, y0, y1, xmid, r) # same as pars.fft$hi$lambda
ls.lo <- sigmoid.x(x.lo, y0, y1, xmid, r) # same as pars.fft$hi$lambda

paste(ls.hi[1:10], collapse=", ")
# "0.100000375000363, 0.100000377351446, 0.100000379717269, 0.100000382097925, 0.100000384493506, 0.100000386904106, 0.10000038932982, 0.100000391770742, 0.100000394226967, 0.100000396698591"

paste(ls.hi[2501:2510], collapse=", ")
# "0.195816352937447, 0.195841335181637, 0.195866174682834, 0.195890872179954, 0.195915428409005, 0.195939844103087, 0.19596411999239, 0.195988256804195, 0.196012255262877, 0.196036116089903"

paste(ls.hi[3989:3999], collapse=", ")
# "0.199999600814288, 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637"



# (7) testDimensions

ext.fft$nx # nXbinsHi = 4096, nXbinslo = 1024
ndat.lo # nUsefulXbinsLo = 999
ndat.hi # nUsefulXbinsHi = 3999
ext.fft$padding # nLeftFlanksLo = nRightFlanksLo = 12, nLeftFlanksHi = nRightFlanksHi = 48
xmin.lo # xMinLo = -4.99
xmin.hi # xMinHi = -4.9975

paste(ext.fft$x[[2]][1:10], collapse=", ") # expectedXLoFirst10 = -4.99, -4.98, -4.97, -4.96, -4.95, -4.94, -4.93, -4.92, -4.91, -4.9
paste(ext.fft$x[[2]][990:999], collapse=", ") # expectedXLoLast10 = 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99
paste(ext.fft$x[[1]][1:10], collapse=", ") # expectedXHiFirst10 = -4.9975, -4.995, -4.9925, -4.99, -4.9875, -4.985, -4.9825, -4.98, -4.9775, -4.975
paste(ext.fft$x[[1]][3990:3999], collapse=", ") # expectedXHiLast10 = 4.975, 4.9775, 4.98, 4.9825, 4.985, 4.9875, 4.99, 4.9925, 4.995, 4.9975

paste(pars.fft$lo$lambda[1:10], collapse=", ") # expectedLambdaLoFirt10 = 0.100000382097925, 0.100000391770742, 0.100000401688425, 0.100000411857174, 0.100000422283344, 0.100000432973452, 0.100000443934178, 0.100000455172374, 0.100000466695064, 0.10000047850945
paste(pars.fft$lo$lambda[989:999], collapse=", ") # expectedLambdaLoLast10 = 0.199999509377086, 0.199999521490551, 0.199999533304936, 0.199999544827626, 0.199999556065822, 0.199999567026548, 0.199999577716656, 0.199999588142826, 0.199999598311575, 0.199999608229258, 0.199999617902075
paste(pars.fft$hi$lambda[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.100000375000363, 0.100000377351446, 0.100000379717269, 0.100000382097925, 0.100000384493506, 0.100000386904106, 0.10000038932982, 0.100000391770742, 0.100000394226967, 0.100000396698591
paste(pars.fft$hi$lambda[3989:3999], collapse=", ") # expectedLambdaHiLast10 = 0.199999600814288, 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637

paste(pars.fft$lo$mu[1:10], collapse=", ") # expectedMuLoFirt10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$lo$mu[989:999], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$hi$mu[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$hi$mu[3989:3999], collapse=", ") # expectedLambdaHiLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03



# (8)

if (getRversion() >= "3.6.0") {
    RNGkind(sample.kind = "Rounding")
}

lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(drift, diffusion)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=3, x0=0, single.lineage=FALSE, verbose=FALSE)
pars <- c(.1, .2, 0, 2.5, .03, 0, .01)
sd <- 1/20
control.C.1 <- list(dt.max=1/20) # dt = 1/20

cache <- make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)
cache$control

make.all.branches.quasse <- function(cache, control) {
  branches <- diversitree:::make.branches.quasse(cache, control)
  initial.conditions <- make.initial.conditions.quasse(control)
  ## TODO: This is where tips.combined goes, *not* in the likelihood
  ## function...

  # FKM: this function does not seem to execute
  function(pars, intermediates, preset=NULL) {
      cache$y = initial.tip.quasse(cache, cache$control, pars[[1]]$x)

      # FKM
      print(cache$y)

      diversitree:::all.branches.list(pars, cache, initial.conditions,
                                      branches, preset)
  }
}

# inside R/model-quasse.R
initial.tip.quasse <- function(cache, control, x) {
  nx = control$nx * control$r
  npad = nx - length(x)
  e0 = 1 - cache$sampling.f

  # what is tips.combined?
  if ( control$tips.combined ) {
    tips = cache$tips
    t = cache$len[tips]
    i = order(t)
    target = tips[i]

    states = cache$states[i]
    states.sd = cache$states.sd[i]

    y = mapply(function(mean, sd)
                c(dnorm(x, mean, sd), rep(0, npad)),
                states, states.sd, SIMPLIFY=FALSE)
    y = matrix(c(rep(e0, nx), unlist(y)), nx, length(target)+1)

    list(target=target, y=y, t=t[i])
  } else {
    # FKM: seems to be the default
      y = mapply(function(mean, sd)
          c(rep(e0, nx), dnorm(x, mean, sd), rep(0, npad)),
          cache$states, cache$states.sd,
          SIMPLIFY=FALSE) # returns list
    diversitree:::dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
  }
}

## taking following lines from R/model-quasse.R, make.quasse()
# all.branches <- make.all.branches.quasse(cache, cache$control) # all.branches is a function
f.pars <- diversitree:::make.pars.quasse(cache)
pars2 <- f.pars(pars)

# understanding tip initialization
sampling.f <- 1
e0 <- 1 - sampling.f
nx <- cache$control$nx * cache$control$r # 1024 * 4 = 4096
npad <- nx - length(pars2[[1]]$x)

# this is how the tip y's are initialized
sp1.y <- c(rep(e0, nx), dnorm(pars2[[1]]$x, cache$states[1], cache$states.sd[1]), rep(0, npad)) # E's and D's (8192 elements)
sp2.y <- c(rep(e0, nx), dnorm(pars2[[1]]$x, cache$states[2], cache$states.sd[2]), rep(0, npad))
sp3.y <- c(rep(e0, nx), dnorm(pars2[[1]]$x, cache$states[3], cache$states.sd[3]), rep(0, npad))

# we can compare it to running the function as is
cache$y <- initial.tip.quasse(cache, cache$control, pars2[[1]]$x) # this is a list
# cache$y$t # branch lengths
# cache$y$target # tip indices (I think)
identical(cache$y$y$sp1, sp1.y)
identical(cache$y$y$sp2, sp2.y)
identical(cache$y$y$sp3, sp3.y)
identical(cache$y$y$sp1, sp2.y) # just as a control

# ans <- all.branches(pars2, NULL)



###################

# inside R/model-quasse.R
make.initial.conditions.quasse <- function(control) {
  tc = control$tc
  r = control$r
  nx.lo = control$nx
  nx.hi = nx.lo * r

  ## There is the chance that we could be slightly off on the depth
  ## by rounding error.  Because of this, I've done the testing
  ## against the *length* of the data, and then checked that the time
  ## is appropriate (to within eps of the correct value).  It is
  ## possible that two different branches with different numbers of
  ## nodes that finish right at the critical interval might have
  ## conflicting lengths.
  eps = 1e-8

  function(init, pars, t, idx) {
    if ( length(init[[1]]) != length(init[[2]]) )
      stop("Data have incompatible length")

    if ( t < tc ) {
      nx <- nx.hi
      lambda <- pars[[1]]$lambda
    } else {
      nx <- nx.lo
      lambda <- pars[[2]]$lambda
    }

    ndat = length(lambda)
    i = seq_len(nx)
    j = seq.int(nx+1, nx + ndat)

    c(init[[1]][i],
      init[[1]][j] * init[[2]][j] * lambda,
      rep.int(0.0, nx - ndat))
  }
}
