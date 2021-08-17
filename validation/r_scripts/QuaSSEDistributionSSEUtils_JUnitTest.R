# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# PropagatesQuaSSE
# (1) testPropagateTimeOneChQuaSSETest
# (2) testMakeNormalKernInPlaceAndFFT
# (3) testConvolve
# (4) testPropagateChOneChQuaSSETest
# (5) testLogistic
# (6) testDimensions
# (7) testInitializationOfTips
# (8) testIntegrateOneBranchLoRes48BinsOutsideClassJustT
# (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX
# (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX and testPropagateChOneCh1024QuaSSETest

library(diversitree)

## First, diversitree's functions

fftR.propagate.t <- function(vars, lambda, mu, dt, ndat) {
  i <- seq_len(ndat) # creates a seq that starts with 1 and ends in ndat, with increments of 1
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  print(paste("mu =", mu, "lambda =", lambda, "e0 =", e0, "z =", z))
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) / (mu + z*(e0 - 1)*lambda - lambda*e0) # E's
  print(paste("e =", vars[i,1]))
  dd <- (z * r * r)/(z * lambda - mu + (1-z)*lambda*e0)^2

  # for debugging
  ## print(paste("numerator", (z * r * r)))
  ## print(paste("denominator", (z * lambda - mu + (1-z)*lambda*e0)))
  ## print(paste("numerator/denominator", dd))

  ## print(paste("scratch", dd))
  ## print(paste("d0", d0))
  vars[i,-1] <- dd * d0 # D's are set in place (note that E's are not!!! can only capture E if we throw the return of this function into a variable)
  ## print(paste("dd * d0", vars[i,-1]))
  vars
}

normalise <- function(x) x / sum(x)

fftR.make.kern <- function(mean, sd, nx, dx, nkl, nkr) {
  kern <- rep(0, nx)
  xkern <- (-nkl:nkr)*dx
  ikern <- c((nx - nkl + 1):nx, 1:(nkr + 1))
  ## print("xkern=")
  ## print(xkern)
  ## print(paste0("mean=", mean, " sd=", sd))
  ## print("dnorm=")
  ## print(dnorm(xkern, mean, sd))
  kern[ikern] <- normalise(dnorm(xkern, mean, sd))
  kern
}

ifft <- function(x) fft(x, inverse=TRUE)

fftR.propagate.x <- function(vars, nx, fy, nkl, nkr) {
  print(paste("nkl=", nkl, " nkr=", nkl, " nx=", nx))
  vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft))/nx
  ndat <- nx - (nkl + 1 + nkr) # this plus 1 here is confusing, need to ask Xia
  i.prev.l <- 1:nkl
  i.prev.r <- (ndat-nkr+1):ndat
  i.zero <- (ndat+1):nx

  # print(c(i.prev.l, i.prev.r))

  vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
  vars.out[i.zero,] <- 0
  vars.out[vars.out < 0] <- 0
  vars.out
}

quasse.extent <- function (control, drift, diffusion) {
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
        print(paste("mean =", mean, "sd =", sd, "w =", w, "r =", r))
        print(paste("nkl =", nkl, "nkr =", nkr))
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

expand.pars.quasse <- function (lambda, mu, args, ext, pars) {
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

# changed from 0.001 -> 0.01, updating!
dx <- 0.01 # (dx * nx) gives the range of x values we're considering

# '-' in -dt because we're going back in time
kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion), nx, dx, nkl, nkr)

paste(kern, collapse=", ") # 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386

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
## nkl = nkr = 4
##
## x is at 10th bin, so we ignore anything beyond nkl and nkr bins away from 10th, meaning
## x will influence from 5th - 15th bins
##
## the first bin gives the probability that x stays where it is
## the 2nd - 6th bin give the probability x changes to 11th - 15th bins
## the 95th - 100th bins give the probability x changes to the 9th - 5th

fy <- fft(kern)
ify <- ifft(fy) # note that this is NOT what we do in the likelihood calculation

paste(Re(fy)[1:10], collapse=", ") # 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113

# when dx = 0.001
# [1]  1.000000000  0.956873917  0.835109753  0.655887828  0.449367562
# [5]  0.248270676  0.081152021 -0.033081907 -0.088189035 -0.090559850

paste(Re(ify)[1:10], collapse=", ") # 47.3617577773852, 0.319121013687545, 9.76198587343688e-08, 7.45931094670027e-17, -4.44089209850063e-16, 1.51614216742138e-16, 4.05992386804242e-16, -1.43599463769065e-16, -4.44089209850063e-16, -5.17654502592123e-16

# when dx = 0.001
# [1]  7.149417e+00  6.800735e+00  5.853447e+00  4.558669e+00  3.212440e+00
# [6]  5.957485e-16  4.549779e-16  1.366926e-15  2.766926e-16  3.584917e-17



### (3) testConvolve

# see vars.fft.just.x below to see where this came from
vars <- cbind(rep(0.0, 48), c(0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0))

paste(Re(apply(apply(vars, 2, fft) * fy, 2, ifft))[,2], collapse=", ") # 0.280447809353859, 0.58934289318654, 1.18632299309589, 2.29444263475988, 4.26373864030165, 7.61277219718136, 13.0597170956171, 21.5259924869188, 34.090305996801, 51.8724486751153, 75.8369022635638, 106.527631278782, 143.774440226047, 186.439841249409, 232.29148807448, 278.077162458051, 319.841447837597, 353.46100325742, 375.306060575871, 382.883752104488, 375.306060575871, 353.46100325742, 319.841447837597, 278.077162458051, 232.29148807448, 186.439841249409, 143.774440226047, 106.527631278782, 75.8369022635638, 51.8724486751154, 34.090305996801, 21.5259924869189, 13.0597170956171, 7.61277219718137, 4.26373864030163, 2.29444263475987, 1.18632299309587, 0.589342893186554, 0.280447809353873, 0.00186332917266441, 5.70054226045613e-10, 0, 2.8421709430404e-14, 0, 5.6843418860808e-14, 0, 5.70025804336183e-10, 0.00186332917272125

vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft)) / nx
paste(vars.out[,1], collapse=", ") # 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(vars.out[,2], collapse=", ") # 0.00584266269487207, 0.0122779769413863, 0.0247150623561643, 0.0478008882241642, 0.0888278883396176, 0.158599420774612, 0.272077439492023, 0.448458176810809, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164128, 2.99530083804266, 3.88416335936269, 4.83940600155166, 5.79327421787606, 6.66336349661661, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661661, 5.79327421787606, 4.83940600155166, 3.88416335936269, 2.99530083804266, 2.21932565164129, 1.57993546382425, 1.0806760140649, 0.710214708266689, 0.44845817681081, 0.272077439492023, 0.158599420774612, 0.0888278883396172, 0.0478008882241639, 0.024715062356164, 0.0122779769413865, 0.00584266269487236, 3.88193577638418e-05, 1.18761297092836e-11, 0, 5.9211894646675e-16, 0, 1.1842378929335e-15, 0, 1.18755375903371e-11, 3.88193577650261e-05



# (4) testPropagateChOneChQuaSSETest
res <- fftR.propagate.x(vars, nx, fy, nkl, nkr)

paste(res[,1], collapse=", ") # 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(res[,2], collapse=", ") # 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0247150623561646, 0.0478008882241639, 0.088827888339617, 0.158599420774613, 0.272077439492024, 0.448458176810809, 0.710214708266688, 1.0806760140649, 1.57993546382425, 2.21932565164128, 2.99530083804266, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661661, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661661, 5.79327421787606, 4.83940600155166, 3.8841633593627, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.0806760140649, 0.710214708266687, 0.44845817681081, 0.272077439492023, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0, 0, 0, 0, 0, 0, 0, 0, 0



### PREPARING THINGS FOR A FEW TESTS BELOW ###

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
len <- 1/20 # Integrate down a branch length of 1, doesn't get to low res (@tc=1.3)

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

pars.fft <- expand.pars.quasse(lambda, mu, args, ext.fft, pars) # adds lambda and mu vectors to ext.fft

# pars.fft
# $hi: high-res stuff
# $hi$x: high-res X bins
# $hi$lambda: high-res lambda values, calculated from the 'lambda' variable (sigmoid.x)
# $hi$mu: high-res mu values, calculated from the 'mu' variable (constant.x)
# $hi$drift, $hi$diffusion, $hi$padding, $hi$ndat, $hi$nx are self-explanatory
# $lo: low-res stuff (counterparts of high-res)



# (5) testLogistic

## See above to get ext.fft done first

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



# (6) testDimensions

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



# (7) testInitializationOfTips

if (getRversion() >= "3.6.0") {
    RNGkind(sample.kind = "Rounding")
}

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(drift, diffusion)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=3, x0=0, single.lineage=FALSE, verbose=FALSE)
pars <- c(.1, .2, 0, 2.5, .03, 0, .001) # 6th and 7th elements are drift and diffusion
sd <- 1/20
control.C.1 <- list(dt.max=0.01) # dt = 0.01

cache <- make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)

## these four lines are to make this cache match the previous unit tests
cache$args$drift <- 6
cache$args$diffusion <- 7
cache$control$dx <- 0.01 ## to match unit test
cache$control$xmid <- 0.0 ## to match unit test
cache$control$w <- 10

# inside R/model-quasse.R
initial.tip.quasse <- function(cache, control, x) {

  print(paste("nx=", control$nx))

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

make.pars.quasse <- function (cache) {
    args <- cache$args
    function(pars) {
        names(pars) <- NULL
        drift <- pars[args$drift]
        diffusion <- pars[args$diffusion]
        ext <- quasse.extent(cache$control, drift, diffusion)
        pars <- expand.pars.quasse(cache$lambda, cache$mu, args,
            ext, pars)
        diversitree:::check.pars.quasse(pars$hi$lambda, pars$hi$mu, drift,
            diffusion)
        pars
    }
}

## taking following lines from R/model-quasse.R, make.quasse()
# all.branches <- make.all.branches.quasse(cache, cache$control) # all.branches is a function
f.pars <- make.pars.quasse(cache)
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
identical(cache$y$y$sp1, sp1.y) # TRUE
identical(cache$y$y$sp2, sp2.y) # TRUE
identical(cache$y$y$sp3, sp3.y) # TRUE
identical(cache$y$y$sp1, sp2.y) # FALSE (just as a control)

sp1.y.ds <- sp1.y[nx+1:(length(sp1.y)-nx)]
sp2.y.ds <- sp2.y[nx+1:(length(sp2.y)-nx)]
sp3.y.ds <- sp3.y[nx+1:(length(sp3.y)-nx)]

sp1.y.ds.exp <- sp1.y.ds[1886:1895]
sp2.y.ds.exp <- sp2.y.ds[1886:1895]
sp3.y.ds.exp <- sp3.y.ds[1886:1895]

paste(sp1.y.ds.exp, collapse=", ") # 8.13623712450468e-08, 1.1005578408194e-07, 1.48496566106919e-07, 1.99863833658288e-07, 2.68328176840003e-07, 3.59345830399246e-07, 4.80035331999228e-07, 6.39658332527801e-07, 8.50231484597812e-07, 1.12730275515906e-06
paste(sp2.y.ds.exp, collapse=", ") # 1.82648393304909e-11, 2.63062707587468e-11, 3.77934883411903e-11, 5.41612820961595e-11, 7.74239202284714e-11, 1.10401669802924e-10, 1.57032805818765e-10, 2.2280216312045e-10, 3.15328103714112e-10, 4.45164187528514e-10
paste(sp3.y.ds.exp, collapse=", ") # 3.51799453441839e-17, 5.49394711042998e-17, 8.55831075496442e-17, 1.32985990766802e-16, 2.06128478871848e-16, 3.18701690685798e-16, 4.91524307682648e-16, 7.56170789433044e-16, 1.16040357229278e-15, 1.77628431652675e-15



# (8) testIntegrateOneBranchLoRes48BinsOutsideClassJustT

quasse.integrate.fftR.2 <- function (vars, lambda, mu, drift, diffusion, nstep, dt, nx,
    ndat, dx, nkl, nkr) {
    kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
        nx, dx, nkl, nkr)
    fy <- fft(kern)
    for (i in seq_len(nstep)) {
        vars <- fftR.propagate.t(vars, lambda, mu, dt, ndat)
        # print(paste("nx = ", nx, " length(fy)=", length(fy)))
        # vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr) # ignoring propagate X
    }

    vars
}

make.pde.quasse.fftR.2 <- function (nx, dx, dt.max, nd) {
    function(y, len, pars, t0) {
        padding <- pars$padding
        ndat <- length(pars$lambda)
        nt <- as.integer(ceiling(len/dt.max))
        dt <- len/nt

        if (!(length(y) %in% (nd * nx)))
            stop("Wrong size y")
        if (length(pars$lambda) != length(pars$mu) || length(pars$lambda) >
            (nx - 3))
            stop("Incorrect length pars")
        if (pars$diffusion <= 0)
            stop("Invalid diffusion parameter")
        if (!is.matrix(y))
            y <- matrix(y, nx, nd)
        ans <- quasse.integrate.fftR.2(y, pars$lambda, pars$mu,
            pars$drift, pars$diffusion, nt, dt, nx, ndat, dx,
            padding[1], padding[2])
        q <- sum(ans[, 2]) * dx

        # print(ans[,1]) # E's
        # print(ans[,2]) ## FKM: D's that are returned from fftR.propagate.t

        # ans[,2] <- ans[,2]/q ## FKM: D's are normalized before returning

        list(log(q), ans)
    }
}

control.fft.just.t <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=48, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x

args.fft.just.t <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.t <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.t <- quasse.extent(control.fft.just.t, drift, diffusion) # prepares X axis stuff
ext.fft.just.t$padding

pars.fft.just.t <- expand.pars.quasse(lambda, mu, args.fft.just.t, ext.fft.just.t, pars.just.t) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.t <- matrix(0, control.fft.just.t$nx, 2) # low resolution
vars.fft.just.t[seq_len(ext.fft.just.t$ndat[2]),2] <- dnorm(ext.fft.just.t$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.t[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0

vars.fft.just.t.tmp <- fftR.propagate.t(vars.fft.just.t, pars.fft.just.t$lo$lambda, pars.fft.just.t$lo$mu, control.fft.just.t$dt.max, ext.fft.just.t$ndat[2])
paste(vars.fft.just.t.tmp[,1], collapse=", ") # 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(vars.fft.just.t.tmp[,2], collapse=", ") # 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0

pde.fftR.just.t <- with(control.fft.just.t, make.pde.quasse.fftR.2(nx, dx, dt.max, 2L))
ans.fftR.just.t <- pde.fftR.just.t(vars.fft.just.t, control.fft.just.t$dt.max, pars.fft.just.t$lo, 0)

## same as above
paste(ans.fftR.just.t[[2]][,1], collapse=", ") # 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.just.t[[2]][,2], collapse=", ") #  0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX

control.fft.just.x <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=48, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x

args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x <- quasse.extent(control.fft.just.x, drift, diffusion) # prepares X axis stuff
ext.fft.just.x$padding

pars.fft.just.x <- expand.pars.quasse(lambda, mu, args.fft.just.x, ext.fft.just.x, pars.just.x) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x <- matrix(0, control.fft.just.x$nx, 2) # low resolution
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[2]),2] <- dnorm(ext.fft.just.x$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.x[,2], collapse=", ") # 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.047681764029297, 0.0886369682387602, 0.1583090316596, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.99454931271489, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924628, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0.0122380386022754, 0.0058389385158292, 0.0026766045152977, 0.0011788613551308, 0, 0, 0, 0, 0

# checking kernel matches (see PropagatesQuaSSETest -> testMakeNormalKernInPlaceAndFFtAndIfft)
kern.just.x <- fftR.make.kern(-control.fft.just.x$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.just.x$dt * pars.fft.just.x$hi$diffusion), control.fft.just.x$nx, control.fft.just.x$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2])
paste(kern.just.x, collapse=", ") # 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 5.09772422059472e-55, 6.62490770689586e-79, 3.90875553004076e-107, 1.04701370374391e-139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.04701370374391e-139, 3.90875553004076e-107, 6.62490770689586e-79, 5.09772422059472e-55, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386

# checking fft-ed fY matches
fy.just.x <- fft(kern.just.x)
paste(Re(fy.just.x), collapse=", ") # 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints

ds.prop.x <- fftR.propagate.x(vars.fft.just.x, control.fft.just.x$nx, fy.just.x, pars.fft.just.x$lo$padding[1], pars.fft.just.x$lo$padding[2])
paste(ds.prop.x[,1], collapse=", ") # 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ds.prop.x[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0888278883396178, 0.158599420774613, 0.272077439492024, 0.44845817681081, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804265, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661662, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.08067601406491, 0.710214708266689, 0.448458176810809, 0.272077439492024, 0.158599420774613, 0.088827888339617, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (10) testPropagateChOneCh1024QuaSSETest and testIntegrateOneBranchLoRes1024BinsOutsideClassJustX

control.fft.just.x.1024 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=1024, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

drift <- 0.0
diffusion <- 0.001
lambda <- sigmoid.x
mu <- constant.x

args.fft.just.x.1024 <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x.1024 <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x.1024 <- quasse.extent(control.fft.just.x.1024, drift, diffusion) # prepares X axis stuff
# ext.fft.just.x.1024$padding

pars.fft.just.x.1024 <- expand.pars.quasse(lambda, mu, args.fft.just.x.1024, ext.fft.just.x.1024, pars.just.x.1024) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x.1024 <- matrix(0, control.fft.just.x.1024$nx, 2) # low resolution
vars.fft.just.x.1024[seq_len(ext.fft.just.x.1024$ndat[2]),2] <- dnorm(ext.fft.just.x.1024$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.1024$lo$padding should be 0.0
paste(vars.fft.just.x.1024[316:326,2], collapse=", ") # 5.07356011714376e-320, 1.07646593078779e-316, 2.19444210413816e-313, 4.29809786779161e-310, 8.0882896186969e-307, 1.46239706910102e-303, 2.54040022766109e-300, 4.24001310304921e-297, 6.79924162752444e-294, 1.04756739392715e-290, 1.55071393737006e-287

kern.just.x.1024 <- fftR.make.kern(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$diffusion), control.fft.just.x.1024$nx, control.fft.just.x.1024$dx, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])

fy.just.x.1024 <- fft(kern.just.x.1024)

ds.prop.x.1024 <- fftR.propagate.x(vars.fft.just.x.1024, control.fft.just.x.1024$nx, fy.just.x.1024, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])
paste(ds.prop.x.1024[1:48,2], collapse=", ") # 0, 0, 0, 0, 1.11022302462516e-16, 0, 0, 0, 2.77555756156289e-17, 2.77555756156289e-17, 0, 8.32667268468867e-17, 0, 1.04083408558608e-16, 5.89805981832114e-17, 0, 0, 4.99275100429575e-17, 1.74936020530536e-16, 0, 5.47996435555642e-17, 1.64060117487791e-16, 2.55997708857483e-16, 0, 8.5609964086127e-17, 9.49329645689794e-17, 2.91487168428854e-16, 3.34092533592135e-17, 0, 2.19845474956207e-16, 3.57227590802606e-16, 3.37917798902023e-16, 0, 4.14102761209603e-17, 0, 0, 0, 0, 0, 0, 0, 3.49655200626575e-17, 2.86229373536173e-17, 0, 0, 0, 1.55257751099924e-16, 0
paste(ds.prop.x.1024[977:1024,2], collapse=", ") # 0, 0, 0, 0, 0, 0, 7.52487427515263e-17, 0, 0, 0, 1.55620691069343e-16, 1.14818872876954e-17, 0, 0, 1.58130728038333e-16, 1.0696554404076e-16, 0, 1.57854754116206e-16, 1.40605775178319e-16, 0, 0, 0, 0, 0, 0, 1.04083408558608e-16, 1.31838984174237e-16, 1.17961196366423e-16, 8.32667268468867e-17, 5.55111512312578e-17, 5.55111512312578e-17, 1.11022302462516e-16, 1.11022302462516e-16, 1.11022302462516e-16, 1.11022302462516e-16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (11) testPropagateChOneCh4096QuaSSETest and testIntegrateOneBranchLoRes4096BinsOutsideClassJustX (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)

control.fft.just.x.4096 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=4096, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

drift <- 0.0
diffusion <- 0.001
lambda <- sigmoid.x
mu <- constant.x

args.fft.just.x.4096 <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x.4096 <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x.4096 <- quasse.extent(control.fft.just.x.4096, drift, diffusion) # prepares X axis stuff
# ext.fft.just.x.4096$padding

pars.fft.just.x.4096 <- expand.pars.quasse(lambda, mu, args.fft.just.x.4096, ext.fft.just.x.4096, pars.just.x.4096) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x.4096 <- matrix(0, control.fft.just.x.4096$nx, 2) # low resolution
vars.fft.just.x.4096[seq_len(ext.fft.just.x.4096$ndat[2]),2] <- dnorm(ext.fft.just.x.4096$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.4096$lo$padding should be 0.0
paste(vars.fft.just.x.4096[1852:1861,2], collapse=", ") # 5.07356011714376e-320, 1.07646593078779e-316, 2.19444210413816e-313, 4.29809786778964e-310, 8.08828961870339e-307, 1.46239706909986e-303, 2.54040022766244e-300, 4.2400131030559e-297, 6.79924162752622e-294, 1.04756739392852e-290

kern.just.x.4096 <- fftR.make.kern(-control.fft.just.x.4096$dt * pars.fft.just.x.4096$hi$drift, sqrt(control.fft.just.x.4096$dt * pars.fft.just.x.4096$hi$diffusion), control.fft.just.x.4096$nx, control.fft.just.x.4096$dx, pars.fft.just.x.4096$hi$padding[1], pars.fft.just.x.4096$hi$padding[2])

fy.just.x.4096 <- fft(kern.just.x.4096)

ds.prop.x.4096 <- fftR.propagate.x(vars.fft.just.x.4096, control.fft.just.x.4096$nx, fy.just.x.4096, pars.fft.just.x.4096$lo$padding[1], pars.fft.just.x.4096$lo$padding[2])
paste(ds.prop.x.4096[1:48,2], collapse=", ") # 0, 0, 0, 0, 0, 1.11022302462516e-16, 0, 8.32667268468867e-17, 0, 1.38777878078145e-16, 3.46944695195361e-17, 3.46944695195361e-17, 6.93889390390723e-17, 9.36750677027476e-17, 9.54097911787244e-17, 4.33680868994202e-17, 7.89299181569447e-17, 7.88214979396962e-17, 0, 0, 2.81181057170538e-17, 0, 3.0465657530369e-18, 1.259640695312e-16, 0, 0, 0, 0, 0, 0, 0, 0, 4.0355134457755e-17, 3.12301577048253e-17, 6.55884080838606e-17, 2.91817144010091e-17, 0, 4.28408088897544e-17, 0, 4.61700718889374e-17, 3.56363701568829e-17, 3.91532509538828e-17, 7.02834058313728e-17, 4.43438688546571e-17, 0, 5.26922255827955e-17, 4.46691295064028e-17, 3.03576608295941e-18
paste(ds.prop.x.4096[4049:4096,2], collapse=", ") # 4.01024466275727e-18, 4.03692777680313e-18, 1.05349692261114e-16, 5.02042488979581e-17, 2.39443578005289e-18, 0, 6.6119663477603e-17, 1.00144116158942e-16, 0, 0, 1.93444933680503e-17, 1.37819341184318e-17, 0, 0, 4.30475908080028e-17, 8.53335931173056e-17, 0, 8.7890679705948e-17, 0, 5.47860910284081e-17, 0, 7.38341679462629e-17, 0, 7.58941520739853e-18, 1.17093834628434e-17, 3.46944695195361e-18, 0, 0, 0, 8.32667268468867e-17, 2.77555756156289e-17, 5.55111512312578e-17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


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
