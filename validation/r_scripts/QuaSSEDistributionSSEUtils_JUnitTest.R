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
# (8) testIntegrateOneBranchHiResOutsideClass

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

  # print(paste("numerator", (z * r * r)))
  # print(paste("denominator", (z * lambda - mu + (1-z)*lambda*e0)))
  # print(paste("numerator/denominator", dd))

  # print(paste("scratch", dd))
  # print(paste("d0", d0))
  vars[i,-1] <- dd * d0 # D's are set in place (note that E's are not!!! can only capture E if we throw the return of this function into a variable)
  # print(paste("dd * d0", vars[i,-1]))
  vars
}

normalise <- function(x) x / sum(x)

fftR.make.kern <- function(mean, sd, nx, dx, nkl, nkr) {
  kern <- rep(0, nx)
  xkern <- (-nkl:nkr)*dx
  ikern <- c((nx - nkl + 1):nx, 1:(nkr + 1))
  print("xkern=")
  print(xkern)
  print(paste0("mean=", mean, " sd=", sd))
  print("dnorm=")
  print(dnorm(xkern, mean, sd))
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
pars <- c(.1, .2, 0, 2.5, .03, 0, .01) # 6th and 7th elements are drift and difussion
sd <- 1/20
control.C.1 <- list(dt.max=1/20) # dt = 1/20

cache <- make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)

## these four lines are to make this cache match the previous unit tests
cache$args$drift <- 6
cache$args$diffusion <- 7
cache$control$dx <- 0.01 ## to match unit test
cache$control$xmid <- 0.0 ## to match unit test

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

sp1.y.ds <- sp1.y[nx+1:(length(sp1.y)-nx)]
sp2.y.ds <- sp2.y[nx+1:(length(sp2.y)-nx)]
sp3.y.ds <- sp3.y[nx+1:(length(sp3.y)-nx)]

sp1.y.ds.exp <- sp1.y.ds[1886:1895]
sp2.y.ds.exp <- sp2.y.ds[1886:1895]
sp3.y.ds.exp <- sp3.y.ds[1886:1895]

paste(sp1.y.ds.exp, collapse=", ")
paste(sp2.y.ds.exp, collapse=", ")
paste(sp3.y.ds.exp, collapse=", ")
# ans <- all.branches(pars2, NULL)



# (8) testIntegrateOneBranchHiResOutsideClassJustT

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

# single branch with length 1/20 = 0.05

control.fft$nx <- 1024 * 4 # TODO: try 1032 (multiple of 12... might work)
vars.fft <- matrix(0, control.fft$nx, 2) # high resolution

# vars.fft
# [,1]: E, [,2]: D
# D comes from a normal distn, with states=dnorm(ext.fft$x[[2]], 0, sd), where sd=1/20
vars.fft[seq_len(ext.fft$ndat[1]),2] <- dnorm(ext.fft$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0)

# vars.fft # to see initial E's and D's
# pars.fft$hi$lambda
# pars.fft$hi$mu

# vars.fft.just.t <- fftR.propagate.t(vars.fft, pars.fft$hi$lambda, pars.fft$hi$mu, control.fft$dt.max, ext.fft$ndat[1]) # gotta capture the result of fftR.propagate.t to get the E's (only D's are set in place by this functio
# vars.fft.just.t[,1] # E's
# vars.fft.just.t[2711:2720,-1] # D's

pde.fftR.just.t <- with(control.fft, make.pde.quasse.fftR.2(nx, dx, dt.max, 2L))
ans.fftR.just.t <- pde.fftR.just.t(vars.fft, len, pars.fft$hi, 0) # calculates answer with R; t0 = 0; E's work, but D's are further normalized so then they don't match with the result of fftR.propagate.t

paste(ans.fftR.just.t[[2]][2711:2720,1], collapse=", ") # E's 0.00149145856502394, 0.00149145829907251, 0.00149145803473995, 0.00149145777201677, 0.00149145751089328, 0.00149145725136005, 0.0014914569934076, 0.00149145673702653, 0.00149145648220743, 0.00149145622894108
paste(ans.fftR.just.t[[2]][2711:2720,2], collapse=", ") # D's 2.92247877978117e-274, 4.93457667574128e-275, 8.31118025653625e-276, 1.39633545736575e-276, 2.34008210599557e-277, 3.91189050109308e-278, 6.52313773491418e-279, 1.0850273169382e-279, 1.80027587020561e-280, 2.97955710585536e-281

## pde.fftC <- with(control.fft, make.pde.quasse.fftC(nx, dx, dt.max, 2L, flags)) # partial differential equation
## ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0) # calculates answer with C



# (9) testIntegrateOneBranchHiResOutsideClassJustX

# TODO: Replicating test in testPropagateChOneChQuaSSETest, but building control variable
# TODO: need to change test in Java to use dx = 0.01 instead of 0.001, also update the R code here to make sure it matches

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

pars.fft.just.x <- expand.pars.quasse(lambda, mu, args, ext.fft.just.x, pars.just.x) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x <- matrix(0, control.fft.just.x$nx, 2) # low resolution
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[2]),2] <- dnorm(ext.fft.just.x$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.x[,2], collapse=", ") # 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.047681764029297, 0.0886369682387602, 0.1583090316596, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.99454931271489, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924628, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0.0122380386022754, 0.0058389385158292, 0.0026766045152977, 0.0011788613551308, 0, 0, 0, 0, 0

# checking kernel matches (see PropagatesQuaSSETest -> testMakeNormalKernInPlaceAndFFtAndIfft)
kern.just.x <- fftR.make.kern(-control.fft.just.x$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.just.x$dt * pars.fft.just.x$hi$diffusion), control.fft.just.x$nx, control.fft.just.x$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2]) ## in different orientation than java code, but below the FFT matches
paste(kern.just.x, collapse=", ") # 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 5.09772422059472e-55, 6.62490770689586e-79, 3.90875553004076e-107, 1.04701370374391e-139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.04701370374391e-139, 3.90875553004076e-107, 6.62490770689586e-79, 5.09772422059472e-55, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386

# checking fft-ed fY matches
fy.just.x <- fft(kern.just.x)
paste(Re(fy.just.x), collapse=", ") # 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x
Re(cstep.2[,2])[1:10] # up to here, matches Java
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # this one bombs, have no clue why

# BELOW doesn't matter


control.fft.just.x$nx <- 1032 * 4 # TODO: try 1032 (multiple of 12... might work)
vars.fft.just.x <- matrix(0, control.fft.just.x$nx, 2) # high resolution

## args <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
## pars <- c(.1, .2, 0, 2.5, .03, 0, .01) # 6th and 7th elements are drift and difussion

## ext.fft.just.x <- quasse.extent(control.fft.just.x, drift, diffusion) # prepares X axis stuff

## pars.fft.just.x <- expand.pars.quasse(lambda, mu, args, ext.fft.just.x, pars) # adds lambda and mu vectors to ext.fft

quasse.integrate.fftR.3 <- function (vars, lambda, mu, drift, diffusion, nstep, dt, nx,
    ndat, dx, nkl, nkr) {
    kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
        nx, dx, nkl, nkr)
    fy <- fft(kern)
    for (i in seq_len(nstep)) {
        # vars <- fftR.propagate.t(vars, lambda, mu, dt, ndat) # ignoring propagate t
        vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
    }

    vars
}

make.pde.quasse.fftR.3 <- function (nx, dx, dt.max, nd) {
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
        ans <- quasse.integrate.fftR.3(y, pars$lambda, pars$mu,
            pars$drift, pars$diffusion, nt, dt, nx, ndat, dx,
            padding[1], padding[2])
        q <- sum(ans[, 2]) * dx

        # print(ans[,1]) # E's
        # print(ans[,2]) ## FKM: D's that are returned from fftR.propagate.t

        # ans[,2] <- ans[,2]/q ## FKM: D's are normalized before returning

        list(log(q), ans)
    }
}

# vars.fft
# [,1]: E, [,2]: D
# D comes from a normal distn, with states=dnorm(ext.fft$x[[2]], 0, sd), where sd=1/20
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[1]),2] <- dnorm(ext.fft.just.x$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0)

# checking initial D's match
paste(vars.fft.just.x[1245:1254,2], collapse=", ")
# paste(vars.fft.just.x[1229:1238,2], collapse=", ") # 1.58101006669199e-322, 1.10176639022598e-321, 7.4109846876187e-321, 5.07356011714376e-320, 3.45643385174078e-319, 2.34868926720012e-318, 1.59204205849798e-317, 1.07646593078779e-316, 7.26038680888007e-316, 4.88465194356752e-315
paste(vars.fft.just.x[2778:2787,2], collapse=", ")
# paste(vars.fft[2762:2771,2], collapse=", ") # 4.88465194356752e-315, 7.26038680888007e-316, 1.07646593078779e-316, 1.59204205849798e-317, 2.34868926720012e-318, 3.45643385174078e-319, 5.07356011714376e-320, 7.4109846876187e-321, 1.10176639022598e-321, 1.58101006669199e-322

# checking kernel matches
kern.just.x <- fftR.make.kern(-control.fft.just.x$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.just.x$dt * pars.fft.just.x$hi$diffusion), control.fft.just.x$nx, control.fft.just.x$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2]) ## in different orientation than java code, but below the FFT matches

kern.just.x[1:49]


# checking fft-ed fY matches
fy.just.x <- fft(kern.just.x) ## matches java code
paste(Re(fy.just.x)[1:50], collapse=", ") # 1, 0.999994208125932, 0.999976832705003, 0.999947874341024, 0.999907334040305, 0.999855213211602, 0.999791513666033, 0.999716237616971, 0.999629387679916, 0.999530966872351, 0.999420978613554, 0.999299426724412, 0.999166315427195, 0.999021649345309, 0.998865433503034, 0.998697673325228, 0.998518374637017, 0.998327543663455, 0.998125187029166, 0.997911311757958, 0.99768592527242, 0.99744903539349, 0.997200650340003, 0.99694077872822, 0.996669429571321, 0.996386612278893, 0.996092336656382, 0.995786612904526, 0.995469451618766, 0.995140863788637, 0.994800860797127, 0.994449454420028, 0.994086656825248, 0.993712480572116, 0.993326938610653, 0.992930044280825, 0.992521811311775, 0.992102253821033, 0.991671386313699, 0.991229223681611, 0.990775781202482, 0.990311074539025, 0.989835119738051, 0.989347933229541, 0.988849531825706, 0.988339932720017, 0.987819153486219, 0.987287212077317, 0.986744126824548, 0.986189916436328

paste(Re(fy.just.x)[4077:4126], collapse=", ") # 0.984460727178075, 0.985048196966617, 0.985624599997174, 0.986189916436327, 0.986744126824547, 0.987287212077316, 0.987819153486218, 0.988339932720017, 0.988849531825705, 0.98934793322954, 0.98983511973805, 0.990311074539025, 0.990775781202481, 0.99122922368161, 0.991671386313699, 0.992102253821032, 0.992521811311774, 0.992930044280824, 0.993326938610652, 0.993712480572116, 0.994086656825248, 0.994449454420027, 0.994800860797126, 0.995140863788636, 0.995469451618766, 0.995786612904525, 0.996092336656382, 0.996386612278893, 0.99666942957132, 0.996940778728219, 0.997200650340003, 0.997449035393489, 0.997685925272419, 0.997911311757958, 0.998125187029166, 0.998327543663455, 0.998518374637017, 0.998697673325228, 0.998865433503033, 0.999021649345309, 0.999166315427194, 0.999299426724412, 0.999420978613553, 0.99953096687235, 0.999629387679916, 0.99971623761697, 0.999791513666033, 0.999855213211602, 0.999907334040304, 0.999947874341023

## paste(Re(fy.just.x)[1:50], collapse=", ") # 1, 0.999994117274659, 0.999976469306275, 0.999947056717749, 0.999905880547209, 0.999852942247951, 0.999788243688349, 0.999711787151749, 0.999623575336331, 0.999523611354957, 0.999411898734979, 0.999288441418039, 0.999153243759831, 0.999006310529851, 0.998847646911114, 0.998677258499846, 0.99849515130516, 0.998301331748702, 0.998095806664269, 0.997878583297413, 0.997649669305013, 0.997409072754826, 0.99715680212501, 0.99689286630363, 0.996617274588134, 0.996330036684809, 0.996031162708207, 0.995720663180558, 0.995398549031146, 0.995064831595672, 0.994719522615589, 0.994362634237409, 0.993994179011997, 0.993614169893832, 0.993222620240249, 0.992819543810655, 0.992404954765725, 0.991978867666571, 0.991541297473892, 0.991092259547098, 0.990631769643406, 0.990159843916929, 0.989676498917723, 0.989181751590822, 0.988675619275247, 0.988158119702999, 0.987629270998015, 0.987089091675116, 0.986537600638924, 0.985974817182761

## paste(Re(fy.just.x)[4057:4096], collapse=", ") # 0.990631769643406, 0.991092259547097, 0.991541297473893, 0.991978867666572, 0.992404954765725, 0.992819543810655, 0.993222620240249, 0.993614169893832, 0.993994179011997, 0.994362634237409, 0.994719522615589, 0.995064831595673, 0.995398549031146, 0.995720663180558, 0.996031162708207, 0.996330036684809, 0.996617274588134, 0.99689286630363, 0.99715680212501, 0.997409072754826, 0.997649669305013, 0.997878583297413, 0.998095806664269, 0.998301331748702, 0.99849515130516, 0.998677258499846, 0.998847646911114, 0.999006310529851, 0.99915324375983, 0.999288441418038, 0.999411898734979, 0.999523611354957, 0.999623575336331, 0.999711787151749, 0.999788243688349, 0.999852942247951, 0.999905880547209, 0.999947056717749, 0.999976469306275, 0.999994117274659

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x, 2, fft) ## TODO: continue debugging here
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x
Re(cstep.2[,2])[1:10] # up to here, matches Java
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # this one bombs, have no clue why

Re(apply(apply(vars.fft.just.x, 2, fft) * fy.just.x, 2, ifft))

/ control.fft.just.x$nx

fftR.propagate.x(vars.fft.just.x, control.fft.just.x$nx, fy.just.x, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2])

vars.fft.just.x <- fftR.propagate.x(vars.fft.just.x, control.fft.just.x$nx, fy.just.x, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2]) # gotta capture the result of fftR.propagate.t to get the E's (only D's are set in place by this function
vars.fft.just.x[,1] # E's
vars.fft.just.x[,-1] # D's

# pde.fftR.just.x <- with(control.fft, make.pde.quasse.fftR.3(nx, dx, dt.max, 2L))
# ans.fftR.just.x <- pde.fftR.just.x(vars.fft, len, pars.fft$hi, 0) # calculates answer with R; t0 = 0; E's work, but D's are further normalized so then they don't match with the result of fftR.propagate.t

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
