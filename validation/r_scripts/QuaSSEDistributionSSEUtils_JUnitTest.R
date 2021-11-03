# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests

## in PropagatesQuaSSETest
## (1) testPropagateTimeOneChQuaSSETest
## (2) testMakeNormalKernInPlaceAndFFT
## (3) testConvolve
## (4) testPropagateChOneCh48QuaSSETest
## (5) testLogistic
## (9.1) testPropagateChOneCh48QuaSSETest
## (10.1) testPropagateChOneCh1024QuaSSETest
## (11.1) testPropagateChOneCh4096QuaSSETest

## in QuaSSEDistributionTest
## (6) testDimensions
## (7) testInitializationOfTips
## (8) testIntegrateOneBranchLoRes48BinsOutsideClassJustT
## (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX
## (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX
## (11) testIntegrateOneBranchLoRes4096BinsOutsideClassJustX (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)
## (12) testBothXandTPropagateMethodsInsideClassOneBranchLoRes48Bins
## (13) testBothXandTPropagateMethodsInsideClassOneBranchLoRes48BinsLargerDt
## (14) testIntegrateOneBranchHiRes48BinsInsideClassBothXandTTwoDt
## (15) testPriorProbAtRootObserved
## (16) testPriorProbAtRootFlat
## (17) testRootCalcProcedure
## (18) testPruneTree32Bins

## IMPORTANT: expectations (outputs from paste(xyz, collapse=", ") will vary depending on machine architecture; when dnorm() gets very small inputs (e.g., dx * nx are both large), the initial values of D will be tiny, and then the result of FFT on different machines will differ
##
## This is the reason for 48 bins having a dx = 0.01, and 1024/4096 bins having a dx = 0.0005


library(diversitree)

## First, diversitree's functions

fftR.propagate.t <- function(vars, lambda, mu, dt, ndat) {
  i <- seq_len(ndat) # creates a seq that starts with 1 and ends in ndat, with increments of 1
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  ## print(paste("mu =", mu, "lambda =", lambda, "e0 =", e0, "z =", z))
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) / (mu + z*(e0 - 1)*lambda - lambda*e0) # E's
  ## print(paste("e =", vars[i,1]))
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
    print("vars="); print(vars)

    vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft))/nx

    print("vars after fft-ing, multiplying by fy kernel, i-ffting, and normalizing"); print(vars.out)
    ndat <- nx - (nkl + 1 + nkr) # this plus 1 here is confusing, need to ask Xia
    i.prev.l <- 1:nkl
    i.prev.r <- (ndat-nkr+1):ndat
    i.zero <- (ndat+1):nx

    # print(c(i.prev.l, i.prev.r))

    vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
    vars.out[i.zero,] <- 0
    vars.out[vars.out < 0] <- 0

    print("vars after rearranging it"); print(vars.out)

    vars.out
}

quasse.extent.debug <- function (control, drift, diffusion) {
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

        print(paste("nx =", nx))
        print(paste("diffusion = ", diffusion, " drift = ", drift))
        print(paste("mean =", mean, "sd =", sd, "w =", w, "r =", r))
        print(paste("nkl =", nkl, "nkr =", nkr))

        ndat <- nx * c(r, 1) - (nkl + 1 + nkr)

        print(paste("ndat =", ndat))

        padding <- cbind(nkl, nkr)
        storage.mode(padding) <- "integer"
    }
    x0.2 <- xmid - dx * ceiling((ndat[2] - 1)/2)
    x0.1 <- x0.2 - dx * (1 - 1/r)
    print(paste("xMinLo =", x0.2, " xMinHi =", x0.1))
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
        pars.use[[i]] <- list(x = ext$x[[i]],
                              lambda = do.call(lambda, c(ext$x[i], pars[args$lambda])),
                              mu = do.call(mu, c(ext$x[i], pars[args$mu])),
                              drift = pars[args$drift],
                              diffusion = pars[args$diffusion],
                              padding = ext$padding[i,],
                              ndat = ext$ndat[i],
                              nx = ext$nx[i])
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
res[1:10,]

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

## plot(kern)

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
paste(res[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0888278883396176, 0.158599420774612, 0.272077439492023, 0.448458176810809, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164128, 2.99530083804266, 3.88416335936269, 4.83940600155166, 5.79327421787606, 6.66336349661661, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661661, 5.79327421787606, 4.83940600155166, 3.88416335936269, 2.99530083804266, 2.21932565164129, 1.57993546382425, 1.0806760140649, 0.710214708266689, 0.44845817681081, 0.272077439492023, 0.158599420774612, 0.0888278883396172, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (5) testLogistic

## The following code is inside quasse.extent
drift <- 0.0; diffusion <- 0.01
nx <- 1024; dx <- 0.01; dt <- 0.05
y0 <- 0.1; y1 <- 0.2; xmid <- 0; r <- 2.5; death <- 0.03
hi.lo.ratio <- 4 # r=4L when making control.fft
w <- 5

mean4test <- drift * dt
sd4test <- sqrt(diffusion * dt)
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2];
ndat.hi <- ndat[1];

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

## Requires some variables from (5)
# modifying some of them though...
dt <- 0.01; dx <- 0.0005; diffusion <- 0.001; w <- 10

sd4test <- sqrt(diffusion * dt);
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2]; ndat.lo
ndat.hi <- ndat[1]; ndat.hi

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

paste(x.lo[1:10], collapse=", ") # expectedXLoFirst10 = -0.2235, -0.223, -0.2225, -0.222, -0.2215, -0.221, -0.2205, -0.22, -0.2195, -0.219
paste(x.lo[886:895], collapse=", ") # expectedXLoLast10 = 0.219, 0.2195, 0.22, 0.2205, 0.221, 0.2215, 0.222, 0.2225, 0.223, 0.2235

paste(x.hi[1:10], collapse=", ") # expectedXHiFirst10 = -0.223875, -0.22375, -0.223625, -0.2235, -0.223375, -0.22325, -0.223125, -0.223, -0.222875, -0.22275
paste(x.hi[3574:3583], collapse=", ") # expectedXHiFirst10 = 0.22275, 0.222875, 0.223, 0.223125, 0.22325, 0.223375, 0.2235, 0.223625, 0.22375, 0.223875

lambda <- sigmoid.x
mu <- constant.x

# lambdas
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaLoFirt10 = 0.13638367349015, 0.136412610857295, 0.13644155805569, 0.136470515067719, 0.136499481875741, 0.136528458462084, 0.136557444809052, 0.13658644089892, 0.136615446713936, 0.136644462236322
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[886:895], collapse=", ") # expectedLambdaLoLast10 = 0.163355537763678, 0.163384553286064, 0.16341355910108, 0.163442555190948, 0.163471541537916, 0.163500518124259, 0.163529484932281, 0.16355844194431, 0.163587389142705, 0.16361632650985

paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.136361976927131, 0.136369208498794, 0.136376440686558, 0.13638367349015, 0.136390906909294, 0.136398140943716, 0.136405375593141, 0.136412610857295, 0.136419846735901, 0.136427083228686
paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[3574:3583], collapse=", ") # expectedLambdaHiLast10 = 0.163572916771314, 0.163580153264099, 0.163587389142705, 0.163594624406859, 0.163601859056284, 0.163609093090706, 0.16361632650985, 0.163623559313442, 0.163630791501206, 0.163638023072869

# mus
paste(do.call(mu, c(list(x.lo), death))[1:10], collapse=", ") # expectedMuLoFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.lo), death))[886:895], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03

paste(do.call(mu, c(list(x.hi), death))[1:10], collapse=", ") # expectedMuHiFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.hi), death))[3574:3583], collapse=", ") # expectedMuHiLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03




# (7) testInitializationOfTips

if (getRversion() >= "3.6.0") {
    RNGkind(sample.kind = "Rounding")
}

drift <- 0.0
diffusion <- 0.001
lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)

char <- make.brownian.with.drift(drift, diffusion)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=FALSE)
pars <- c(.1, .2, 0, 2.5, .03, 0, .001) # 6th and 7th elements are drift and diffusion
sd <- 1/20
control.C.1 <- list(dt.max=0.01) # dt = 0.01
phy$tip.state[1] <- 0.0
phy$tip.state[2] <- 0.1

cache <- make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)
cache$height[1] <- cache$height[2] <- 0.01 # setting things manually
cache$depth[3] <- 0.01
cache$edge.length <- c(0.01, 0.01)
cache$len <- c(0.01, 0.01, NA)

## these four lines are to make this cache match the previous unit tests
cache$args$drift <- 6
cache$args$diffusion <- 7
cache$control$dx <- 0.0005 ## to match unit test
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
        ext <- quasse.extent.debug(cache$control, drift, diffusion)
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

# we can compare it to running the function as is
cache$y <- initial.tip.quasse(cache, cache$control, pars2[[1]]$x) # this is a list

# cache$y$t # branch lengths
# cache$y$target # tip indices (I think)
identical(cache$y$y$sp1, sp1.y) # TRUE
identical(cache$y$y$sp2, sp2.y) # TRUE
identical(cache$y$y$sp1, sp2.y) # FALSE (just as a control)

sp1.y.ds <- sp1.y[nx+1:(length(sp1.y)-nx)]
sp2.y.ds <- sp2.y[nx+1:(length(sp2.y)-nx)]

sp1.y.ds.exp <- sp1.y.ds[2001:2048]
sp2.y.ds.exp <- sp2.y.ds[2001:2048]

paste(sp1.y.ds.exp, collapse=", ") # 6.96077358436849, 6.95166528584696, 6.94252551477923, 6.93335442671583, 6.92415217758908, 6.91491892370871, 6.90565482175746, 6.89636002878667, 6.88703470221182, 6.87767899980818, 6.86829307970631, 6.85887710038768, 6.84943122068018, 6.83995559975374, 6.83045039711584, 6.82091577260705, 6.81135188639661, 6.80175889897795, 6.79213697116422, 6.78248626408384, 6.772806939176, 6.76309915818623, 6.75336308316186, 6.74359887644761, 6.73380670068105, 6.72398671878815, 6.71413909397875, 6.70426398974212, 6.69436156984244, 6.68443199831428, 6.67447543945815, 6.66449205783599, 6.65448201826663, 6.64444548582134, 6.63438262581928, 6.62429360382306, 6.61417858563415, 6.60403773728847, 6.59387122505179, 6.58367921541529, 6.57346187509105, 6.5632193710075, 6.55295187030495, 6.54265954033109, 6.53234254863645, 6.52200106296994, 6.51163525127429, 6.50124528168164
paste(sp2.y.ds.exp, collapse=", ") # 2.67859021074856, 2.68849414736158, 2.69841783805887, 2.70836123148143, 2.71832427571076, 2.72830691826807, 2.73830910611356, 2.74833078564564, 2.75837190270027, 2.76843240255022, 2.77851222990441, 2.78861132890721, 2.79872964313779, 2.80886711560949, 2.81902368876922, 2.82919930449678, 2.83939390410431, 2.84960742833573, 2.85983981736613, 2.87009101080125, 2.88036094767694, 2.89064956645866, 2.90095680504094, 2.91128260074695, 2.92162689032799, 2.93198960996305, 2.9423706952584, 2.95277008124712, 2.96318770238874, 2.97362349256884, 2.9840773850987, 2.9945493127149, 3.00503920757903, 3.01554700127736, 3.02607262482052, 3.03661600864323, 3.04717708260405, 3.05775577598507, 3.06835201749176, 3.07896573525268, 3.0895968568193, 3.10024530916586, 3.11091101868917, 3.12159391120842, 3.13229391196515, 3.14301094562307, 3.15374493626797, 3.16449580740766



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

control.fft.48 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
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
diffusion <- 0.001
sd <- 0.05

args.fft.just.t <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.t <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.t <- quasse.extent.debug(control.fft.48, drift, diffusion) # prepares X axis stuff
ext.fft.just.t$padding

pars.fft.just.t <- expand.pars.quasse(lambda, mu, args.fft.just.t, ext.fft.just.t, pars.just.t) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.t <- matrix(0, control.fft.48$nx, 2) # low resolution
vars.fft.just.t[seq_len(ext.fft.just.t$ndat[2]),2] <- dnorm(ext.fft.just.t$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.t[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0

vars.fft.just.t.tmp <- fftR.propagate.t(vars.fft.just.t, pars.fft.just.t$lo$lambda, pars.fft.just.t$lo$mu, control.fft.48$dt.max, ext.fft.just.t$ndat[2])
paste(vars.fft.just.t.tmp[,1], collapse=", ") # 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(vars.fft.just.t.tmp[,2], collapse=", ") # 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0

pde.fftR.just.t <- with(control.fft.48, make.pde.quasse.fftR.2(nx, dx, dt.max, 2L))
ans.fftR.just.t <- pde.fftR.just.t(vars.fft.just.t, control.fft.48$dt.max, pars.fft.just.t$lo, 0)

## same as above
paste(ans.fftR.just.t[[2]][,1], collapse=", ") # 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.000299744087825142, 0.000299743181686865, 0.000299742271619522, 0.000299741357883741, 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0.000299719056361739, 0.000299718142720333, 0.000299717232754363, 0.000299716326724287, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.just.t[[2]][,2], collapse=", ") #  0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.0884858018341865, 0.158038086959484, 0.271192794627236, 0.447118602442875, 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052117, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0.447052058075915, 0.271149126249059, 0.158010719780681, 0.0884694089104536, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX

args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x <- quasse.extent.debug(control.fft.48, drift, diffusion) # prepares X axis stuff
ext.fft.just.x$padding # 16/16, 4/4

pars.fft.just.x <- expand.pars.quasse(lambda, mu, args.fft.just.x, ext.fft.just.x, pars.just.x) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x <- matrix(0, control.fft.48$nx, 2) # low resolution
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[2]),2] <- dnorm(ext.fft.just.x$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.x[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0

# checking kernel matches (see PropagatesQuaSSETest -> testMakeNormalKernInPlaceAndFFtAndIfft)
kern.just.x <- fftR.make.kern(-control.fft.48$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.48$dt * pars.fft.just.x$hi$diffusion), control.fft.48$nx, control.fft.48$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2])
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

ds.prop.x <- fftR.propagate.x(vars.fft.just.x, control.fft.48$nx, fy.just.x, pars.fft.just.x$lo$padding[1], pars.fft.just.x$lo$padding[2])
paste(ds.prop.x[,1], collapse=", ") # 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ds.prop.x[,2], collapse=", ") # 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0888278883396178, 0.158599420774613, 0.272077439492024, 0.44845817681081, 0.710214708266689, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804265, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661662, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.08067601406491, 0.710214708266689, 0.448458176810809, 0.272077439492024, 0.158599420774613, 0.088827888339617, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX (and 10.1 testPropagateChOneCh1024QuaSSETest)

control.fft.just.x.1024 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=1024, # number of X bins
                    dx=0.0005, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

sd <- 0.05
drift <- 0.0
diffusion <- 0.001
lambda <- sigmoid.x
mu <- constant.x

args.fft.just.x.1024 <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x.1024 <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x.1024 <- quasse.extent.debug(control.fft.just.x.1024, drift, diffusion) # prepares X axis stuff
# ext.fft.just.x.1024$padding

pars.fft.just.x.1024 <- expand.pars.quasse(lambda, mu, args.fft.just.x.1024, ext.fft.just.x.1024, pars.just.x.1024) # adds lambda and mu vectors to ext.fft.just.x

## initialization
vars.fft.just.x.1024 <- matrix(0, control.fft.just.x.1024$nx, 2) # low resolution
vars.fft.just.x.1024[seq_len(ext.fft.just.x.1024$ndat[2]),2] <- dnorm(ext.fft.just.x.1024$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.1024$lo$padding should be 0.0
paste(vars.fft.just.x.1024[1:48,2], collapse=", ") # 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

kern.just.x.1024 <- fftR.make.kern(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$diffusion), control.fft.just.x.1024$nx, control.fft.just.x.1024$dx, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])

fy.just.x.1024 <- fft(kern.just.x.1024)
paste(Re(fy.just.x.1024)[1:48], collapse=", ") # 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

## checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x.1024, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x.1024
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints (in Dell XPS + Ubuntu, this differs!)

ds.prop.x.1024 <- fftR.propagate.x(vars.fft.just.x.1024, control.fft.just.x.1024$nx, fy.just.x.1024, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])
paste(ds.prop.x.1024[1:48,2], collapse=", ") # 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

paste(ds.prop.x.1024[848:895,2], collapse=", ") # 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948

## pde.fftC <- with(control.fft.just.x.1024$nx, diversitree:::make.pde.mosse.fftC(control.fft.just.x.1024$nx, control.fft.just.x.1024$nx*control.fft.just.x.1024$r, control.fft.just.x.1024$dt.max, 4L, flags))
## ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0, control.fft$dt.max)



## (11) testIntegrateOneBranchLoRes4096BinsOutsideClassJustX (note that I'm using 4096 as low res here, but I'll co-opt 4096 high res in Java)  (and 11.1 testPropagateChOneCh1024QuaSSETest)

## initialization
vars.fft.just.x.4096 <- matrix(0, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, 2) # low resolution
vars.fft.just.x.4096[seq_len(ext.fft.just.x.1024$ndat[1]),2] <- dnorm(ext.fft.just.x.1024$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.4096$lo$padding should be 0.0
paste(vars.fft.just.x.4096[1:48,2], collapse=", ") # 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846

kern.just.x.4096 <- fftR.make.kern(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$diffusion), control.fft.just.x.1024$nx * control.fft.just.x.1024$r, control.fft.just.x.1024$dx, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])

fy.just.x.4096 <- fft(kern.just.x.4096)
paste(Re(fy.just.x.4096)[1:48], collapse=", ") # 1, 0.999952939166244, 0.999811769952891, 0.999576532217434, 0.999247292368191, 0.998824143333056, 0.998307204515792, 0.997696621739869, 0.996992567179915, 0.996195239280797, 0.995304862664403, 0.994321688024177, 0.993245992007481, 0.992078077085863, 0.990818271413312, 0.989466928672585, 0.988024427909734, 0.986491173356906, 0.984867594243561, 0.983154144596211, 0.981351303026827, 0.979459572510042, 0.977479480149295, 0.975411576932071, 0.973256437474405, 0.971014659754794, 0.968686864837709, 0.966273696586883, 0.963775821368548, 0.961193927744829, 0.958528726157476, 0.955780948602156, 0.952951348293497, 0.950040699321108, 0.947049796296786, 0.943979453993153, 0.940830506973933, 0.937603809216109, 0.934300233724213, 0.930920672136968, 0.92746603432656, 0.92393724799076, 0.920335258238176, 0.916661027166888, 0.912915533436717, 0.909099771835414, 0.905214752839018, 0.90126150216667

ds.prop.x.4096 <- fftR.propagate.x(vars.fft.just.x.4096, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, fy.just.x.4096, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])
paste(ds.prop.x.4096[1:48,2], collapse=", ") # 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846
paste(ds.prop.x.4096[1001:1048,2], collapse=", ") # 1.12964601442883, 1.13523957985681, 1.14085371385815, 1.14648844781881, 1.1521438128928, 1.15781984000018, 1.16351655982519, 1.16923400281429, 1.17497219917421, 1.18073117887008, 1.18651097162341, 1.19231160691021, 1.19813311395904, 1.20397552174904, 1.20983885900802, 1.21572315421049, 1.22162843557575, 1.22755473106588, 1.23350206838385, 1.23947047497152, 1.24545997800775, 1.25147060440637, 1.25750238081428, 1.26355533360948, 1.26962948889911, 1.27572487251749, 1.28184151002416, 1.28797942670193, 1.29413864755491, 1.30031919730656, 1.30652110039772, 1.31274438098464, 1.31898906293703, 1.32525516983611, 1.33154272497262, 1.33785175134486, 1.34418227165675, 1.35053430831583, 1.35690788343134, 1.36330301881221, 1.36971973596514, 1.37615805609261, 1.38261800009092, 1.38909958854822, 1.39560284174258, 1.40212777963999, 1.40867442189243, 1.41524278783588



## (12) testBothXandTPropagateMethodsInsideClassOneBranchLoRes48Bins

## need to control.fft.48 from (8)

## version for debugging
quasse.integrate.fftR.3 <- function (vars, lambda, mu, drift, diffusion, nstep, dt, nx,
    ndat, dx, nkl, nkr) {

    print("vars"); print(vars)
    print("lambda"); print(lambda)
    print("mu"); print(mu)
    print("drift"); print(drift)
    print("diffusion"); print(diffusion)
    print("nstep"); print(nstep)
    print("dt"); print(dt)
    print("nx"); print(nx)
    print("dx"); print(dx)
    print("nkl"); print(nkl)
    print("nkr"); print(nkr)

    kern = diversitree:::fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
        nx, dx, nkl, nkr)

    fy = fft(kern)

    print("fy")
    print(paste0(Re(fy), collapse=", "))

    for (i in seq_len(nstep)) {
        print("E's before propagate in t")
        print(paste(vars[,1], collapse=", "))
        print("D's before propagate in t")
        print(paste(vars[,2], collapse=", "))

        vars = diversitree:::fftR.propagate.t(vars, lambda, mu, dt, ndat)

        print("E's after propagate in t")
        print(paste(vars[,1], collapse=", "))
        print("D's after propagate in t")
        print(paste(vars[,2], collapse=", "))

        vars = diversitree:::fftR.propagate.x(vars, nx, fy, nkl, nkr)

        print("E's after propagate in t and x")
        print(paste(vars[,1], collapse=", "))
        print("D's after propagate in t and x")
        print(paste(vars[,2], collapse=", "))
    }

    vars
}

## version for debugging
make.pde.quasse.fftR.3 <- function (nx, dx, dt.max, nd) {
    function(y, len, pars, t0) {
        padding <- pars$padding

        ndat <- length(pars$lambda)
        nt <- as.integer(ceiling(len/dt.max)) # number of time intervals, rounded up
        dt <- len/nt # size of the interval;

        ## branch length = 10
        ## dt.max = 0.45
        ## 0.1 length left after 2 * dt.max
        ## when he ceiling, he gets nt = 3
        ## then dt = 10 / 3 = 0.33333333
        ## then there's no left over, and dt is by definition < dt.max

        print("y"); print(y)
        print("length(y)"); print(length(y))
        print("nd"); print(nd)
        print("nd * nx"); print(nd * nx)
        print("ndat"); print(ndat)

        if (!(length(y) %in% (nd * nx)))
            stop("Wrong size y")
        if (length(pars$lambda) != length(pars$mu) || length(pars$lambda) >
            (nx - 3))
            stop("Incorrect length pars")
        if (pars$diffusion <= 0)
            stop("Invalid diffusion parameter")
        if (!is.matrix(y))
            y <- matrix(y, nx, nd)

        ## y comes in normalized from combine.branches.quasse
        ## line: y[, 2] <- y[, 2]/q0 (the y here is normalized by q0 in that function)
        ans = quasse.integrate.fftR.3(y, pars$lambda, pars$mu,
            pars$drift, pars$diffusion, nt, dt, nx, ndat, dx,
            padding[1], padding[2])

        q = sum(ans[, 2]) * dx

        print("Normalization factor in make.pde.quasse.fftR = ")
        print(q)

        # print(ans[,1]) # E's
        # print(ans[,2]) ## FKM: D's that are returned from fftR.propagate.t

        ans[,2] <- ans[,2]/q ## FKM: D's are normalized again before returning

        list(log(q), ans)
    }
}

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.48.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.48.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.48.both <- quasse.extent.debug(control.fft.48, drift, diffusion) # prepares X axis stuff

vars.fft.48.both <- matrix(0, control.fft.48$nx, 2) # low resolution
vars.fft.48.both[seq_len(ext.fft.48.both$ndat[2]),2] <- dnorm(ext.fft.48.both$x[[2]], 0.0, sd)
paste(c(vars.fft.48.both[seq_len(ext.fft.48.both$ndat[2]),2], rep(0, (48-length(ext.fft.48.both$x[[2]])))), collapse=", ") # expectedInitialDs = 0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694604, 0.0122380386022755, 0.0058389385158292, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.48.both <- expand.pars.quasse(lambda, mu, args.fft.48.both, ext.fft.48.both, pars.48.both)

pde.fftR.48.both <- with(control.fft.48, make.pde.quasse.fftR.3(nx, dx, dt.max, 2L)) # creates function

ans.fftR.48.both <- pde.fftR.48.both(vars.fft.48.both, control.fft.48$dt.max, pars.fft.48.both$lo, 0)

## fy = 1, 0.999886244673461, 0.999546925086093, 0.998987847110852, 0.998218576759891, 0.997252276505192, 0.996105480062091, 0.994797809489411, 0.993351639446935, 0.991791714355113, 0.990144725007753, 0.988438851882212, 0.986703282961364, 0.984967714317709, 0.983261842004857, 0.981614853950298, 0.980054930543287, 0.978608762462816, 0.977301093995625, 0.976154299658014, 0.97518800136532, 0.97441873269917, 0.97385965601673, 0.973520337242051, 0.973406582192705, 0.973520337242051, 0.97385965601673, 0.97441873269917, 0.97518800136532, 0.976154299658014, 0.977301093995625, 0.978608762462816, 0.980054930543287, 0.981614853950298, 0.983261842004857, 0.984967714317709, 0.986703282961364, 0.988438851882212, 0.990144725007753, 0.991791714355113, 0.993351639446935, 0.994797809489411, 0.996105480062091, 0.997252276505192, 0.998218576759891, 0.998987847110852, 0.999546925086093, 0.999886244673461

## E's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)
paste(ans.fftR.48.both[[2]][,1], collapse=", ") # 0.00029974766804671, 0.000299746780132305, 0.000299745887296174, 0.000299744989778572, 0.00029974408779732, 0.000299743181660744, 0.000299742271595132, 0.000299741357861114, 0.00029974044072366, 0.000299739520451864, 0.000299738597318595, 0.000299737671600252, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507617, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742340082, 0.000299721816686999, 0.00029972089362645, 0.000299719973434657, 0.000299719056384414, 0.000299718142744768, 0.00029971723278053, 0.000299716326752154, 0.000299715424885881, 0.000299714527489882, 0.000299713634781839, 0.00029971274700187, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)
paste(ans.fftR.48.both[[2]][,2], collapse=", ") # 0.00582911973804044, 0.0122173866829305, 0.0246026489190146, 0.0476007314229174, 0.08867639188041, 0.15832797168282, 0.271610119667664, 0.447685177240542, 0.708986186416951, 1.07880005021062, 1.57718311562593, 2.21544576526305, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066457, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21522515007474, 1.57700658374745, 1.07866601799939, 0.708889397846809, 0.447618584199003, 0.271566407577905, 0.158300569088559, 0.0886599725440681, 0.0475913399553493, 0.0245975002358219, 0.0122146843473713, 0.00582776134335783, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (13) testBothXandTPropagateMethodsInsideClassOneBranchLoRes48BinsLargerDt

## need to define quasse.integrate.fftR.3 and make.pde.quasse.fftR.3 from (12)

control.fft.48.2dt <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.02, # dt
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
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.48.2dt.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.48.2dt.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.48.2dt.both <- quasse.extent.debug(control.fft.48.2dt, drift, diffusion) # prepares X axis stuff

vars.fft.48.2dt.both <- matrix(0, control.fft.48.2dt$nx, 2) # low resolution
vars.fft.48.2dt.both[seq_len(ext.fft.48.2dt.both$ndat[2]),2] <- dnorm(ext.fft.48.2dt.both$x[[2]], 0.0, sd)
paste(c(vars.fft.48.2dt.both[seq_len(ext.fft.48.2dt.both$ndat[2]),2], rep(0, (48-length(ext.fft.48.2dt.both$x[[2]])))), collapse=", ") # expectedInitialDs = 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.1583090316596, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0.447890605896858, 0.271659384673712, 0.158309031659599, 0.08863696823876, 0.0476817640292968, 0.0246443833694603, 0.0122380386022755, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.48.2dt.both <- expand.pars.quasse(lambda, mu, args.fft.48.2dt.both, ext.fft.48.2dt.both, pars.48.2dt.both)

pde.fftR.48.2dt.both <- with(control.fft.48.2dt, make.pde.quasse.fftR.3(nx, dx, dt.max, 2L)) # creates function

ans.fftR.48.2dt.both <- pde.fftR.48.2dt.both(vars.fft.48.2dt.both, control.fft.48.2dt$dt.max, pars.fft.48.2dt.both$lo, 0)

## fy = 1, 0.998791000226283, 0.995184822859046, 0.989243568247244, 0.981069525480477, 0.970803379084231, 0.958621745686731, 0.944734086634627, 0.92937905384229, 0.91282033618305, 0.895342082276453, 0.87724398244744, 0.858836097834951, 0.840433528062905, 0.822351010557628, 0.80489754454845, 0.788371131106066, 0.773053717371028, 0.759206428540079, 0.747065165363734, 0.736836638024787, 0.728694899474875, 0.722778432760492, 0.719187837716992, 0.717984152783717, 0.719187837716992, 0.722778432760492, 0.728694899474875, 0.736836638024787, 0.747065165363734, 0.759206428540079, 0.773053717371028, 0.788371131106066, 0.80489754454845, 0.822351010557628, 0.840433528062905, 0.858836097834951, 0.87724398244744, 0.895342082276454, 0.91282033618305, 0.92937905384229, 0.944734086634627, 0.958621745686731, 0.970803379084231, 0.981069525480477, 0.989243568247244, 0.995184822859046, 0.998791000226283

## E's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)
## E's after t = 0.000598987856474924, 0.000598984289866083, 0.000598980704570308, 0.000598977101569278, 0.000598973481865553, 0.000598969846481545, 0.000598966196458183, 0.000598962532854178, 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940327231526, 0.000598936599425321, 0.000598932866991505, 0.000598929131089762, 0.000598925392884107, 0.000598921653541243, 0.000598917914229581, 0.000598914176117241, 0.000598910440370952, 0.000598906708154486, 0.000598902980627182, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0.000598880777947923, 0.000598877115040413, 0.000598873465770468, 0.000598869831195672, 0.000598866212356192, 0.000598862610273333, 0.000598859025948505, 0.000598855460362313, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.48.2dt.both[[2]][,1], collapse=", ") # after t and x = 0.000598987856474924, 0.000598984289866083, 0.000598980704570308, 0.000598977101569278, 0.000598973481865553, 0.00059896984544713, 0.000598966195498576, 0.000598962531970537, 0.00059895885593801, 0.000598955168491078, 0.000598951470733874, 0.000598947763783104, 0.000598944048766786, 0.000598940326823012, 0.000598936599098334, 0.000598932866746461, 0.000598929130926968, 0.000598925392803752, 0.000598921653543447, 0.000598917914314325, 0.000598914176284426, 0.000598910440620369, 0.000598906708485822, 0.000598902981040027, 0.000598899259436428, 0.000598895544821229, 0.000598891838331867, 0.00059888814109588, 0.000598884454229337, 0.000598880778835635, 0.000598877116004042, 0.000598873466808826, 0.000598869831195672, 0.000598866212356192, 0.000598862610273333, 0.000598859025948505, 0.000598855460362313, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)
## D's after t = 0.0121967797643638, 0.0245610056719514, 0.0475198764133255, 0.0883349677000203, 0.15776773954185, 0.270727236150491, 0.446348310766164, 0.707040094708608, 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82233340368848, 5.77330938669844, 6.64080371901192, 7.33913154877464, 7.79286078082761, 7.95018769794391, 7.79266608851656, 7.33876489743091, 6.64030620874454, 5.77273291052126, 4.82173179361502, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0.706847016883081, 0.446215500932528, 0.270640082416696, 0.157713119885317, 0.0883022505584118, 0.0475011328396031, 0.0245507298846628, 0.0121913864192371, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.48.2dt.both[[2]][,2], collapse=", ") # after t and x = 0.0121967797643638, 0.0245610056719514, 0.0475198764133255, 0.0883349677000203, 0.15776773954185, 0.275155027270175, 0.452359435761905, 0.714695461673121, 1.08514694006539, 1.58338053088545, 2.22029371795898, 2.99201095412974, 3.87474275820035, 4.82224125131593, 5.76741044168543, 6.62885082364537, 7.32184914756387, 7.77191831137883, 7.92794195899393, 7.77172522522963, 7.32148540005609, 6.62835697979005, 5.76683776896088, 4.8216430122779, 3.8741662641039, 2.9914919608128, 2.21985391738919, 1.5830280032418, 1.08487876446248, 0.714501388583886, 0.4522255935314, 0.275066946374052, 0.157713119885317, 0.0883022505584118, 0.0475011328396031, 0.0245507298846628, 0.0121913864192371, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (14) testIntegrateOneBranchHiRes48BinsInsideClassBothXandTTwoDt

## need to control.fft.48.2dt from (13)
## need to define quasse.integrate.fftR.3 and make.pde.quasse.fftR.3 from (12)

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.48.2dt.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.48.2dt.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.48.2dt.both <- quasse.extent.debug(control.fft.48.2dt, drift, diffusion) # prepares X axis stuff

vars.fft.48.2dt.both <- matrix(0, control.fft.48.2dt$nx * 4, 2) # high resolution
vars.fft.48.2dt.both[seq_len(ext.fft.48.2dt.both$ndat[1]),2] <- dnorm(ext.fft.48.2dt.both$x[[1]], 0.0, sd)

print(paste(vars.fft.48.2dt.both[1:10,2], collapse=", ")) # expectedInitialDsFirst10 = 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755, 0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349, 0.0344513787810736
print(paste(vars.fft.48.2dt.both[147:156,2], collapse=", ")) # expectedInitialDsLater10 = 0.0146332892566062, 0.0122380386022755, 0.0102092994868837, 0.00849560541101504, 0.00705191364734891, 0, 0, 0, 0, 0

## have to normalize it to match behavior inside likelihood class
vars.fft.48.2dt.both[,2] <- vars.fft.48.2dt.both[,2] / (sum(vars.fft.48.2dt.both[,2]) * 0.01 / 4) #
print(paste(vars.fft.48.2dt.both[1:10,2], collapse=", ")) # 0.00705304040908636, 0.00849696284725042, 0.0102109307388535, 0.0122399940081877, 0.0146356273780371, 0.0174564426629985, 0.0207689438654755, 0.0246483210734104, 0.0291792776423502, 0.0344568834564503
print(paste(vars.fft.48.2dt.both[147:156,2], collapse=", ")) # 0.0146356273780371, 0.0122399940081877, 0.0102109307388535, 0.00849696284725042, 0.00705304040908636, 0, 0, 0, 0, 0

pars.fft.48.2dt.both <- expand.pars.quasse(lambda, mu, args.fft.48.2dt.both, ext.fft.48.2dt.both, pars.48.2dt.both)

pde.fftR.48.2dt.both <- with(control.fft.48.2dt, make.pde.quasse.fftR.3(48 * 4, dx, dt.max, 2L)) # note the 48 * 4 to do high res

ans.fftR.48.2dt.both <- pde.fftR.48.2dt.both(vars.fft.48.2dt.both, control.fft.48.2dt$dt.max, pars.fft.48.2dt.both$hi, 0)
ans.fftR.48.2dt.both[[2]][,2] <- ans.fftR.48.2dt.both[[2]][,2] * exp(ans.fftR.48.2dt.both[[1]]) # unnormalizing it from what make.pde.quasse.fftR.3 does
ans.fftR.48.2dt.both[[2]][,2] <- ans.fftR.48.2dt.both[[2]][,2] / (sum(ans.fftR.48.2dt.both[[2]][,2]) * 0.0025) # re-normalize it to match the java class

## E's
paste(ans.fftR.48.2dt.both.hi[[2]][1:10,1], collapse=", ") # first 10 after t and x = 0.000598990518590012, 0.000598989632470972, 0.000598988745094367, 0.000598987856474924, 0.000598986966627531, 0.000598986075566877, 0.000598985183308093, 0.000598984289866083, 0.000598983395256048, 0.000598982499493234
paste(ans.fftR.48.2dt.both.hi[[2]][147:156,1], collapse=", ") # first 10 after t and x = 0.000598856349949292, 0.000598855460362313, 0.000598854572006534, 0.000598853684896709, 0.000598852799047445, 0, 0, 0, 0, 0

## D's
paste(ans.fftR.48.2dt.both[[2]][1:10,2], collapse=", ") # expectedSp1DsAfterPropTandXFirst10 = 0.00705472178566115, 0.00849896328584349, 0.0102133044244718, 0.0122428030418449, 0.0146389426900442, 0.0174603449881642, 0.0207735247765771, 0.0246536840674864, 0.0291855392723177, 0.0344641744982804
paste(ans.fftR.48.2dt.both[[2]][147:156,2], collapse=", ") # expectedSp1DsAfterPropTandXLater10 = 0.0146325564177363, 0.0122373893454924, 0.0102087275584534, 0.0084951043048786, 0.0070514768252866, 0, 0, 0, 0, 0



## (15)

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=TRUE)
# simplifying tree for unit test
tr$tip.state[1] <- 0.0 # ch state
tr$tip.state[2] <- 0.1 # ch state
tr$edge.length <- c(0.01, 0.01)  # branch lengths
tr$orig[,2] <- c(0.01, 0.01) # branch lengths
tr$orig[,4] <- c(0.0, 0.1) # ch states

pars <- c(.1, .2, 0, 2.5, .03, drift, diffusion)

make.cache.quasse.debug <- function (tree, states, states.sd, lambda, mu, control, sampling.f, for.split = FALSE) {
    tree <- diversitree:::check.tree(tree)
    tmp <- diversitree:::check.states.quasse(tree, states, states.sd)
    states <- tmp$states
    states.sd <- tmp$states.sd
    control <- diversitree:::check.control.quasse(control, tree, states)
    cache <- diversitree:::make.cache(tree)
    cache$states <- states
    cache$states.sd <- states.sd
    cache$control <- control
    if (!for.split) {
        n.lambda <- diversitree:::check.f.quasse(lambda)
        n.mu <- diversitree:::check.f.quasse(mu)
        n.args <- n.lambda + n.mu + 2
        args <- list(lambda = seq_len(n.lambda), mu = seq_len(n.mu) +
            n.lambda, drift = n.lambda + n.mu + 1, diffusion = n.lambda +
            n.mu + 2)
        cache$lambda <- lambda
        cache$mu <- mu
        cache$args <- args
        sampling.f <- diversitree:::check.sampling.f(sampling.f, 1)
        cache$sampling.f <- sampling.f
    }
    cache$info <- diversitree:::make.info.quasse(lambda, mu, tree)
    cache
}

combine.branches.quasse.debug <- function (f.hi, f.lo, control) {

    nx <- control$nx
    dx <- control$dx
    tc <- control$tc
    r <- control$r
    eps <- log(control$eps)
    dt.max <- control$dt.max

    careful <- function(f, y, len, pars, t0, dt.max) {

        ## if low res, f() == f.lo
        ## if high res, f() == f.hi
        ans <- f(y, len, pars, t0)

        if (ans[[1]] > eps) {
            ans
        }
        else {
            if (control$method == "fftC" || control$method ==
                "fftR")
                dt.max <- dt.max/2
            len2 <- len/2
            ans1 <- Recall(f, y, len2, pars, t0, dt.max)
            ans2 <- Recall(f, ans1[[2]], len2, pars, t0 + len2,
                dt.max)
            ans2[[1]][[1]] <- ans1[[1]][[1]] + ans2[[1]][[1]]
            ans2
        }
    }

    function(y, len, pars, t0, idx) {

        print(paste0("len = ", len))

        if (t0 < tc) {
            dx0 <- dx/r
            nx0 <- nx * r

            ## print(paste("dx0 = ", dx0))
            ## print(paste("nx0 = ", nx0))

        }
        else {
            dx0 <- dx
            nx0 <- nx
        }
        if (any(y < -1e-08))
            stop("Actual negative D value detected -- calculation failure")

        y[y < 0] <- 0
        y <- matrix(y, nx0, 2) # if in low-res, this re-sizes it

        print("y inside combine.branches before normalization = ")
        print(y)

        q0 <- sum(y[, 2]) * dx0 # calculating normalization factor from D's
        if (q0 <= 0)
            stop("No positive D values")

        y[, 2] <- y[, 2]/q0 # normalizing D's, before we integrate (NOT LOG-fied!!!!)

        print("y inside combine.branches after normalization = ")
        print(paste(y[,2], collapse=", "))

        lq0 <- log(q0) # will be added at the very end
        print(paste0("log-normalization factor starts at ", log(q0)))

        ## option 1: initial time is already past tc, low-res
        if (t0 >= tc) {

            print("y inside combine.branches for lo = ")
            print(y)

            ## ans[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans <- careful(f.lo, y, len, pars$lo, t0, dt.max) # integrate

            print("ans at lo, prior to adding log-normalization factor")
            print(ans)
        }

        ## option 2: initial time + dt is still lower than tc, high-res
        else if (t0 + len < tc) {

            print("y inside combine.branches for hi = ")
            print(y)

            ans <- careful(f.hi, y, len, pars$hi, t0, dt.max) # integrate

            print("option 2: ans at hi, prior to adding log-normalization factor")
            print(ans)
        }

        ## option 3: first segment in high-resolution (up until tc), then remaining in low-res
        ## note that option 3 will happen before option 1 -- in option 1, y has already been
        ## resized to low-res
        else {

            print("option 3: y inside combine.branches for hi = ")
            print(y)

            len.hi <- tc - t0 # this first stretch in high-res

            print(paste0("len.hi = ", len.hi))

            ## ans.hi[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans.hi <- careful(f.hi, y, len.hi, pars$hi, t0, dt.max) # integrate first stretch at high-res

            ## but then here we grab just every r-th element in ans.hi (as given by a vector of indices separated by r in pars$tr)
            ## y.lo will be the first elements inside the low-res y -- the rest will be 0's (see below)
            y.lo <- ans.hi[[2]][pars$tr, ]

            ## I believe ans.hi[[1]] here was normalized and log-fied inside the function that generated ans.hi
            lq0 <- lq0 + ans.hi[[1]] # we add the high-res normalization factor to the total normalization factor
            print(paste0("added ", ans.hi[[1]], " to log-normalization factor"))
            print(paste0("log-normalization factor is now ", lq0))

            ## why would nrow(y.lo) < nx (nx here is the low-res nx)?
            ## we're putting 0's after the y.lo
            print("pars$tr = ")
            print(paste(pars$tr, collapse=", "))
            print("expectedHiLoIdxs4Transfer = ")
            print(paste(pars$tr-1, collapse=", "))
            if (nrow(y.lo) < nx) y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 2))

            print(paste0("len - len.hi = ", len - len.hi))

            ## ans[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans <- careful(f.lo, y.lo, len - len.hi, pars$lo, tc, dt.max) # integrate remaining stretch at low-res

            print("ans at hi, but interval < dt, prior to adding log-normalization factor")
            print(ans)
        }

        print("ans at return of combine.branches.quasse, prior to adding log-normalization factor")
        print(ans)

        print(paste0("log-normalization factor = ", lq0))

        print("ans at return of combine.branches.quasse, after adding log-normalization factor")
        print(c(ans[[1]] + lq0, ans[[2]])) # lq0 here can be either high or low-res depending on which of the if/else if/else block executed, i.e., whether the branch is < >= tc

        ## ans[[1]] is the log((sum over all D's) * dx), the input of which was normalized; we need to keep adding to it; so at the very end of the pruning (at the root) we can unnormalize the likelihood
        ## ans[[2]] is the dataframe with the E's and D's
        c(ans[[1]] + lq0, ans[[2]])
    }
}

make.branches.quasse.fftR.debug <- function(control) {
    nx <- control$nx
    dx <- control$dx
    r <- control$r
    dt.max <- control$dt.max
    tc <- control$tc

    ## make.pde.quasse.fftR returns a function, which in turn runs quasse.integrate.fftR, returning c(log(q), ans)
    f.hi <- make.pde.quasse.fftR.3(nx * r, dx/r, dt.max, 2L)
    f.lo <- make.pde.quasse.fftR.3(nx, dx, dt.max, 2L)

    combine.branches.quasse.debug(f.hi, f.lo, control)
}

initial.tip.quasse.debug <- function(cache, control, x) {
    nx <- control$nx * control$r
    npad <- nx - length(x)
    e0 <- 1 - cache$sampling.f
    if (control$tips.combined) {
        tips <- cache$tips
        t <- cache$len[tips]
        i <- order(t)
        target <- tips[i]
        states <- cache$states[i]
        states.sd <- cache$states.sd[i]
        y <- mapply(function(mean, sd) c(dnorm(x, mean, sd),
            rep(0, npad)), states, states.sd, SIMPLIFY = FALSE)
        y <- matrix(c(rep(e0, nx), unlist(y)), nx, length(target) + 1)

        ## print("if: y inside initial.tip.quasse = ");  print(y)

        list(target = target, y = y, t = t[i])
    }
    else {

        print("cache$states = "); print(cache$states)
        print("cache$states.sd = "); print(cache$states.sd)

        y <- mapply(function(mean, sd) c(rep(e0, nx), dnorm(x,
            mean, sd), rep(0, npad)), cache$states, cache$states.sd,
            SIMPLIFY = FALSE)

        print("else: y inside initial.tip.quasse = ");
        print(paste(y, collapse=", "))
        ## print("esDsHiAtNodeInitial0[1] = ")
        ## print(paste(y$sp1[129:231], collapse=", "))
        ## print("esDsHiAtNodeInitial1[1] = ")
        ## print(paste(y$sp2[129:231], collapse=", "))

        diversitree:::dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
    }
}

all.branches.list.debug <- function (pars, cache, initial.conditions, branches, preset) {
    len <- cache$len
    depth <- cache$depth
    children <- cache$children
    order <- cache$order[-length(cache$order)]
    root <- cache$root
    n <- length(len)
    lq <- rep(0, n)
    n.tip <- cache$n.tip
    y <- cache$y

    branch.init <- branch.base <- vector("list", n)

    ## not sure what this is for
    if (!is.null(preset)) {
        lq[preset$target] <- preset$lq
        branch.base[preset$target] <- preset$base
    }

    ## doing all tips first
    if (is.null(names(y))) {
        ## x is a tip
        for (x in y) {
            if (!is.null(x)) {
                idx <- x$target
                branch.init[idx] <- list(x$y)
                ans <- branches(x$y, x$t, pars, 0, idx) # return from combine.branches.quasse()
                lq[idx] <- unlist(lapply(ans, "[[", 1)) # getting the log((sum over all D's) * dt) for this node
                branch.base[idx] <- lapply(ans, "[", -1)
            }
        }
    }

    ## tip (not used)
    else {
        tip.t <- y$t
        tip.target <- y$target
        tip.y <- branch.init[tip.target] <- y$y
        for (i in seq_along(tip.t)) {
            idx <- tip.target[i]
            ans <- branches(tip.y[[i]], tip.t[i], pars, 0, idx) # return from combine.branches.quasse()
            lq[idx] <- ans[[1]] # getting the log((sum over all D's) * dt) for this node
            branch.base[[idx]] <- ans[-1]
        }
    }

    ## internal nodes done here
    for (i in order) {
        y.in <- initial.conditions(branch.base[children[i, ]], pars, depth[i], i) # compute initial D for i-th internal node from children
        if (!is.list(y.in) && !all(is.finite(y.in)))
            stop("Bad initial conditions: calculation failure along branches?")
        branch.init[[i]] <- y.in # book keeping, and MLE at each internal node
        ans <- branches(y.in, len[i], pars, depth[i], i) # return from combine.branches.quasse()
        lq[i] <- ans[[1]] # getting the log((sum over all D's) * dt) for this node
        branch.base[[i]] <- ans[-1] # E's and normalized D's
    }

    y.in <- initial.conditions(branch.base[children[root, ]],
        pars, depth[root], root)

    branch.init[[root]] <- y.in

    print(paste0("Log-normalizing factors = ", paste(lq, collapse=", ")))

    ## return with answer we care: lq (normalizing factors for all nodes in the tree) and y.in (normalized D of root)
    list(init = branch.init, base = branch.base, lq = lq, vals = y.in)
}

make.initial.conditions.quasse.debug <- function(control) {
    tc <- control$tc
    r <- control$r
    nx.lo <- control$nx
    nx.hi <- nx.lo * r
    eps <- 1e-08
    function(init, pars, t, idx) {
        if (length(init[[1]]) != length(init[[2]]))
            stop("Data have incompatible length")
        if (t < tc) {
            nx <- nx.hi
            lambda <- pars[[1]]$lambda
        }
        else {
            nx <- nx.lo
            lambda <- pars[[2]]$lambda
        }
        ndat <- length(lambda)
        i <- seq_len(nx)
        j <- seq.int(nx + 1, nx + ndat)

        ## return
        ## init[[1]][i] is the E from child 1 (= E from child 2)
        ## init[[1]][j] is the D from child 1
        ## init[[2]][j] is the D from child 2
        ## and 0's for padding
        print(paste0("nx = ", nx))
        print(paste0("ndat = ", ndat))
        print("merging:")
        print(paste0("E's at merge = ", paste(init[[1]][i], collapse=", ")))
        print(paste0("Left D's at merge = ", paste(init[[1]][j], collapse=", ")))
        print(paste0("Right D's at merge = ", paste(init[[2]][j], collapse=", ")))
        print(paste0("Lambdas at merge = ", paste(lambda, collapse=", ")))
        print(c(init[[1]][i], init[[1]][j] * init[[2]][j] * lambda,
              rep.int(0, nx - ndat)))

        c(init[[1]][i], init[[1]][j] * init[[2]][j] * lambda,
            rep.int(0, nx - ndat))
    }
}

make.all.branches.quasse.debug <- function (cache, control)  {

    branches <- make.branches.quasse.fftR.debug(control)

    initial.conditions <- make.initial.conditions.quasse.debug(control)

    function(pars, intermediates, preset = NULL) {
        cache$y <- initial.tip.quasse.debug(cache, cache$control, pars[[1]]$x)

        ## prunes tree and returns list with likelihood
        all.branches.list.debug(pars, cache, initial.conditions, branches,
            preset)
    }
}

make.pars.quasse.debug <- function (cache) {
    args <- cache$args
    function(pars) {
        names(pars) <- NULL
        drift <- pars[args$drift]
        diffusion <- pars[args$diffusion]
        ext <- quasse.extent.debug(cache$control, drift, diffusion)
        pars <- diversitree:::expand.pars.quasse(cache$lambda, cache$mu, args,
            ext, pars)
        diversitree:::check.pars.quasse(pars$hi$lambda, pars$hi$mu, drift,
            diffusion)
        pars
    }
}

make.quasse.debug <- function (tree, states, states.sd, lambda, mu, control = NULL, sampling.f = NULL) {
    cache <- make.cache.quasse.debug(tree, states, states.sd, lambda, mu, control, sampling.f)

    print("cache$control")
    print(cache$control)

    ## We go function to function to function to function...
    ## make.all.branches.quasse (1) initializes things with make.initial.conditions.quasse(),
    ## and then (2) computes the likelihood with all.branches.list(), which in turn calls make.branches.quasse.fftR()
    ## then make.branches.quasse.fftR() calls combine.branches.quasse()
    ## combine.branches.quasse() makes use of make.pde.quasse.fftR(), which returns a function that actually does the integration
    ## in quasse.integrate.fftR()
    all.branches <- make.all.branches.quasse.debug(cache, cache$control)

    rootfunc <- make.rootfunc.quasse.debug(cache)
    f.pars <- make.pars.quasse.debug(cache)
    ll <- function(pars, condition.surv = TRUE, root = ROOT.OBS,
        root.f = NULL, intermediates = FALSE) {
        pars2 <- f.pars(pars)
        ans <- all.branches(pars2, intermediates) # list containing numbers we care about
        rootfunc(ans, pars2, condition.surv, root, root.f, intermediates)
    }
    class(ll) <- c("quasse", "dtlik", "function")
    ll
}

root.p.quasse.debug <- function(d.root, pars, root, root.f) {
    if (!is.null(root.f) && root != ROOT.GIVEN)
        warning("Ignoring specified root state")
    x <- pars$x
    dx <- x[2] - x[1]
    if (root == ROOT.FLAT) {
        print("root d's =")
        print(paste(d.root, collapse=", "))
        p <- 1/((pars$nx - 1) * dx)
        print("root p =")
        print(paste(p, collapse=", "))
    }
    else if (root == ROOT.OBS) {
        print("getting prior at root")
        print(paste0("denominator = ", (sum(d.root) * dx)))
        print(paste0("dx = ", dx))
        p <- d.root/(sum(d.root) * dx)
        print(paste0("root d's (size=", length(d.root), ") ="))
        print(paste(d.root, collapse=", "))
        print("root p =")
        print(paste(p, collapse=", "))
    }
    else if (root == ROOT.GIVEN) {
        p <- root.f(x)
    }
    else {
        stop("Unsupported root mode")
    }
    p
}

make.rootfunc.quasse.debug <- function(cache) {
    root.idx <- cache$root
    nx <- cache$control$nx
    dx <- cache$control$dx
    function(res, pars, condition.surv, root, root.f, intermediates) {
        ## vals and lq come from all.branches.list
        vals <- matrix(res$vals, nx, 2)[seq_len(pars$lo$ndat),] # just making matrix
        lq <- res$lq # all normalizing factors for all nodes in tree
        d.root <- vals[, 2] # normalized D's at root
        root.p <- root.p.quasse.debug(d.root, pars$lo, root, root.f) # prior at root for each trait value bin (several options, user-defined); the default is the weighted sum, which is (d.root/sum(d.root) ==> returns a vector)

        if (condition.surv) {
            lambda <- pars$lo$lambda
            e.root <- vals[, 1]

            print("root e's = ")
            print(paste(e.root, collapse=", "))

            ## param = lambda, mu
            ## data D = tree (fixed)
            ## tip states are boundary conditions
            ## P(D|param) = L(param;D)
            ## at this point, d.root = P(D,tree_exists|params) = P(D|params,tree_exists) * P(tree_exists)
            ## but what we want is P(D|params,tree_exists)
            ## so we need to d.root / P(tree_exists)

            print(paste0("lambda (size=", length(lambda), ") = "))
            print(paste(lambda, collapse=", "))
            d.root <- d.root /
                sum(root.p * lambda * (1 - e.root)^2) * dx ## normalizing by this sum here accounts for conditioning for survival, P(tree_exists)
            ## e.root is the probability of extinction of a lineage starting at the root
            ## (1 - e.root) is the probability that lineage survives until today
            ## we square it because there are two lineages coming from the root (i.e., the definition of a tree)
            ## multiply by lambda because we need a sp'n event at the root generating those two lineages
            ## multiply by root.p to account for the trait value states (those influence the extinction probability)
            ## * dx: multiplying by bin size converts densities (median of bin is a density) into probabilities

            print(paste0("Denominator sum for conditioning on survival = ", sum(root.p * lambda * (1 - e.root)^2)))
            print(paste0("dx = ", dx))
            print(paste0("D's after conditioning on survival (size=", length(d.root), ") = "))
            print(paste(d.root, collapse=", "))
        }

        print(paste0("sum(lq) = ", sum(lq)))

        loglik <- log(sum(root.p * d.root) * dx) + sum(lq) # goes back to actual likelihood ("unnormalizing it")

        ## doesn't seem important (?)
        if (intermediates) {
            attr(loglik, "intermediates") <- res
            attr(loglik, "vals") <- vals
        }

        ## returns final result
        loglik
    }
}


##
##
##

control.C.1 <- list(tc=0.005,
                    dt.max=0.005,
                    nx=32,
                    dx=0.01,
                    r=4L,
                    xmid=0,
                    w=10,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC") # single dt covering entire bifurcating tree height

lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars))

## (ll.C.1 <- lik.C.1(pars, root=ROOT.FLAT))
## paste(rep("3.2258064516129", 25), collapse=", ")

## the above command will print a lot of things, some of those things are the unit tests expectations
## at the top, $sp1 and $sp2 correspond to the expectedEsDsAtTips
##
## sp1, elements 129-138, expected = 0.5043644, 0.5665408, 0.6347930, 0.7094919, 0.7910008, 0.8796719, 0.9758404, 1.0798193, 1.1918941, 1.3123163

## in R
## control.R.1 <- list(dt.max=0.01, method="fftR")
## lik.R.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
## (ll.R.1 <- lik.R.1(pars)) # -28.43912

## diversitree:::make.quasse

## testPruneTree32Bins

## esDsHiAtNodeInitial0[1] = 0.308986942687903, 0.350566009871371, 0.396747087835907, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.3217049979751, 5.55969772261993, 5.79383105522966, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768, 6.66449205783599, 6.45944719335828, 6.24507866733522, 6.02274864309609, 5.79383105522965, 5.55969772261993, 5.32170499797509, 5.08118112938378, 4.83941449038287, 4.59764281368466, 4.35704354065101, 4.11872537439949, 3.88372109966426, 3.65298170778044, 3.42737184095615, 3.20766654683839, 2.9945493127149, 2.7886113289072, 2.59035191331783, 2.40018001393971, 2.21841669358911, 2.04529849127956, 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583655, 0.879671919608544, 0.791000831787404, 0.709491856924628, 0.634793036713349, 0.566540754832024, 0.504364398303888, 0.447890605896858, 0.396747087835907, 0.350566009871371, 0.308986942687903

## esDsHiAtNodeInitial1[1] = 0.000254946647636669, 0.000319674822138109, 0.000399835934138456, 0.000498849425801072, 0.000620828141157003, 0.000770703934841743, 0.000954372730824099, 0.0011788613551308, 0.00145251860604505, 0.00178523314354266, 0.00218868086879601, 0.00267660451529771, 0.00326512817532484, 0.00397310942785545, 0.00482253160451986, 0.0058389385158292, 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755, 0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349, 0.0344513787810736, 0.0405809611459954, 0.0476817640292969, 0.0558851682975889, 0.0653363811239983, 0.0761952419644361, 0.08863696823876, 0.102852818461079, 0.119050648395517, 0.137455333812279, 0.158309031659599, 0.181871250031821, 0.208418696288452, 0.238244872152104, 0.271659384673712, 0.308986942687903, 0.350566009871371, 0.396747087835906, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.32170499797509, 5.55969772261993, 5.79383105522965, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768

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
