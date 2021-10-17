# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# PropagatesQuaSSE
# (1) testPropagateTimeOneChQuaSSETest
# (2) testMakeNormalKernInPlaceAndFFT
# (3) testConvolve
# (4) testPropagateChOneCh48QuaSSETest
# (5) testLogistic
# (6) testDimensions
# (7) testInitializationOfTips
# (8) testIntegrateOneBranchLoRes48BinsOutsideClassJustT
# (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX
# (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX and testPropagateChOneCh1024QuaSSETest
# (11) testIntegrateOneBranchLoRes4096BinsOutsideClassJustX and testPropagateChOneCh4096QuaSSETest (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)

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
        print(paste("nx =", nx))
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
paste(res[,2], collapse=", ") # 0.0011788613551308, 0.00267660451529771, 0.0058389385158292, 0.0122380386022755, 0.0247150623561646, 0.0478008882241639, 0.088827888339617, 0.158599420774613, 0.272077439492024, 0.448458176810809, 0.710214708266688, 1.0806760140649, 1.57993546382425, 2.21932565164128, 2.99530083804266, 3.88416335936269, 4.83940600155166, 5.79327421787607, 6.66336349661661, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119626, 6.66336349661661, 5.79327421787606, 4.83940600155166, 3.8841633593627, 2.99530083804265, 2.21932565164129, 1.57993546382425, 1.0806760140649, 0.710214708266687, 0.44845817681081, 0.272077439492023, 0.158309031659599, 0.0886369682387602, 0.0476817640292969, 0.0246443833694604, 0, 0, 0, 0, 0, 0, 0, 0, 0



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
dx <- 0.0005; diffusion <- 0.001; w <- 10

sd4test <- sqrt(diffusion * dt)
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2]; ndat.lo
ndat.hi <- ndat[1];

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

paste(x.lo[1:10], collapse=", ") # expectedXLoFirst10 = -4.99, -4.98, -4.97, -4.96, -4.95, -4.94, -4.93, -4.92, -4.91, -4.9
paste(x.lo[990:999], collapse=", ") # expectedXLoLast10 = 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99

paste(x.hi[1:10], collapse=", ") # expectedXHiFirst10 = -4.9975, -4.995, -4.9925, -4.99, -4.9875, -4.985, -4.9825, -4.98, -4.9775, -4.975
paste(x.hi[3990:3999], collapse=", ") # expectedXHiFirst10 = 4.975, 4.9775, 4.98, 4.9825, 4.985, 4.9875, 4.99, 4.9925, 4.995, 4.9975-4.9975, -4.995, -4.9925, -4.99, -4.9875, -4.985, -4.9825, -4.98, -4.9775, -4.975

lambda <- sigmoid.x
mu <- constant.x

# lambdas
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaLoFirt10 = 0.100000382097925, 0.100000391770742, 0.100000401688425, 0.100000411857174, 0.100000422283344, 0.100000432973452, 0.100000443934178, 0.100000455172374, 0.100000466695064, 0.10000047850945
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[989:999], collapse=", ") # expectedLambdaLoLast10 = 0.199999509377086, 0.199999521490551, 0.199999533304936, 0.199999544827626, 0.199999556065822, 0.199999567026548, 0.199999577716656, 0.199999588142826, 0.199999598311575, 0.199999608229258, 0.199999617902075

paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.100000375000363, 0.100000377351446, 0.100000379717269, 0.100000382097925, 0.100000384493506, 0.100000386904106, 0.10000038932982, 0.100000391770742, 0.100000394226967, 0.100000396698591
paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[3989:3999], collapse=", ") # expectedLambdaHiLast10 = 0.199999600814288, 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637

# mus
paste(do.call(mu, c(list(x.lo), death))[1:10], collapse=", ") # expectedMuLoFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.lo), death))[989:999], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03

paste(pars.fft$lo$mu[1:10], collapse=", ") # expectedMuLoFirt10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$lo$mu[989:999], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$hi$mu[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(pars.fft$hi$mu[3989:3999], collapse=", ") # expectedLambdaHiLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03



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

ext.fft.just.t <- quasse.extent(control.fft.48, drift, diffusion) # prepares X axis stuff
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



# (8.5)

pde.fftR.t.x <- with(control.fft.48, diversitree:::make.pde.quasse.fftR(control.fft.48$nx, control.fft.48$dx, control.fft.48$dt.max, 2L))
ans.fftR.t.x <- pde.fftR.t.x(vars.fft.just.t, control.fft.48$dt.max, pars.fft.just.t$lo, 0)


# (9) testIntegrateOneBranchLoRes48BinsOutsideClassJustX

args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x <- quasse.extent(control.fft.48, drift, diffusion) # prepares X axis stuff
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



# (10) testPropagateChOneCh1024QuaSSETest and testIntegrateOneBranchLoRes1024BinsOutsideClassJustX

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
paste(vars.fft.just.x.1024[1:48,2], collapse=", ") # 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

kern.just.x.1024 <- fftR.make.kern(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$diffusion), control.fft.just.x.1024$nx, control.fft.just.x.1024$dx, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])

fy.just.x.1024 <- fft(kern.just.x.1024)
paste(Re(fy.just.x.1024)[1:48], collapse=", ") # 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x.1024, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x.1024
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints (in Dell XPS + Ubuntu, this differs!)

ds.prop.x.1024 <- fftR.propagate.x(vars.fft.just.x.1024, control.fft.just.x.1024$nx, fy.just.x.1024, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])
paste(ds.prop.x.1024[1:48,2], collapse=", ") # 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

paste(ds.prop.x.1024[848:895,2], collapse=", ") # 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948

pde.fftC <- with(control.fft.just.x.1024$nx, diversitree:::make.pde.mosse.fftC(control.fft.just.x.1024$nx, control.fft.just.x.1024$nx*control.fft.just.x.1024$r, control.fft.just.x.1024$dt.max, 4L, flags))
ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0, control.fft$dt.max)



# (11) testPropagateChOneCh4096QuaSSETest and testIntegrateOneBranchLoRes4096BinsOutsideClassJustX (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)

# initialization
vars.fft.just.x.4096 <- matrix(0, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, 2) # low resolution
vars.fft.just.x.4096[seq_len(ext.fft.just.x.1024$ndat[1]),2] <- dnorm(ext.fft.just.x.1024$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.4096$lo$padding should be 0.0
paste(vars.fft.just.x.4096[1:48,2], collapse=", ") # 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846

kern.just.x.4096 <- fftR.make.kern(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$diffusion), control.fft.just.x.1024$nx * control.fft.just.x.1024$r, control.fft.just.x.1024$dx, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])

fy.just.x.4096 <- fft(kern.just.x.4096)
paste(Re(fy.just.x.4096)[1:48], collapse=", ") # 1, 0.999952939166244, 0.999811769952891, 0.999576532217434, 0.999247292368191, 0.998824143333056, 0.998307204515792, 0.997696621739869, 0.996992567179915, 0.996195239280797, 0.995304862664403, 0.994321688024177, 0.993245992007481, 0.992078077085863, 0.990818271413312, 0.989466928672585, 0.988024427909734, 0.986491173356906, 0.984867594243561, 0.983154144596211, 0.981351303026827, 0.979459572510042, 0.977479480149295, 0.975411576932071, 0.973256437474405, 0.971014659754794, 0.968686864837709, 0.966273696586883, 0.963775821368548, 0.961193927744829, 0.958528726157476, 0.955780948602156, 0.952951348293497, 0.950040699321108, 0.947049796296786, 0.943979453993153, 0.940830506973933, 0.937603809216109, 0.934300233724213, 0.930920672136968, 0.92746603432656, 0.92393724799076, 0.920335258238176, 0.916661027166888, 0.912915533436717, 0.909099771835414, 0.905214752839018, 0.90126150216667

ds.prop.x.4096 <- fftR.propagate.x(vars.fft.just.x.4096, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, fy.just.x.4096, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])
paste(ds.prop.x.4096[1:48,2], collapse=", ") # 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846
paste(ds.prop.x.4096[1001:1048,2], collapse=", ") # 1.12964601442883, 1.13523957985681, 1.14085371385815, 1.14648844781881, 1.1521438128928, 1.15781984000018, 1.16351655982519, 1.16923400281429, 1.17497219917421, 1.18073117887008, 1.18651097162341, 1.19231160691021, 1.19813311395904, 1.20397552174904, 1.20983885900802, 1.21572315421049, 1.22162843557575, 1.22755473106588, 1.23350206838385, 1.23947047497152, 1.24545997800775, 1.25147060440637, 1.25750238081428, 1.26355533360948, 1.26962948889911, 1.27572487251749, 1.28184151002416, 1.28797942670193, 1.29413864755491, 1.30031919730656, 1.30652110039772, 1.31274438098464, 1.31898906293703, 1.32525516983611, 1.33154272497262, 1.33785175134486, 1.34418227165675, 1.35053430831583, 1.35690788343134, 1.36330301881221, 1.36971973596514, 1.37615805609261, 1.38261800009092, 1.38909958854822, 1.39560284174258, 1.40212777963999, 1.40867442189243, 1.41524278783588

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
