source("helper-diversitree.R")

context("MoSSE")

test_that("mosse", {

  ## requires R packages diversitree, mosse and C library fftw3 to be installed
  ## 1. test calculation of Q matrix
  ## 2. test pde solver on a single branch
  ## 3. fftC should give the same answer to fftR


  ## Imports that are generally hidden.
  make.pars.mosse <- diversitree:::make.pars.mosse
  make.pde.mosse.fftC <- diversitree:::make.pde.mosse.fftC
  make.pde.mosse.fftR <- diversitree:::make.pde.mosse.fftR

  make.branches.mosse.fftC <- diversitree:::make.branches.mosse.fftC
  make.branches.mosse.fftR <- diversitree:::make.branches.mosse.fftR

  ## Basic control list.
  control.fft <- list(tc=1.3,
                      dt.max=0.01,
                      nx=1024,
                      dx=10^-4,
                      r=4L,
                      w=5,
                      method="fftR",
                      flags=0L, # fftC only
                      atol=1e-6, # mol only
                      rtol=1e-6, # mol only
                      eps=1e-6,  # perhaps scale with dx?
                      verbose=0L,
                      ntypes=4L)

  lambda <- sigmoid.x
  mu <- constant.x
  diffusion <- 0.001
  sd <- 0.001

  ## test pde solver on a single branch
  len <- 2 # Integrate down a branch length of 2

  Q <- t(matrix(c(-0.5,0.2,0.1,0.1,
                  0.1,-0.8,0.3,0.4,
                  0.2,0.2,-0.6,0.2,
                  0.3,0.4,0.2,-0.7),4,4))

  # test case 1
  drift <- 0.0
  args <- list(lambda=1:4, mu=5, drift=6, diffusion=7)
  pars <- c(0, 0.1, 0.1, 0.01, 0.01, drift, diffusion)

  control.fft$method <- "fftR"
  cache <- list(args=args,control=control.fft,lambda=lambda,mu=mu,Q=Q)
  f.pars <- make.pars.mosse(cache)
  pars.fft <- f.pars(pars)

  ## initial conditions:
  vars.fft <- matrix(0, control.fft$nx, 5)
  vars.fft[seq_len(pars.fft$lo$ndat),2] <- dnorm(pars.fft$lo$x, 0.1, sd)
  # vars fft is data array passed into fft

  vars.hi.fft <- matrix(0, control.fft$nx*control.fft$r, 5)
  vars.hi.fft[seq_len(pars.fft$hi$ndat),2] <-
    dnorm(pars.fft$hi$x, 0.1, sd)

  ## run fftR
  pde.fftR <- with(control.fft, make.pde.mosse.fftR(nx, dx*r, dt.max, 5L))
  ans.fftR <- pde.fftR(vars.fft, len, pars.fft$lo, 0, control.fft$dt.max)

  ## run fftC
  control.fft$method <- "fftC"
  cache <- list(args=args,control=control.fft,lambda=lambda,mu=mu,Q=Q)
  f.pars <- make.pars.mosse(cache) # set mosse parameters
  pars.fft <- f.pars(pars)

  ## TEST
  # every 4th Q.hi matrix (or every 16 elements) should equal every Q.lo matrix
  expect_that(pars.fft$lo$Q[seq(1,pars.fft$lo$ndat*16,16)], equals(pars.fft$hi$Q[seq((control.fft$r-1)*16+1,pars.fft$hi$ndat*16,control.fft$r*16)]))

  ## Bail here if no FFTW support
  if (!check.fftC(FALSE)) {
    next
  }
  pde.fftC <- with(control.fft, make.pde.mosse.fftC(nx, dx*r, dt.max, 5L, flags))
  ans.fftC <- pde.fftC(vars.fft, len, pars.fft$lo, 0, control.fft$dt.max)

  ## TEST
  # FFTC method should give the same answer to the FFTR method
  print("ans.fftR")
  print(ans.fftR[[1]])
  print("ans.fftC")
  print(ans.fftC[[1]])
  # C ans array
  # print(ans.fftC[[2]])

  # test case 1
  expect_that(ans.fftC, equals(ans.fftR))

})
