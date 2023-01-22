source("helper-diversitree.R")
library(jsonlite)
library(phytools)

context("MoSSE")

test_that("mosse", {
  
  ##1. test calculation of Q matrix
  ##2. test pde solver on a single branch
  ##3. test mosse on a whole tree
  ##4. fftC should give the same answer to fftR
  
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
  
  Q <- t(matrix(c(-0.5,0.2,0.1,0.1,
                  0.1,-0.8,0.3,0.4,
                  0.2,0.2,-0.6,0.2,
                  0.3,0.4,0.2,-0.7),4,4))
  
  
  ##3. test mosse on a whole tree
  #load a tree 
  load("phy.Rdata")
  types <- c(1, 1, 1, 3, 4, 1, 2, 1, 2, 3, 1, 3, 3, 4, 1)
    names(types) <- phy$tip.label
    states <- c(0.010707127, 0.009700834, 0.010450199, 0.009184059, 0.010065188, 0.010411688, 0.010025265, 0.010276946, 0.009451115, 0.010506482, 0.009989017, 0.009375481, 0.008813563, 0.009812874, 0.010246622)
  names(states) <- phy$tip.label
  states.sd <- rep(0.001,length(phy$tip.label))
  names(states.sd) <- phy$tip.label
  
  pars <- c(0,0.1,0.1,0.01,0.01,0,0.001)
  
  
  
  control.C.1 <- list(xmax=1)
  control.R.1 <- list(xmax=1, method="fftR")
  
  print(phy)
  
  writeNexus(phy, "test.tree")
  
  print(states)
  
  print(states.sd)
  
  print(sigmoid.x)
  
  print(constant.x)
  
  # lik.C.1 <- make.mosse(phy, types,ntypes=4,states,states.sd,sigmoid.x, constant.x, Q,
  #                     control.C.1)
  lik.R.1 <- make.mosse(phy, types,ntypes=4,states,states.sd,sigmoid.x, constant.x, Q,
                       control.R.1)
  # 
  # 
  # 
  # expect_that(round(lik.C.1(pars),4), equals(-307.0883)) 
  expect_that(round(lik.R.1(pars),4), equals(-307.0883)) 
  # 
  # ## Now, test root treatment:
  # expect_that(round(lik.C.1(pars, root=ROOT.FLAT),4), equals(-314.0188)) 
  # expect_that(round(lik.C.1(pars, condition.surv=FALSE),4), equals(-310.2858))
  # root.f <- function(x)
  #   dnorm(x, mean(states), sd(states))
  # expect_that(round(lik.C.1(pars, root=ROOT.GIVEN, root.f=root.f),4),
  #             equals(-379.7512)) 
  # 
  # ## With drift:
  # pars2 <- pars
  # pars2[6] <- 0.1
  # expect_that(round(lik.C.1(pars2),4), equals(-234.6243))
  # 
})
