## author: Fabio K. Mendes
##
## This R script gives us expected values for the JUnit tests
## inside QuaSSEFunctionsTest.java
##
## (1) testLogistic
## (2) testConstant
## (3) testQu2MacroevolFailLogistic
## (4) testQu2MacroevolFailConstant

## Auxiliary functions
## r: logistic growth rate
## xmid (or x0): sigmoid midpoint
## y0: curve y base value
## y1: curve maximum
sigmoid.x <- function (x, y0, y1, xmid, r) {
    curve.max.minus.base.value = y1 - y0
    to.add = curve.max.minus.base.value / (1 + exp(r * (xmid - x)))
    y0 + to.add # add base value back
}

## Test expectations
##
## (1) testLogistic

## The following code is inside diversitree's quasse.extent
drift <- 0.0
diffusion <- 0.01
nx <- 1024
dx <- 0.01
dt <- 0.05
y0 <- 0.1
y1 <- 0.2
xmid <- 0
r <- 2.5
death <- 0.03
hi.lo.ratio <- 4 # r=4L when making control.fft
w <- 5

mean4test <- drift * dt
sd4test <- sqrt(diffusion * dt)
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright) # number of useful bins
ndat.lo <- ndat[2]; # low res
ndat.hi <- ndat[1]; # high res

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01; xmin.hi < xmin.lo
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

## (2) testConstant (this is a trivial test, and requires no expectation)

## (3) testQu2MacroevolFailLogistic (this is a trivial test, and requires no expectation)

## (4) testQu2MacroevolFailConstant (this is a trivial test, and requires no expectation)

## NOTE: Tests 3 and 4 check if an Exception is raised (when the number of x and y bins are different)
