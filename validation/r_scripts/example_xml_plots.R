library(ggplot2)
library(gridExtra)
library(RColorBrewer)

lighten <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col*factor
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
    col
}

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

make.post.plot <- function(a.df, var.col, param.name, x.min.max.vec, a.color.fill, a.color.line, a.truth, a.mle, a.mean) {
    post.samples.plot = ggplot(a.df, aes(x=a.df[,var.col])) +
        geom_histogram(aes(y=stat(density)), color=a.color.line, fill=a.color.fill, alpha=.01, bins=40, size=2) +
        geom_histogram(aes(y=stat(density)), fill=a.color.fill, bins=40) +
        annotate("point", y=0, x=a.truth, size=3) +
        geom_vline(xintercept=a.mle, color="red") +
        geom_vline(xintercept=a.mean) +
        xlab(param.name) + ylab("Density") + 
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            plot.title = element_text(hjust=0.5),
            axis.line = element_line(),
            axis.ticks = element_line(color="black"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12)
        ) 
    return(post.samples.plot)
}


# ----- Plotting BiSSE_fixed_tree_SDSEP and BiSSE_fixed_tree_HSDSEP ----- #
pal <- colorRampPalette(pal)(7)

log.df <- read.table("../examples/BiSSE_fixed_tree_HSDSEP.log", header=TRUE)
log.df <- log.df[12:nrow(log.df),] # removing burnin and extra generations

## -- MLE from diversitree -- ##
##      lambda0      lambda1          mu0          mu1          q01          q10 
## 2.304392e-08 3.622767e-01 4.535375e-01 6.484517e-02 2.127918e-01 1.526163e-01

l0.truth <- 0.15
l0.mle <- 2.304392e-08
l0.mean <- mean(log.df$"Lambda1")
l1.truth <- 0.3
l1.mle <- 3.622767e-01
l1.mean <- mean(log.df$"Lambda2")
m0.truth <- 0.1
m0.mle <- 4.535375e-01
m0.mean <- mean(log.df$"Mu1")
m1.truth <- 0.1
m1.mle <- 6.484517e-02
m1.mean <- mean(log.df$"Mu2")
q01.truth <- 0.1
q01.mle <- 2.127918e-01
q01.mean <- mean(log.df$"FlatQMatrix1")
q10.truth <- 0.1
q10.mle <- 1.526163e-01
q10.mean <- mean(log.df$"FlatQMatrix2")

l0.min.max <- range(log.df$Lambda1)
fill.color = lighten(pal[1], 1.3)
l0.plot <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, fill.color, pal[1], l0.truth, l0.mle, l0.mean)
## l0.plot

l1.min.max <- range(log.df$Lambda2)
fill.color = lighten(pal[2], 1.3)
l1.plot <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, fill.color, pal[2], l1.truth, l1.mle, l1.mean)
## l1.plot

m0.min.max <- range(log.df$Mu1)
fill.color = lighten(pal[3], 1.3)
m0.plot <- make.post.plot(log.df, 7, expression(mu[0]), m0.min.max, fill.color, pal[3], m0.truth, m0.mle, m0.mean)
## m0.plot

m1.min.max <- range(log.df$Mu2)
fill.color = lighten(pal[4], 1.3)
m1.plot <- make.post.plot(log.df, 8, expression(mu[1]), m1.min.max, fill.color, pal[4], m1.truth, m1.mle, m1.mean)
## m1.plot

q01.min.max <- range(log.df$FlatQMatrix1)
fill.color = lighten(pal[5], 4.5)
q01.plot <- make.post.plot(log.df, 9, expression(q[0][1]), q01.min.max, fill.color, pal[5], q01.truth, q01.mle, q01.mean)
## q01.plot

q10.min.max <- range(log.df$FlatQMatrix2)
fill.color = lighten(pal[6], 2.2)
q10.plot <- make.post.plot(log.df, 10, expression(q[1][0]), q10.min.max, fill.color, pal[6], q10.truth, q10.mle, q10.mean)
## q10.plot

pdf("posts_BiSSE_HDSEP_60spp.pdf", width=6, height=6)
grid.arrange(l0.plot, l1.plot,
             m0.plot, m1.plot,
             q01.plot, q10.plot)
dev.off()

# ----- Plotting BiSSE_fixed_tree_on_HiSSE_HSDSEP  ----- #

## pal <- brewer.pal(10, "Accent")
pal <- colorRampPalette(pal)(11)

log.df <- read.table("../examples/BiSSE_fixed_tree_on_HiSSE_HSDSEP.log", header=TRUE)
log.df <- log.df[22:nrow(log.df),] # removing burnin and extra generations

## -- MLE from hisse -- ##
##                       lambda         mu
## rate0A             0.1594821 0.08712492
## rate1A             0.09142577 0.07480436
## rate0B                  0  0
## rate1B             0.2472443 0.1383931

## Transition Rates
##           (0A)       (1A) (0B)      (1B)
## (0A)        NA 0.29523936    0 0.0000000
## (1A) 0.4509993         NA    0 0.1355034
## (0B) 0.0000000 0.00000000   NA 0.0000000
## (1B) 0.0000000 0.01983649    0        NA

l0.truth <- 0.1
l0.mle <- 0.1594821
l0.mean <- mean(log.df$"Lambda1")
l1.truth <- 0.1
l1.mle <- 0.09142577
l1.mean <- mean(log.df$"Lambda2")
l2.truth <- 0.3
l2.mle <- 0.2472443
l2.mean <- mean(log.df$"Lambda3")
m0.truth <- 0.05
m0.mle <- 0.08712492
m0.mean <- mean(log.df$"Mu1")
m1.truth <- 0.05
m1.mle <- 0.07480436
m1.mean <- mean(log.df$"Mu2")
m2.truth <- 0.05
m2.mle <- 0.1383931
m2.mean <- mean(log.df$"Mu3")
q01.truth <- 0.1
q01.mle <- 0.29523936
q01.mean <- mean(log.df$"FlatQMatrix1")
q10.truth <- 0.1
q10.mle <- 0.4509993
q10.mean <- mean(log.df$"FlatQMatrix2")
q12.truth <- 0.1
q12.mle <- 0.1355034
q12.mean <- mean(log.df$"FlatQMatrix3")
q21.truth <- 0.1
q21.mle <- 0.01983649
q21.mean <- mean(log.df$"FlatQMatrix4")

l0.min.max <- range(log.df$Lambda1)
fill.color = lighten(pal[1], 1.3)
l0.plot <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, fill.color, pal[1], l0.truth, l0.mle, l0.mean)
## l0.plot

l1.min.max <- range(log.df$Lambda2)
fill.color = lighten(pal[2], 1.3)
l1.plot <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, fill.color, pal[2], l1.truth, l1.mle, l1.mean)
## l1.plot

l2.min.max <- range(log.df$Lambda3)
fill.color = lighten(pal[3], 1.3)
l2.plot <- make.post.plot(log.df, 7, expression(lambda[2]), l2.min.max, fill.color, pal[3], l2.truth, l2.mle, l2.mean)
## l2.plot

m0.min.max <- range(log.df$Mu1)
fill.color = lighten(pal[4], 1.3)
m0.plot <- make.post.plot(log.df, 8, expression(mu[0]), m0.min.max, fill.color, pal[4], m0.truth, m0.mle, m0.mean)
## m0.plot

m1.min.max <- range(log.df$Mu2)
fill.color = lighten(pal[5], 1.4)
m1.plot <- make.post.plot(log.df, 9, expression(mu[1]), m1.min.max, fill.color, pal[5], m1.truth, m1.mle, m1.mean)
## m1.plot 

m2.min.max <- range(log.df$Mu3)
fill.color = lighten(pal[6], 1.3)
m2.plot <- make.post.plot(log.df, 10, expression(mu[2]), m2.min.max, fill.color, pal[6], m2.truth, m2.mle, m2.mean)
## m2.plot

q01.min.max <- range(log.df$FlatQMatrix1)
fill.color = lighten(pal[7], 2.0)
q01.plot <- make.post.plot(log.df, 11, expression(q[0][1]), q01.min.max, fill.color, pal[7], q01.truth, q01.mle, q01.mean)
## q01.plot

q10.min.max <- range(log.df$FlatQMatrix2)
fill.color = lighten(pal[8], 8.0)
q10.plot <- make.post.plot(log.df, 12, expression(q[1][0]), q10.min.max, fill.color, pal[8], q10.truth, q10.mle, q10.mean)
## q10.plot

q12.min.max <- range(log.df$FlatQMatrix3)
fill.color = lighten(pal[9], 2.0)
q12.plot <- make.post.plot(log.df, 13, expression(q[1][2]), q12.min.max, fill.color, pal[9], q12.truth, q12.mle, q12.mean)
## q12.plot

q21.min.max <- range(log.df$FlatQMatrix4)
fill.color = lighten(pal[10], 1.4)
q21.plot <- make.post.plot(log.df, 14, expression(q[2][1]), q21.min.max, fill.color, pal[10], q21.truth, q21.mle, q21.mean)
## q21.plot

pdf("posts_BiSSE_on_HiSSE_HDSEP_60spp.pdf", width=9, height=9)
grid.arrange(l0.plot, l1.plot, l2.plot,
             m0.plot, m1.plot, m2.plot,
             q01.plot, q10.plot, q12.plot, q21.plot)
dev.off()
