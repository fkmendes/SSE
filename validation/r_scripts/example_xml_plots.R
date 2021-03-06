library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(diversitree)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

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
        annotate("point", y=0, x=a.mean, size=3) +
        annotate("point", y=0, x=a.mle, size=3, color="red") +
        geom_vline(xintercept=a.truth, lty="dashed") +
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

make.post.plot.no.mle <- function(a.df, var.col, param.name, x.min.max.vec, a.color.fill, a.color.line, a.truth, a.mean) {
    post.samples.plot = ggplot(a.df, aes(x=a.df[,var.col])) +
        geom_histogram(aes(y=stat(density)), color=a.color.line, fill=a.color.fill, alpha=.01, bins=40, size=2) +
        geom_histogram(aes(y=stat(density)), fill=a.color.fill, bins=40) +
        annotate("point", y=0, x=a.mean, size=3) +
        geom_vline(xintercept=a.truth, lty="dashed") +
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

make.post.dens.plot <- function(a.df, param.name, param.xlab, some.colors, some.means, a.truth) {
    the.plot = ggplot(a.df, aes(x=a.df[,param.name], fill=Model)) + geom_density(alpha=.25, color=NA) +
        xlab(param.xlab) + ylab("Density") +
        annotate("point", y=0, x=some.means[1], color=some.colors[1], size=3) +
        annotate("point", y=0, x=some.means[2], color=some.colors[2], size=3) +
        annotate("point", y=0, x=some.means[3], color=some.colors[3], size=3) +
        annotate("point", y=0, x=some.means[4], color=some.colors[4], size=3) +
        annotate("point", y=0, x=some.means[5], color=some.colors[5], size=3) +
        geom_vline(xintercept=a.truth, lty="dashed") +

    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(),
        axis.line = element_line(),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)
    ) +
    scale_fill_manual(values=some.colors[1:5])

    return(the.plot)
}
    
# ----- Plotting BiSSE_fixed_tree_SDSEP ----- #
pal <- brewer.pal(8, "Accent")
pal <- colorRampPalette(pal)(14)

log.df <- read.table("../examples/BiSSE_fixed_tree_SDSEP.log", header=TRUE)
log.df <- log.df[102:nrow(log.df),] # removing burnin and extra generations

## -- MLE from diversitree -- ##
##      lambda0      lambda1          mu0          mu1          q01          q10 
## 9.148731e-02 5.415631e-01 2.056091e-05 3.561561e-05 2.588266e-01 8.936368e-02 

l0.truth <- 0.1
l0.mle <- 9.148731e-02
l0.mean <- mean(log.df$"Lambda1")
l1.truth <- 0.5
l1.mle <- 5.415631e-01
l1.mean <- mean(log.df$"Lambda2")
m0.truth <- 0.05
m0.mle <- 2.056091e-05
m0.mean <- mean(log.df$"Mu1")
m1.truth <- 0.05
m1.mle <- 3.561561e-05
m1.mean <- mean(log.df$"Mu2")
q01.truth <- 0.1
q01.mle <- 2.588266e-01
q01.mean <- mean(log.df$"FlatQMatrix1")
q10.truth <- 0.1
q10.mle <- 8.936368e-02
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
fill.color = lighten(pal[4], 1.1)
m1.plot <- make.post.plot(log.df, 8, expression(mu[1]), m1.min.max, fill.color, pal[4], m1.truth, m1.mle, m1.mean)
## m1.plot

q01.min.max <- range(log.df$FlatQMatrix1)
fill.color = lighten(pal[5], 1.3)
q01.plot <- make.post.plot(log.df, 9, expression(q[0][1]), q01.min.max, fill.color, pal[5], q01.truth, q01.mle, q01.mean)
## q01.plot

q10.min.max <- range(log.df$FlatQMatrix2)
fill.color = lighten(pal[6], 1.3)
q10.plot <- make.post.plot(log.df, 10, expression(q[1][0]), q10.min.max, fill.color, pal[6], q10.truth, q10.mle, q10.mean)
## q10.plot

pdf("plots/BiSSE_fixed_tree_SDSEP_posteriors_120spp.pdf", width=6, height=5)
grid.arrange(l0.plot, l1.plot,
             m0.plot, m1.plot,
             q01.plot, q10.plot)
dev.off()

# ----- Plotting ClaSSE_fixed_tree_SDSEP ----- #
log.df <- read.table("../examples/ClaSSE_fixed_tree_SDSEP.log", header=TRUE)
log.df <- log.df[102:nrow(log.df),] # removing burnin and extra generations

l.s.truth <- 0.1
l.s.mean <- mean(log.df$SympatricRate)
l.s.min.max <- range(log.df$SympatricRate)

l.ss.truth <- 0.3
l.ss.mean <- mean(log.df$SubsympatricRate)
l.ss.min.max <- range(log.df$SubsympatricRate)

l.v.truth <- 0.5
l.v.mean <- mean(log.df$VicariantRate)
l.v.min.max <- range(log.df$VicariantRate)

m0.truth <- 0.05
m0.mean <- mean(log.df$Mu1)
m0.min.max <- range(log.df$Mu1)

m1.truth <- 0.05
m1.mean <- mean(log.df$Mu2)
m1.min.max <- range(log.df$Mu2)

m2.truth <- 0.05
m2.mean <- mean(log.df$Mu3)
m2.min.max <- range(log.df$Mu3)

q01.truth <- 0.1
q01.mean <- mean(log.df$FlatQMatrix1)
q01.min.max <- range(log.df$FlatQMatrix1)

q02.truth <- 0.1
q02.mean <- mean(log.df$FlatQMatrix2)
q02.min.max <- range(log.df$FlatQMatrix2)

q10.truth <- 0.1
q10.mean <- mean(log.df$FlatQMatrix3)
q10.min.max <- range(log.df$FlatQMatrix3)

q12.truth <- 0.1
q12.mean <- mean(log.df$FlatQMatrix4)
q12.min.max <- range(log.df$FlatQMatrix4)

q20.truth <- 0.1
q20.mean <- mean(log.df$FlatQMatrix5)
q20.min.max <- range(log.df$FlatQMatrix5)

q21.truth <- 0.1
q21.mean <- mean(log.df$FlatQMatrix6)
q21.min.max <- range(log.df$FlatQMatrix6)

fill.color = lighten(pal[1], 1.3)
l.s.plot <- make.post.plot.no.mle(log.df, 5, expression(lambda[sym]), l.s.min.max, fill.color, pal[1], l.s.truth, l.s.mean)
l.s.plot

fill.color = lighten(pal[2], 1.3)
l.ss.plot <- make.post.plot.no.mle(log.df, 6, expression(lambda[subsym]), l.ss.min.max, fill.color, pal[2], l.ss.truth, l.ss.mean)
l.ss.plot

fill.color = lighten(pal[3], 1.3)
l.v.plot <- make.post.plot.no.mle(log.df, 7, expression(lambda[vic]), l.v.min.max, fill.color, pal[3], l.v.truth, l.v.mean)
l.v.plot

fill.color = lighten(pal[4], 1.1)
m0.plot <- make.post.plot.no.mle(log.df, 8, expression(mu[0]), m0.min.max, fill.color, pal[4], m0.truth, m0.mean)
m0.plot

fill.color = lighten(pal[5], 1.1)
m1.plot <- make.post.plot.no.mle(log.df, 9, expression(mu[1]), m1.min.max, fill.color, pal[5], m1.truth, m1.mean)
m1.plot

fill.color = lighten(pal[6], 1.3)
m2.plot <- make.post.plot.no.mle(log.df, 10, expression(mu[3]), m2.min.max, fill.color, pal[6], m2.truth, m2.mean)
m2.plot

fill.color = lighten(pal[7], 1.3)
q01.plot <- make.post.plot.no.mle(log.df, 11, expression(q[0][1]), q01.min.max, fill.color, pal[7], q01.truth, q01.mean)
q01.plot

fill.color = lighten(pal[8], 1.3)
q02.plot <- make.post.plot.no.mle(log.df, 12, expression(q[0][2]), q02.min.max, fill.color, pal[8], q02.truth, q02.mean)
q02.plot

fill.color = lighten(pal[9], 2.0)
q10.plot <- make.post.plot.no.mle(log.df, 13, expression(q[1][0]), q10.min.max, fill.color, pal[9], q10.truth, q10.mean)
q10.plot

fill.color = lighten(pal[10], 1.6)
q12.plot <- make.post.plot.no.mle(log.df, 14, expression(q[1][2]), q12.min.max, fill.color, pal[10], q12.truth, q12.mean)
q12.plot

fill.color = lighten(pal[11], 2.2)
q20.plot <- make.post.plot.no.mle(log.df, 15, expression(q[2][0]), q20.min.max, fill.color, pal[11], q20.truth, q20.mean)
q20.plot

fill.color = lighten(pal[12], 1.3)
q21.plot <- make.post.plot.no.mle(log.df, 16, expression(q[2][1]), q21.min.max, fill.color, pal[12], q21.truth, q21.mean)
q21.plot

pdf("plots/ClaSSE_fixed_tree_SDSEP_posteriors_120spp.pdf", width=9, height=7.5)
grid.arrange(l.s.plot, l.ss.plot, l.v.plot,
             m0.plot, m1.plot, m2.plot,
             q01.plot, q02.plot, q10.plot,
             q12.plot, q20.plot, q21.plot)
dev.off()

# ----- Plotting BiSSE_fixed_tree_on_HiSSE_HSDSEP ----- #
log.df <- read.table("../examples/BiSSE_fixed_tree_on_HiSSE_HSDSEP.log", header=TRUE)
log.df <- log.df[102:nrow(log.df),] # removing burnin and extra generations

## -- MLE from diversitree -- ##
##      lambda0      lambda1          mu0          mu1          q01          q10 
## 5.490635e-02 5.094458e-01 1.254943e-06 1.949697e-01 4.718035e-02 5.817885e-03 

l0.truth <- 0.1
l0.mle <- 5.490635e-02
l0.mean <- mean(log.df$Lambda1)
l0.min.max <- range(log.df$Lambda1)

l1.truth <- 0.1
l1.mle <- 5.094458e-01
l1.mean <- mean(log.df$Lambda2)
l1.min.max <- range(log.df$Lambda2)

m0.truth <- 0.05
m0.mle <- 1.254943e-06
m0.mean <- mean(log.df$Mu1)
m0.min.max <- range(log.df$Mu1)

m1.truth <- 0.05
m1.mle <- 1.949697e-01
m1.mean <- mean(log.df$Mu2)
m1.min.max <- range(log.df$Mu2)

q01.truth <- 0.1
q01.mle <- 4.718035e-02
q01.mean <- mean(log.df$FlatQMatrix1)
q01.min.max <- range(log.df$FlatQMatrix1)

q10.truth <- 0.1
q10.mle <- 5.817885e-03
q10.mean <- mean(log.df$FlatQMatrix2)
q10.min.max <- range(log.df$FlatQMatrix2)

fill.color = lighten(pal[1], 1.3)
l0.plot <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, fill.color, pal[1], l0.truth, l0.mle, l0.mean)
l0.plot.bisse <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, "gray90", "gray60", l0.truth, l0.mle, l0.mean)
## l0.plot

fill.color = lighten(pal[2], 1.3)
l1.plot <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, fill.color, pal[2], l1.truth, l1.mle, l1.mean)
l1.plot.bisse <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, "lightskyblue", "steelblue2", l1.truth, l1.mle, l1.mean)
## l1.plot

fill.color = lighten(pal[3], 1.3)
m0.plot <- make.post.plot(log.df, 7, expression(mu[0]), m0.min.max, fill.color, pal[3], m0.truth, m0.mle, m0.mean)
## m0.plot

fill.color = lighten(pal[4], 1.1)
m1.plot <- make.post.plot(log.df, 8, expression(mu[1]), m1.min.max, fill.color, pal[4], m1.truth, m1.mle, m1.mean)
## m1.plot 

fill.color = lighten(pal[5], 1.3)
q01.plot <- make.post.plot(log.df, 9, expression(q[0][1]), q01.min.max, fill.color, pal[5], q01.truth, q01.mle, q01.mean)
## q01.plot

fill.color = lighten(pal[6], 1.3)
q10.plot <- make.post.plot(log.df, 10, expression(q[1][0]), q10.min.max, fill.color, pal[6], q10.truth, q10.mle, q10.mean)
## q10.plot

pdf("plots/BiSSE_fixed_tree_on_HiSSE_HSDSEP_posteriors_120spp.pdf", width=6, height=5)
grid.arrange(l0.plot, l1.plot,
             m0.plot, m1.plot,
             q01.plot, q10.plot)
dev.off()

# ----- Plotting HiSSE_fixed_tree_on_HiSSE_HSDSEP  ----- #
log.df <- read.table("../examples/HiSSE_fixed_tree_on_HiSSE_HSDSEP.log", header=TRUE)
log.df <- log.df[102:nrow(log.df),] # removing burnin and extra generations

## -- MLE from hisse -- ##

##                       lambda         mu
## rate0A             0.02913626 0.005376561
## rate1A             0.0636953 1.312858e-10
## rate1B             0.6093937 0.2787272

## Transition Rates
##            (0A)       (1A) (0B)       (1B)
## (0A)         NA 0.05197767    0 0.00000000
## (1A) 0.07369927         NA    0 0.03730514
## (0B) 0.00000000 0.00000000   NA 0.00000000
## (1B) 0.00000000 0.04861138    0         NA

l0.truth <- 0.1
l0.mle <- 0.02913626
l0.mean <- mean(log.df$Lambda1)
l0.min.max <- range(log.df$Lambda1)

l1.truth <- 0.1
l1.mle <- 0.0636953
l1.mean <- mean(log.df$Lambda2)
l1.min.max <- range(log.df$Lambda2)

l2.truth <- 0.5
l2.mle <- 0.6093937
l2.mean <- mean(log.df$Lambda3)
l2.min.max <- range(log.df$Lambda3)

m0.truth <- 0.05
m0.mle <- 0.005376561
m0.mean <- mean(log.df$Mu1)
m0.min.max <- range(log.df$Mu1)

m1.truth <- 0.05
m1.mle <- 1.312858e-10
m1.mean <- mean(log.df$Mu2)
m1.min.max <- range(log.df$Mu2)

m2.truth <- 0.05
m2.mle <- 0.2787272
m2.mean <- mean(log.df$Mu3)
m2.min.max <- range(log.df$Mu3)

q01.truth <- 0.1
q01.mle <- 0.05197767
q01.mean <- mean(log.df$FlatQMatrix1)
q01.min.max <- range(log.df$FlatQMatrix1)

q10.truth <- 0.1
q10.mle <- 0.07369927
q10.mean <- mean(log.df$FlatQMatrix2)
q10.min.max <- range(log.df$FlatQMatrix2)

q12.truth <- 0.1
q12.mle <- 0.03730514
q12.mean <- mean(log.df$FlatQMatrix3)
q12.min.max <- range(log.df$FlatQMatrix3)

q21.truth <- 0.1
q21.mle <- 0.04861138
q21.mean <- mean(log.df$FlatQMatrix4)
q21.min.max <- range(log.df$FlatQMatrix4)

fill.color = lighten(pal[1], 1.3)
l0.plot <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, fill.color, pal[1], l0.truth, l0.mle, l0.mean)
l0.plot.hisse <- make.post.plot(log.df, 5, expression(lambda[0]), l0.min.max, "gray90", "gray60", l0.truth, l0.mle, l0.mean)
## l0.plot

fill.color = lighten(pal[2], 1.3)
l1.plot <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, fill.color, pal[2], l1.truth, l1.mle, l1.mean)
l1.plot.hisse <- make.post.plot(log.df, 6, expression(lambda[1]), l1.min.max, "lightskyblue", "steelblue2", l1.truth, l1.mle, l1.mean)
## l1.plot

fill.color = lighten(pal[3], 1.3)
l2.plot <- make.post.plot(log.df, 7, expression(lambda[1][H]), l2.min.max, fill.color, pal[3], l2.truth, l2.mle, l2.mean)
l2.plot.hisse <- make.post.plot(log.df, 7, expression(lambda[1][H]), l2.min.max, "hotpink", "deeppink3", l2.truth, l2.mle, l2.mean)
## l2.plot

fill.color = lighten(pal[4], 1.3)
m0.plot <- make.post.plot(log.df, 8, expression(mu[0]), m0.min.max, fill.color, pal[4], m0.truth, m0.mle, m0.mean)
## m0.plot

fill.color = lighten(pal[5], 1.4)
m1.plot <- make.post.plot(log.df, 9, expression(mu[1]), m1.min.max, fill.color, pal[5], m1.truth, m1.mle, m1.mean)
## m1.plot 

fill.color = lighten(pal[6], 1.4)
m2.plot <- make.post.plot(log.df, 10, expression(mu[2]), m2.min.max, fill.color, pal[6], m2.truth, m2.mle, m2.mean)
## m2.plot 

fill.color = lighten(pal[7], 1.3)
q01.plot <- make.post.plot(log.df, 11, expression(q[0][1]), q01.min.max, fill.color, pal[7], q01.truth, q01.mle, q01.mean)
## q01.plot

fill.color = lighten(pal[8], 1.3)
q10.plot <- make.post.plot(log.df, 12, expression(q[1][0]), q10.min.max, fill.color, pal[8], q10.truth, q10.mle, q10.mean)
## q10.plot

fill.color = lighten(pal[9], 2.0)
q12.plot <- make.post.plot(log.df, 13, expression(q[1][2]), q12.min.max, fill.color, pal[9], q12.truth, q12.mle, q12.mean)
## q12.plot

fill.color = lighten(pal[10], 1.6)
q21.plot <- make.post.plot(log.df, 14, expression(q[2][1]), q21.min.max, fill.color, pal[10], q21.truth, q21.mle, q21.mean)
## q21.plot

pdf("plots/HiSSE_fixed_tree_on_HiSSE_HSDSEP_posteriors_120spp.pdf", width=9, height=7)
grid.arrange(l0.plot, l1.plot, l2.plot,
             m0.plot, m1.plot, m2.plot,
             q01.plot, q10.plot,
             q12.plot, q21.plot)
dev.off()

lay <- rbind(c(1,2,NA),c(4,5,6))
pdf("plots/BiSSE_HiSSE_lambdas_fixed_tree_120spp.pdf", width=8, height=4)
grid.arrange(l0.plot.bisse, l1.plot.bisse, l0.plot.hisse, l1.plot.hisse, l2.plot.hisse, ncol=3, layout_matrix=lay)
dev.off()

# ----- Plotting BiSSE_fixed_tree_SDSEP_SCM ----- #
# (below) if we want to compare it to diversitree... we see results are very close
# note that the difference is due to the fact that l0 and l1 share a prior in SSE that has
# its mean=.2, that is, in between the true value of l0(=.15) and l1(=.3),
                                        # while diversitree's SCM in this case uses the true values
pars <- c(.15, .3, .1, .1, .5, .5) # lambdas, mus, qs
set.seed(12345)
phy <- tree.bisse(pars, max.taxa=60, include.extinct=FALSE, x0=NA)
## node.truth.df <- as.data.frame(cbind(phy$node.label, phy$node.state))
## names(node.truth.df) <- c("ndname", "truestate")
## sampling.f <- c(1,1)
## lik <- make.bisse(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE)
## anc.states <- asr.marginal(lik, pars)
## anc.states <- as.data.frame(t(rbind(phy$node.label, anc.states)))
## names(anc.states) <- c("ndname", "p0", "p1")
## anc.states$"ndname" <- as.character(anc.states$"ndname")
## anc.states$"p0" <- as.numeric(as.character(anc.states$"p0"))
## anc.states$"p1" <- as.numeric(as.character(anc.states$"p1"))
## merged.df <- merge(log.df, anc.states, by="ndname")
## merged.df <- merge(merged.df, node.truth.df, by="ndname")
## plot(p1.y~p1.x, data=merged.df)

log.df <- read.table("BiSSE_fixed_tree_SDSEP_SCM_parsed.txt", header=FALSE, fill=TRUE)
names(log.df) <- c("ndname", "s0", "s1")
log.df$"s1"[is.na(log.df$"s1")] <- 0
rs <- rowSums(log.df[,c(2,3)])
log.df$"p0" <- log.df$"s0"/rs
log.df$"p1" <- log.df$"s1"/rs
log.df$"ndname" <- as.character(log.df$"ndname")
log.df$"p0" <- as.numeric(as.character(log.df$"p0"))
log.df$"p1" <- as.numeric(as.character(log.df$"p1"))
scm.df <- unname(as.matrix(log.df[order(match(log.df$ndname, phy$node.label)),c("p0", "p1")]))

pal <- brewer.pal(8, "Set1")
pal <- colorRampPalette(pal)(8)

pdf("plots/BiSSE_SDSEP_fixed_tree_60spp_SCM.pdf", width=5, height=5)
plot(phy, show.tip.label=FALSE, type="fan")
## nodelabels(pie=t(anc.states), cex=.5, piecol=c(pal[2],pal[6])) # from example_xml_input
nodelabels(pie=scm.df, cex=.5, piecol=c(pal[2],pal[6]))
tiplabels(phy$tip.state, frame="none", cex=.6, offset=1)
dev.off()

# ----- Plotting ModelAveraging_fixed_tree_BSVSSSDSEP ----- #
mod.avg.log <- read.table("../examples/ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP.log", header=TRUE)

bisse <- mod.avg.log[mod.avg.log$StateMask1==0 & mod.avg.log$StateMask2==0,]
bisse$Model <- "BiSSE"

hisse1 <- mod.avg.log[mod.avg.log$StateMask1==0 & (mod.avg.log$StateMask2==1 | mod.avg.log$StateMask2==2),]
hisse1$Model <- "HiSSE 1 (truth)"

hisse2 <- mod.avg.log[mod.avg.log$StateMask2==0 & (mod.avg.log$StateMask1==1 | mod.avg.log$StateMask1==2),]
hisse2$Model <- "HiSSE 1 (wrong)"

hisse3 <- mod.avg.log[(mod.avg.log$StateMask1==1 | mod.avg.log$StateMask1==2) & (mod.avg.log$StateMask2==1 | mod.avg.log$StateMask2==2),]
hisse3$Model <- "HiSSE 2"

avg.model <- mod.avg.log
avg.model$Model <- "Mean"

density.df <- rbind(bisse,hisse1,hisse2,hisse3,avg.model)

pal <- brewer.pal(8, "Set1")
pal <- colorRampPalette(pal)(8)

l0.density <- make.post.dens.plot(density.df, "Lambda1", expression(lambda[0]), pal, c(mean(density.df[density.df$Model=="BiSSE","Lambda1"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (truth)","Lambda1"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (wrong)","Lambda1"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 2","Lambda1"]),                                                                                       mean(density.df[density.df$Model=="Mean","Lambda1"])), 0.1)
## l0.density

l1.density <- make.post.dens.plot(density.df, "Lambda2", expression(lambda[1]), pal, c(mean(density.df[density.df$Model=="BiSSE","Lambda2"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (truth)","Lambda2"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (wrong)","Lambda2"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 2","Lambda2"]),
                                                                                       mean(density.df[density.df$Model=="Mean","Lambda2"])), 0.1)
## l1.density

l0h.density <- make.post.dens.plot(density.df, "Lambda3", expression(lambda[0][H]), pal, c(mean(density.df[density.df$Model=="BiSSE","Lambda3"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (truth)","Lambda3"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 1 (wrong)","Lambda3"]),
                                                                                       mean(density.df[density.df$Model=="HiSSE 2","Lambda3"]),
                                                                                       mean(density.df[density.df$Model=="Mean","Lambda3"])), 0.1)
## l0h.density

l1h.density <- make.post.dens.plot(density.df, "Lambda4", expression(lambda[1][H]), pal, c(mean(density.df[density.df$Model=="BiSSE", "Lambda4"]),
                                                                                                      mean(density.df[density.df$Model=="HiSSE 1 (truth)","Lambda4"]),
                                                                                                      mean(density.df[density.df$Model=="HiSSE 1 (wrong)","Lambda4"]),
                                                                                                      mean(density.df[density.df$Model=="HiSSE 2","Lambda4"]),
                                                                                                      mean(density.df[density.df$Model=="Mean","Lambda4"])), 0.5)
## l1h.density

pdf("plots/BSSVS_fixed_tree_120spp.pdf", width=6, height=7)
grid.arrange(l0.density, l1.density, l1h.density)
dev.off()
