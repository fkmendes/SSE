library(ggplot2)
library(sjPlot)
library(RColorBrewer)
library(stats)
library(gridExtra)

options(digits=2)

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])

csv.dir <- args[2]

make.regression.plot <- function(true.param.name, beast.param.mean, beast.param.lower, beast.param.upper, true.df, beast.df, hpd=FALSE, param.label) {
    true.name = as.character(true.param.name)
    beast.name = as.character(beast.param.mean)
    beast.lower = as.character(beast.param.lower)
    beast.upper = as.character(beast.param.upper)
    x = as.numeric(true.df[,names(true.df)==true.name])
    min.x = 0
    ## min.x = min(x)
    max.x = max(x)
    n = length(x)
    ## max.x = sort(x, partial=n-1)[n-2] # throwing out one outlier
    y = as.numeric(beast.df[,names(beast.df)==beast.name])
    lower = as.numeric(beast.df[,names(beast.df)==beast.lower])
    upper = as.numeric(beast.df[,names(beast.df)==beast.upper])
    ## min.y = min(lower)
    ## max.y = max(upper)
    ## min.y = min(y)
    min.y = 0
    max.y = max(y)
    reg.df = data.frame(cbind(x,y,lower,upper))
    print(reg.df)

    plot = ggplot() + geom_point(data=reg.df, mapping=aes(x=x, y=y), shape=20) + coord_cartesian(ylim=c(min.y, max.y)) +
        xlab(param.label) + ylab("Posterior mean") + geom_abline(slope=1, linetype="dotted") +
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
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), limits=c(min.x,max.y)) 

    if (hpd) { plot = plot + geom_linerange(data=reg.df, mapping=aes(x=x, ymax=upper, ymin=lower), color="lightgray", alpha=.4) }
    return(plot)
}

# BiSSE
true.names <- unlist(strsplit("l0,l1,m0,m1,q01,q10", split=","))
beast.names <- unlist(strsplit("Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix1,FlatQMatrix2", split=","))
param.labels <- c(expression(paste("Simulated ", lambda[0])),
                  expression(paste("Simulated ",lambda[1])),
                  expression(paste("Simulated ",mu[0])),
                  expression(paste("Simulated ",mu[1])),
                  expression(paste("Simulated ",q[0][1])),
                  expression(paste("Simulated ",q[1][0])))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names), paste0("lower",beast.names), paste0("upper",beast.names)))
names(name.df) <- c("true.names", "beast.names", "beast.lower", "beast.upper")

load(paste0(csv.dir, "hpds_bisse.RData"))
true.df <- read.table(paste0(csv.dir, "data_param_tree_bisse.csv"), sep="|", head=TRUE)

large.idxs <- true.df[,"ntips"]>=200
true.df <- true.df[large.idxs,]
df <- df[large.idxs,]

all.plots <- vector("list", nrow(name.df))
for (r in 1:nrow(name.df)) {
    all.plots[[r]] = make.regression.plot(name.df$true.names[r], name.df$beast.names[r], name.df$beast.lower[r], name.df$beast.upper[r], true.df, df, hpd=TRUE, param.labels[r])
}
## plot_grid(all.plots)

pdf("plots/BiSSE_calibrated_validation_200+spp.pdf", width=6, height=7)
plot_grid(all.plots)
dev.off()

# ClaSSE
true.names <- unlist(strsplit("l_111,l_313,l_312,m1,m2,m3,q01,q02,q10,q12,q20,q21", split=","))
beast.names <- unlist(strsplit("SympatricRate,SubsympatricRate,VicariantRate,Mu1,Mu2,Mu3,FlatQMatrix1,FlatQMatrix2,FlatQMatrix3,FlatQMatrix4,FlatQMatrix5,FlatQMatrix6", split=","))
name.df <- data.frame(cbind(true.names, paste0("mean",beast.names), paste0("lower",beast.names), paste0("upper",beast.names)))
names(name.df) <- c("true.names", "beast.names", "beast.lower", "beast.upper")

load(paste0(csv.dir, "hpds_classe.RData"))
true.df <- read.table(paste0(csv.dir, "data_param_tree_classe.csv"), sep="|", head=TRUE)

param.labels <- c(expression(paste("Simulated ", lambda["sym"])),
                  expression(paste("Simulated ", lambda["subsym"])),
                  expression(paste("Simulated ", lambda["vic"])),
                  expression(paste("Simulated ", mu[0])),
                  expression(paste("Simulated ", mu[1])),
                  expression(paste("Simulated ", mu[2])),
                  expression(paste("Simulated ", q[0][1])),
                  expression(paste("Simulated ", q[0][2])),
                  expression(paste("Simulated ", q[1][0])),
                  expression(paste("Simulated ", q[1][2])),
                  expression(paste("Simulated ", q[2][0])),
                  expression(paste("Simulated ", q[2][1])))
                  
large.idxs <- true.df[,"ntips"]>=200
true.df <- true.df[large.idxs,]
df <- df[large.idxs,]

all.plots <- vector("list", nrow(name.df))
for (r in 1:nrow(name.df)) {
    all.plots[[r]] = make.regression.plot(name.df$true.names[r], name.df$beast.names[r], name.df$beast.lower[r], name.df$beast.upper[r], true.df, df, hpd=TRUE, param.labels[r])
}
plot_grid(all.plots)

pdf("plots/ClaSSE_calibrated_validation_200+spp.pdf", width=9.5, height=10.5)
plot_grid(all.plots)
dev.off()
