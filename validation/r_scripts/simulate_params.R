library(ggplot2)

sample.parameter <- function(output.dir, prefix, n.sim, param.name) {
    # Setting up arbitrary simulation values, and making objects
    ps <- rlnorm(20000, meanlog=0, sdlog=1) # plotting samples: default mean and sd in log scale, and =0 and 1, respectively
    s <- rlnorm(n.sim, meanlog=0, sdlog=1) # actual samples                                        
    p.df <- data.frame(ps); names(p.df) <- "value"                                         
    df <- data.frame(s); names(df) <- "value"
    r <- range(ps) # min and max from samples
    xaxis.max <- 15

    plot.1 <- ggplot(df, aes(x=value)) +
        geom_histogram(aes(y=stat(density)), alpha=.4, bins=100) +
        stat_function(fun=dlnorm,
                      args=list(r[1]:r[2], meanlog=mean(log(p.df$value)), sdlog=sd(log(p.df$value)))
                      ) +
        scale_x_continuous(breaks=seq(0, xaxis.max, by=1), limits=c(0, xaxis.max)) +
        xlab("Parameter value") + ylab("Density") + 
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.line = element_line(),
            axis.ticks = element_line(color="black"),
            axis.text.x = element_text(color="black", size=10),
            axis.text.y = element_text(color="black", size=10),
            axis.title.x = element_text(size=12),
            axis.title.y = element_text(size=12)
        )
  
    png(paste0(output.dir, prefix, "_", param.name, ".png"), width=10, height=6, unit="cm", res=300)
    plot.1
    dev.off()

    return(df)
}

sample.parameters <- function(output.dir, prefix, n.sim, param.names) {
    num_params = length(param.names)
    params = vector("list", length = num_params)

    for (i in 1:num_params) {
        param.name = param.names[i]
        p = sample.parameter(output.dir, prefix, n.sim, param.name)
        params[[i]] = p
    }
    params.df = as.data.frame(params)
    save.path = paste0(output.dir, prefix)
    save.path = paste(save.path, "all_params.csv", sep="_")

    write.table(params.df, file=save.path, sep=",",
                row.names=FALSE, col.names=FALSE) # This csv has all of the parameters for all the simulations together to be used by next script.
  # It should not be used for any other purposes
}

# Running simulations

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
    print("Not enough args. Output directory, prefix, num simulations, parameter names.")
}

output.dir = args[1]
prefix = args[2]
n.sim = as.numeric(args[3])
param.names = args[4:length(args)]  

print(output.dir)
print(prefix)
print(n.sim)
print(param.names)

sample.parameters(output.dir, prefix, n.sim, param.names)
