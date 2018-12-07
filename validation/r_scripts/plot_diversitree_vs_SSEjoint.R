library(ggplot2)

args <- commandArgs(TRUE)

# making sure all dirs are OK
output.dir <- args[3]
if (!dir.exists(output.dir)) {
    cat("Could not find directory to put graph in. Exiting...\n")
    q()
}

# making sure files exist
diversitree.asr.file <- args[1]
sse.asr.file <- args[2]
if (!(file.exists(diversitree.asr.file) && file.exists(sse.asr.file))) {
    cat("Could not find one of the stochastic mapping files (diversitree's or SSE's). Exiting...\n")
    q()
}

sse <- data.frame(t(read.csv(sse.asr.file, header=FALSE)))
div <- data.frame(t(read.csv(diversitree.asr.file, header=FALSE)))
div.0 <- div[,c(1,2)]
h <- c("node", "asr")
names(sse) <- h
names(div.0) <- h

merged.df <- merge(sse, div.0, by="node")
names(merged.df)[c(2,3)] <- c("sse", "div")

merged.df$div <- as.numeric(as.character(merged.df$div))
merged.df$sse <- as.numeric(as.character(merged.df$sse))

prefix <- args[4]

# graph
sse.div.plot <- ggplot(merged.df, aes(x=div, y=sse)) +
  geom_point(size=2) + ylab("Posterior probability in SSE") + xlab("Likelihood in diversitree") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(color="black"),
    axis.text.x = element_text(color="black", size=10),
    axis.text.y = element_text(color="black", size=10),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.key = element_blank()
  ) + geom_abline(slope=1, intercept=0, lty="dotted")

pdf(paste0(output.dir, prefix, "_sse_vs_diversitree_asr.pdf"), width=4.2, height=4)
                                        # plot(sse~div, data=merged.df, xlab="diversitree", ylab="SSE", main="Each dot is a node")
sse.div.plot
dev.off()
