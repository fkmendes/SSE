library(stringr)
library(ape)
library(diversitree)
set.seed(1234)

args <- commandArgs(TRUE)

# making sure all dirs are OK
validation.dir <- args[1]
if (!dir.exists(validation.dir)) {
    cat("Could not find validate directory. Exiting...\n")
    q()
}
asr.dir <- args[2]
asr.full.path <- paste0(validation.dir, asr.dir)
if (!dir.exists(asr.full.path)) {
    cat("\nCreated asr/ directory.\n")
    dir.create(asr.full.path)
}

# set parameter values
pars = c(0.2, 0.4, 0.01, 0.1, 0.1, 0.4)
names(pars) = c("l0", "l1", "m0", "m1", "q01", "q10")
n.taxa = 22

phy = tree.bisse(pars, max.taxa=n.taxa, x0=0)
if (is.null(phy)) {
    cat("You did not get a tree for some reason (maybe because of seed generators in your OS). Exiting...\n")
    q()
}

## tree.str.beast <- str_replace_all(write.tree(phy), "nd[0-9]*", "")
tree.str.beast <- write.tree(phy)
cat("\nTree to do ASR in:")
cat(paste0("\n", tree.str.beast, "\n"))

## [1] "(((sp15:0.5701922606,(sp22:0.1174274481,sp23:0.1174274481):0.4527648125):5.46955786,((sp4:2.913008462,(sp16:0.4790358056,sp17:0.4790358056):2.433972656):1.72680138,sp2:4.639809842):1.399940278):8.039087646,((sp1:5.262858931,((((sp10:1.936988093,sp11:1.936988093):0.8700699862,((sp20:0.1813602217,sp21:0.1813602217):2.59756285,sp6:2.778923072):0.02813500652):0.1038009358,(sp14:1.103215563,(sp18:0.2976700868,sp19:0.2976700868):0.805545476):1.807643452):0.5229591127,sp3:3.433818127):1.829040804):1.760591904,((((sp8:1.951198056,sp9:1.951198056):0.153294648,sp7:2.104492704):0.5588707339,sp12:2.663363438):0.2401874525,sp5:2.90355089):4.119899945):7.055386931):0.0;"

# tip states
tips <- phy$tip.state + 1  # recoding states 0 and 1 into 1 and 2 for SSE

##  sp1  sp2  sp3  sp4  sp5  sp6  sp7  sp8  sp9 sp10 sp11 sp12 sp14 sp15 sp16 sp17 
##    1    1    1    2    1    1    2    1    1    1    1    1    1    2    2    1 
## sp18 sp19 sp20 sp21 sp22 sp23 
##    1    1    2    2    1    1 

cat("\nTip states:\n")
tip.string <- paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")
cat(paste0(tip.string, "\n"))

# node states (truth)
node.truth <- phy$node.state + 1  # recoding states 0 and 1 into 1 and 2 for comparison

##  nd1  nd2  nd6 nd22  nd7  nd9 nd11  nd3  nd4  nd8 nd10 nd12 nd15 nd16 nd17 nd13 
##    1    2    2    1    2    2    2    1    1    1    1    1    1    1    2    1 
## nd21  nd5 nd14 nd18 nd20 
##    1    1    1    1    1 

# computing likelihood
sampling.f <- c(1,1)
lik <- make.bisse(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE)
tree.lik <- lik(pars=pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE)
                                        # -63.001
cat(paste0("\ndiversitree's -lnL = ", tree.lik[1], "\n"))

# diversitree's ASR
anc.states <- asr.marginal(lik, pars)
anc.states.df <- as.data.frame(anc.states)
names(anc.states.df) <- phy$node.label

## sanity check
asr.df <- as.data.frame(matrix(nrow=phy$Nnode, ncol=4))
anc.states.t <- t(anc.states)
asr.df[,1] <- node.truth-1
asr.df[,2] <-
    node.depth.edgelength(phy)[(length(phy$tip.label)+1):length(node.depth.edgelength(phy))] # I believe I'm getting the heights in the right order from the piecharts (below)...
asr.df[asr.df$V1==0,3] <- anc.states.t[asr.df$V1==0,1]
asr.df[asr.df$V1==1,3] <- anc.states.t[asr.df$V1==1,2]
asr.df[anc.states[1,] > anc.states[2,], 4] <- 0 # when post prob of 0 > that of 1, state 0
asr.df[anc.states[1,] < anc.states[2,], 4] <- 1 # when post prob of 0 < that of 0, state 1
names(asr.df) <- c("truth", "depth", "truthLk", "highestLk")

acc <- sum(asr.df[,1]==asr.df[,4])/nrow(asr.df)
cat(paste0("\nAccuracy of diversitree's ASR by taking ancestral state with highest likelihood: ", acc, "\n"))

# a couple of graphs
cat("\nPlotting sanity check graph (node height by likelihood of true ancestral state) inside asr/\n")
pdf(paste0(asr.full.path, "/depth_by_truelk.pdf"), width=5, height=5)
plot(truthLk~depth, data=asr.df, xlab="Node height", ylab="Likelihood of true state", pch=20, bty="n", xlim=c(0,15), ylim=c(0,1.1), xaxs="i", yaxs="i")
dev.off()

cat("\nPlotting diversitree's ASR piecharts inside asr/.\n")
pdf(paste0(asr.full.path, "/asr_piechart.pdf"), width=5, height=5)
plot(phy, cex=.5, label.offset=0.2)
nodelabels(pie=t(anc.states), cex=.5)
dev.off()

# printing!
cat(tree.str.beast, file=paste0(asr.full.path, "/simulated_tree.txt"))
cat(tip.string, file=paste0(asr.full.path, "/simulated_tip_states.txt"))
write.csv(t(node.truth), file=paste0(asr.full.path, "/true_ancestral_states.txt"), quote=FALSE, row.names=FALSE)
write.csv(anc.states.df, file=paste0(asr.full.path, "/diversitree_ancestral_states.txt"), quote=FALSE, row.names=FALSE)
