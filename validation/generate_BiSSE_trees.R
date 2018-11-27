library(ape)
library(diversitree)
set.seed(1234)

args = commandArgs(TRUE)
exp.name = args[1]
dir = "diversitree_data/"
root = "/Users/jeff/Documents/Research/Phylogenetics/SSE/validation/"
ifelse(!dir.exists(dir), dir.create(dir), FALSE)

## HELPERS
# Compares the states in the truth and the reconstruction for given states/indices of interest
compare.states = function(truth, pred, indices) {
    num.correct = 0
    for (i in indices) {
        if (truth[i] == pred[i]) {
            num.correct  = num.correct + 1
        }
    }
    num.samples = length(indices)
    return (num.correct / num.samples)
}


# Set parameters
lambda = 1/20
mu = 1/80
q = 1/20
sim.time = 50
# pars = c(lambda, lambda, mu, mu, q, q)

#pars = c(0.5, 0.4, 0.02, 0.1, 0.1, 0.02)  # from RevBayes exp
#pars = c(0.04, 0.08, 0.01, 0.02, 0.04, 0.01)
#pars = c(0.08, 0.08, 0.01, 0.01, 0.01, 0.01)
pars = c(0.2, 0.4, 0.01, 0.1, 0.1, 0.4)
n.taxa = 22
#n.taxa = 5
names(pars) = c("l0", "l1", "m0", "m1", "q01", "q10")
pars

# Set up a handcrafted tree 
#trstr = "((Human:2.0,Chimp:2.0)nd1:2.0,(Gorilla:2.0,Pizza:2.0)nd2:2.0)nd3:0.0;"
#phy = read.tree(file="", text=trstr)
#plot(phy)
#states = c(1, 1, 0, 0)
#names(states) = phy$tip.label

# Get a ground truth tree
phy = tree.bisse(pars, max.taxa=n.taxa, x0=0)
plot(phy)
if (is.null(phy)) {
    print("bad tree")
    q()
}
tips = phy$tip.state + 1  # to make the two states 1 and 2
print("Tree tips")
print(tips)
print("Tree in Newick format")
write.tree(phy)
tree.file.name = paste0(dir, exp.name, ".tree")
write.tree(phy, file=tree.file.name)
ntips = length(tips)
node.truth = phy$node.state + 1 # to make the two states 1 and 2

print("Node truths")

print(node.truth)
node.truth.save = node.truth
node.truth.names = names(node.truth)
node.truth.save = t(node.truth.save)
node.truth.save = rbind(node.truth.names, node.truth.save)
node.file.name = paste0(dir, exp.name, "-node_truth.csv")
write.table(node.truth.save, file=node.file.name, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)


tips.save = tips
tips.save.names = names(tips.save)
tips.save = t(tips.save)
tips.save = rbind(tips.save.names, tips.save)
tips.file.name = paste0(dir, exp.name, "-tips.csv")
write.table(tips.save, file=tips.file.name, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)


beast.data = ""
for (i in 1:length(tips)) {
    beast.data = paste0(beast.data, names(tips[i]), "=", tips[i], ",")
}   
beast.data = substr(beast.data, 0, nchar(beast.data) - 1) # to remove the extra comma added
beast.data.file.name = paste0(root, dir, exp.name, "-beast_str.txt")
cat(beast.data, file=beast.data.file.name)


# Calculate likelyhood on tree
sampling.f = c(1,1)
lik = make.bisse(tree=phy, states=phy$tip.state, sampling.f=sampling.f, strict=FALSE)
tree.lik = lik(pars=pars, root=ROOT.FLAT, root.p=NULL, intermediates=TRUE, condition.surv=FALSE) 
# prior on root is the weighted average of D0 and D1, i.e., ROOT.OBS = D = D0 * (D0/(D0+D1)) + D1 * (D1/(D0+D1))
print("Tree log likeyhood")
print(tree.lik[1])

# Calculate ancestral MLE
asr.marginal = asr.marginal(lik, pars)


# Compare truth tips and MLE tips for accuracy
node.marginal = vector("list", phy$Nnode)
for (i in 1:phy$Nnode) {
    if (asr.marginal[1,i] > asr.marginal[2,i]) {
        node.marginal[i] = 1 - 1
    } else {
        node.marginal[i] = 2 - 1
    }
}

acc = compare.states(node.truth, node.marginal, 1:phy$Nnode)
print("Compare asr to ground truth")
cat("Total accuracy: ", acc, "\n")


# Compare accuracy on the more recent and ancient nodes
all.node.depth = node.depth.edgelength(phy)  # Note this includes tips too... time from beginning of time 
node.sorted = order(all.node.depth)[1 : phy$Nnode] # gives us the order of internal nodes from most ancient to least
node.ancient = node.sorted[1 : as.integer(phy$Nnode/4)] - ntips # quartile of most ancient 
node.recent = node.sorted[as.integer(phy$Nnode * 3/4) : phy$Nnode] - ntips # quartile of least ancient

acc.ancient = compare.states(node.truth, node.marginal, node.ancient)
cat("Ancient accuracy: ", acc.ancient, "\n")
acc.recent = compare.states(node.truth, node.marginal, node.recent)
cat("Recent accuracy: ", acc.recent, "\n")


# Write the results out... merge the marginal calculates with their labels
print("asr marginal likelyhoods")
asr.marginal.labeled = rbind(phy$node.label, asr.marginal)
print(asr.marginal.labeled)
asr.file.name = paste0(dir, exp.name, "-div_anc_states.csv")
write.table(asr.marginal.labeled, file = asr.file.name, row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)


#attributes(phy)
#print(phy$edge)
