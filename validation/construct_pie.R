# NOTE: THIS DOES NOT MATCH THE NODES PROPERLY
library(ape)

args = commandArgs(TRUE)
anc.state.csv = args[1]
phy.name = args[2]
pie.name = args[3]
#eg = "/Users/jeff/Documents/Research/Phylogenetics/SSE/validation/diversitree_data/rb-div_anc_states.csv"
#egg = read.csv(file=eg, header=TRUE, sep=",")

anc.states.0 = read.csv(file=anc.state.csv, header=TRUE, sep=",")
phy = read.tree(phy.name)

# See if we need to fill out the second row 
if (nrow(anc.states.0) == 1) {
    anc.states.1 = 1 - anc.states.0
    anc.states = rbind(anc.states.0, anc.states.1)  
} else {
    anc.states = anc.states.0
}

pdf(pie.name)
plot(phy, cex=.5, label.offset=0.2)
col = c("#004165", "#eaab00")
nodelabels(pie=t(anc.states), piecol=col, cex=.5)
