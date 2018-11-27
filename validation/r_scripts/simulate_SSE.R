# This file has a few important outputs
#  1. sim_data/params.csv
#  2. sim_data/classe.tree
#  3. sim_data/tips.csv - list of node to character mappings

# Use diversitree:::default.argnames.classe(3) to see arg names
# simulates and returns a phylogeny given parameters and 
# initial state for max.t amount of time 
# No cap on the number of taxa
# Do not need information on extinct taxa
# Sample initial state from equilibrium distribution

library(diversitree)
library(ape)
set.seed(1234)

# write tip data file
write.tips <- function(states, file.name){
  num_taxa = length(states)
  a = FALSE
  for (i in 1:num_taxa) {
    cat(c(names(states[i]), states[i]), file=file.name, append=a, sep=",")
    a = TRUE
    cat("\n", append=a, file=file.name)
  }
}

sim.sse <- function(output.dir, prefix, pars, is.classe=TRUE) {
  # Note: if the extinction rates exceed the speciations rates significantly, 
  #       the tree will die before speciation 
  if (is.classe) {
    phy = tree.classe(pars, max.t=22, max.taxa=1000, 
                      include.extinct=FALSE, x0=NA)  
  } else {
    phy = tree.bisse(pars, max.t=22, max.taxa=1000, 
               include.extinct=FALSE, x0=NA)  
  }
  if (is.null(phy)) {
    print(pars)
    print("Probably extinction rates are too high compared to speciation rates.")
  }
  else {
    states = phy$tip.state
    
    save.path = paste(output.dir, prefix, sep = "/")
    write.tips(states, file.name=paste(save.path, "tips.csv", sep = "_"))
    write.tree(phy, file=paste(save.path, ".tree", sep = ""))
    write.table(pars, file=paste(save.path, "params.csv", sep = "_"), 
                row.names=FALSE, col.names=FALSE, eol=",")
  }
}



# ------------- Test Run -----------------------
# output.dir = "sim_data"
# prefix = "foo"
# lambda_111 = 0.2
# lambda_112 = 0.000001
# lambda_122 = 0.000001
# lambda_211 = 0.000001
# lambda_212 = 0.000001
# lambda_222 = 0.4
# 
# mu_1 = 0.01
# mu_2 = 0.1
# 
# q_12 = 0.1
# q_21 = 0.4
# pars = c(lambda_111, lambda_112, lambda_122,
#          lambda_211, lambda_212, lambda_222,
#          mu_1, mu_2, q_12, q_21)
# pars.bisse = c(lambda_111, lambda_222, mu_1, mu_2, q_12, q_21)
# 
# sim.sse(output.dir, prefix, pars, is.classe = TRUE)
# sim.sse(output.dir, prefix, pars.bisse, is.classe = FALSE)

# ------------- True Run --------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  print("Not enough args.")
}

output.dir = args[1]
prefix = args[2]
pars = c(args[3:length(args)])
pars = as.numeric(pars)
is.classe = length(pars) != 6

sim.sse(output.dir, prefix, pars, is.classe = is.classe)
