# Generating input of example .xmls    

This will guide you through the preparation of input for the .xmls in the examples/ folder.    

## Running simulations (see examples_xml_inputs.R)

```
library(ape)
library(diversitree)
library(hisse)

pars <- c(.15, .3, .1, .1, .1, .1) # lambdas, mus, qs

set.seed(12345)
phy <- tree.bisse(pars, max.taxa=60, include.extinct=FALSE, x0=NA)

# tree
write.tree(phy) # see tree

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

# find MLE
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)
```
