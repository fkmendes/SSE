# Generating input of example .xmls    

This will guide you through the preparation of input for the .xmls in the examples/ folder.
You are going to need the following R packages:

* ape
* diversitree
* hisse

## Running simulations (see examples_xml_inputs.R)

Preparing inputs for BiSSE_fixed_tree_SDSEP.xml and BiSSE_fixed_tree_HDSEP.xml (they amount to the same analysis, we're doing them just to show the more general Java class with hidden states reduces to the original one).

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
