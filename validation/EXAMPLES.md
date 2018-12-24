# Generating input of example .xmls    

This will guide you through the preparation of input for the .xmls in the examples/ folder.
You are going to need the following R packages:

* ape
* diversitree
* hisse
* ggplot2
* gridExtra
* RColorBrewer

## Running simulations in R (see examples_xml_inputs.R)    
### BiSSE_fixed_tree_SDSEP.xml and BiSSE_fixed_tree_HSDSEP.xml    
SDSEP and HSDSEP amount to the same analysis, we're doing them just to show the more general Java class with hidden states reduces to the original one.
Simulating and computing MLEs below.    

```
library(ape)
library(diversitree)
library(hisse)

pars <- c(.1, .5, .05, .05, .1, .1) # lambdas, mus, qs

set.seed(12345)
phy <- tree.bisse(pars, max.taxa=120, include.extinct=FALSE, x0=NA)

# tree (for .xml)
write.tree(phy) # see tree

# taxa (for .xml)
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data (for .xml)
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

# find MLE
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)
```

### HiSSE_fixed_tree_on_HiSSE_HSDSEP.xml and BiSSE_fixed_tree_on_HiSSE_HSDSEP.xml    
Conducting simulations with hidden trait (one hidden state), and computing MLEs.    

```
pars <- c(.1,  .1,  .5,  # lambda 0, 1, 2
.05, .05, .05, # mu 0, 1, 2
.1, 0.0, # q01, q02
.1, .1, # q10, q12
0.0, .1 # q20, q21
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=120, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree (for .xml)
write.tree(phy) # see tree

# taxa (for .xml)
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data (for .xml)
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

# finding MLE under BiSSE
pars.to.estimate <- pars <- c(.1, .3, .05, .05, .1, .1) # lambdas, mus, qs
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

# finding MLE under HiSSE
turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates,c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.no0B, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL)
```

### ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP.xml
Conducting simulations with hidden trait (one hidden trait with one hidden state just like *_on_HiSSE_HSDSEP.xml examples above).    

```
pars <- c(.1,  .1,  .5, .05, .05, .05, .1, 0.0, .1, .1, 0.0, .1) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=120, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree (for .xml)
write.tree(phy)

# taxa (for .xml)
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data (for .xml)
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")
```

After running .xml, remove header lines (starting with '#') from .log file, and save new file as 'ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP_noheader.log'.

### Stochastic character mapping on 60-sp tree under BiSSE

```
cd validation/
python scripts/parse_asm_treesfile.py ../examples/BiSSE_fixed_tree_SDSEP_SCM.trees 101 BiSSE_fixed_tree_SDSEP_SCM_parsed.txt
```

## Plotting all graphs in R (see examples_xml_plots.R)    
If you executed the previous steps (and ran BEAST), you will have produced

examples/BiSSE_fixed_tree_SDSEP.log
examples/BiSSE_fixed_tree_on_HiSSE_HSDSEP.log
examples/HiSSE_fixed_tree_on_HiSSE_HSDSEP.log
examples/ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP.log

You can find all the code to plot the posterior distributions and other graphs in examples_xml_plots.R.
We will just call the script and let it do all the work.    

```
mkdir /path/to/SSE/validation/plots
Rscript r_scripts/example_xml_plots.R /path/to/SSE/validation
```

