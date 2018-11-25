# Generating input of example .xmls    

This will guide you through the preparation of input for the .xmls in the examples/ folder.
You are going to need the following R packages:

* ape
* diversitree
* hisse

## Running simulations in R (see examples_xml_inputs.R)    
### BiSSE_fixed_tree_SDSEP.xml and BiSSE_fixed_tree_HSDSEP.xml    
SDSEP and HSDSEP amount to the same analysis, we're doing them just to show the more general Java class with hidden states reduces to the original one.    

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

### BiSSE_fixed_tree_on_HiSSE_HSDSEP.xml    
Conducting simulations with hidden trait (one hidden state).    

```
pars <- c(.1,  .1,  .3,  # lambda 1, 2, 3
.05, .05, .05, # mu 1, 2, 3
.1, 0.0, # q12, q13
.1, .1, # q21, q23
0.0, .1 # q31, q32
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=60, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree
write.tree(phy)

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data
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

## Plotting all graphs in R (see examples_xml_plots.R)    

You can find all the code to plot the posterior distributions (of all parameters of all .xmls we produced in the steps above) in examples_xml_plots.R.
We will just call the script and let it do all the work.    

```
Rscript example_xml_plots.R /path/to/SSE/validation
```