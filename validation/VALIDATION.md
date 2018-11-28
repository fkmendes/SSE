# SSE validation

This will guide you through the validation procedures we conducted (as well as generate all figures).
You are going to need the following R packages:

* ape
* diversitree

## Ancestral state reconstruction via stochastic character mapping    
In the SSE package, we implement the same stochastic character mapping approach proposed in Freyman and H&ouml;hna (2017).
We compare our implementation to diversitree's like they do (and produce the same figure), and also compare the two mapping methods (sampling at internal nodes only vs. sampling along the entire tree in fixed-size chunks).

### Preparing input for stochastic character mapping (tree, tip data, internal node data, diversitree reconstructions) in R

The following command will simulate a tree and tip states given parameters λ0=0.2, λ1=0.4, µ0=0.01, µ<sub>1</sub>=0.1, q01=0.1, q10=0.4 (we use a fixed seed to obtain the same graph in Freyman and H&ouml;hna, 2017). The R script this command executes also performs ancestral state reconstruction using diversitree.

Below, "/path/to/validation" should be replaced with your local path to the validation folder that comes with this git repository.

```
Rscript /path/to/validation/r_scripts/prepare_BiSSE_ASR_input.R /path/to/validation/ asr
```

A few files will be produced and put into folder /asr, and some information will be printed on the screen (and used in the following validation steps).