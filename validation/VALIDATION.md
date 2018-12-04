# SSE validation

This will guide you through the validation procedures we conducted (as well as generate all figures).
You are going to need the following R packages:

* ape
* diversitree

## Ancestral state reconstruction via stochastic character mapping    
In the SSE package, we implement the same stochastic character mapping approach proposed in Freyman and H&ouml;hna (2017).
We compare our implementation to diversitree's like they do (and produce the same figure), and also compare the two mapping methods (sampling at internal nodes only vs. sampling along the entire tree in fixed-size chunks).

### Preparing input for stochastic character mapping (tree, tip data, internal node data, diversitree reconstructions) in R

The following command will simulate a tree and tip states given parameters λ<sub>0</sub>=0.2, λ<sub>1</sub>=0.4, µ<sub>0</sub>=0.01, µ<sub>1</sub>=0.1, q<sub>01</sub>=0.1, <sub>q10</sub>=0.4 (we use a fixed seed to obtain the same graph in Freyman and H&ouml;hna, 2017). The R script this command executes also performs ancestral state reconstruction using diversitree.

Below, "/path/to/validation" should be replaced with your local path to the validation folder that comes with this git repository.

```
Rscript /path/to/validation/r_scripts/prepare_BiSSE_ASR_input.R /path/to/validation/ asr
```

A few files will be produced and put into folder /asr, and some information will be printed on the screen (and used in the following validation steps).

### Running stochastic character mapping (both only on internal nodes and along branches as well) in Java

Below, "/path/to/SSE.jar" should be replaced with the full path to the SSE.jar file.

```
java -cp /path/to/SSE.jar validation.SDSEPValidationBiSSEASR asr/
```

### Plotting diversitree's ancestral state reconstruction vs. SSE packages's
Now we show that our ancestral state reconstructions matches diversitree's for BiSSE. The R script we call below produces two graphs in the specified output folder (the graph file names are prefixed by "bisse_joint" and "bisse_stochastic").

```
Rscript /path/to/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/diversitree_ancestral_states.txt /path/to/SSE_asr_joint.csv /path/to/output/folder/ bisse_joint
Rscript /path/to/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/diversitree_ancestral_states.txt /path/to/SSE_asr_stoc.csv /path/to/output/folder/ bisse_stoc
```