# SSE validation

This will guide you through the validation procedures we conducted (as well as generate all figures).
You are going to need the following R packages:

* ape
* diversitree

## (2) Ancestral state reconstruction via stochastic character mapping    
In the SSE package, we implement the same stochastic character mapping approach proposed in Freyman and H&ouml;hna (2017).
We compare our implementation to diversitree's like they do (and produce the same figure), and also compare the two mapping methods (sampling at internal nodes only vs. sampling along the entire tree in fixed-size chunks).

### (2.1) Preparing input for stochastic character mapping (tree, tip data, internal node data, diversitree reconstructions) in R

The following command will simulate a tree and tip states given parameters λ<sub>0</sub>=0.2, λ<sub>1</sub>=0.4, µ<sub>0</sub>=0.01, µ<sub>1</sub>=0.1, q<sub>01</sub>=0.1, <sub>q10</sub>=0.4 (we use a fixed seed to obtain the same graph in Freyman and H&ouml;hna, 2017).
The R script this command executes also performs ancestral state reconstruction using diversitree.

Below, "/path/to/SSE/validation" should be replaced with your local path to the validation folder that comes with this git repository.

```
Rscript /path/to/SSE/validation/r_scripts/prepare_BiSSE_ASR_input.R /path/to/validation/ asr
```

The following files will be produced:

* asr/asr_piechart.pdf
* asr/simulated_tree.txt
* asr/simulated_tip_states.txt
* asr/true_ancestral_states.txt
* asr/diversitree_ancestral_states.txt

### (2.2) Running stochastic character mapping (both only on internal nodes and along branches as well) in Java

Below, "/path/to/SSE.jar" should be replaced with the full path to the SSE.jar file.

```
java -cp /path/to/SSE.jar validation.SDSEPValidationBiSSEASR asr/
```

The following files will be produced:

* asr/SSE_asr_joint.csv
* asr/SSE_asr_stoc.csv

### (2.3) Plotting diversitree's ancestral state reconstruction vs. SSE packages's
Now we show that our ancestral state reconstructions matches diversitree's for BiSSE.

```
Rscript /path/to/SSE/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/SSE/validation/asr/diversitree_ancestral_states.txt /path/to/SSE/validation/asr/SSE_asr_joint.csv /path/to/SSE/validation/plots/ bisse_joint
Rscript /path/to/SSE/validation/r_scripts/plot_diversitree_vs_SSEjoint.R /path/to/SSE/validation/asr/diversitree_ancestral_states.txt /path/to/SSE/validation/asr/SSE_asr_stoc.csv /path/to/SSE/validation/plots/ bisse_stoc
```

The R script below will produce two graph in the "plots/" folder:

* plots/bisse_joint_sse_vs_diversitree_asr.pdf
* plots/bisse_stoc_sse_vs_diversitree_asr.pdf
