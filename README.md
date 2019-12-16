### Introduction
Code and data for analyses of ribosomal engagement estimated with parallel TRAP-seq and RNA-seq experiments during stem cell based neuronal differentiation. Look up X for details.

### Data
Before you run the code in the notebooks, unpack the compressed data.
First, concatenate data chunks:
```bash
cd /downloaded_github_folder/data
cat data.tar.gz.part_* > data.tar.gz
```
Then, unpack:
```bash
tar -xzf data.tar.gz
```
Data required to run the R notebooks is fully available in the *data* folder.

### R notebooks
1. **estimate_RE.Rmd**

   Ribosomal engagement is calculated from the TRAP-seq and RNA-seq. Quantile normalization is applied to correct for scaling bias between replicates. Stats between replicates is obtained for downstream analyses.
2. **analyze_RBPs.Rmd**

   Enrichments of RBPs binding elements in 3'UTRs are found for mRNAs with high and low ribosomal engagement. Code in *binding_elements* folder is used to scan 3'UTRs with RBP motifs downloaded from ATtRACT database. Identified binding elements are further filtered based on the conservation score among 100 vertebrates to enrich for functional binding elements.
3. **train_RF_classifier.Rmd**

   Random forest classifier is developed for probing the predictive power of convserved 3'UTR binding elements. Each mRNA is classified into either low or high ribosomal engagement groups. 
4. **get_seq_features.Rmd**

   Sequence features like GC content or length are obtained for 5'-UTR, CDS and 3'-UTR regions.
5. **analyze_seq_features.Rmd**

   Ribosomal engagement is correlated to sequence features.
6. **targets_of_TFs.Rmd**

   Shifts in mRNA steady state abundance of transcription factors targets defined in iRegulon database is calculated for translationally regulated transcription factors.
7. **analyze_lncRNAs.Rmd**

   Ribosomal engagement of lncRNAs is estimated based on TRAP-seq. Then, lncRNAs defined to be translated with TRAP-seq are compared to Ribo-seq measurements in heart and other stem cell derived neurons.
