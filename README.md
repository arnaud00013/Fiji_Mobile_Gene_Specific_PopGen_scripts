# Fiji_SNVs_and_metrics
This pipeline allows to call mobile genes SNVs and compute population genetics metrics of interest from mobile gene alignment with samples metagenomic reads.

## Pipeline
The user only needs to run Data_Exploration.py which will automatically run subsample_bam_v_df.r and curve_Cov_dnds_Ne_.r on every samples. It will also run Samples_comparisons.r at the end of the pipeline to compare the samples data. 

* subsample_bam_v_df.r


## References
1. A language and environment for statistical computing (R Foundation for Statistical Computing, Vienna, Austria, 2019). 
2. Brito, I. L. et al. Mobile genes in the human microbiome are structured from global to individual scales. Nature 535, 435-439, doi:10.1038/nature18927 (2016).
3. Eren, A. M. et al. Anvi'o: an advanced analysis and visualization platform for 'omics data. PeerJ 3, e1319, doi:10.7717/peerj.1319 (2015).
