# Fiji_SNVs_and_metrics
This pipeline allows to call mobile genes SNVs and compute population genetics metrics of interest from mobile gene alignment with samples metagenomic reads. The scripts are built to handle specifically to the FijiComp project files (http://fijicomp.bme.cornell.edu//data.html) and perform gene-specific molecular evolution analysis based on metagenomics mapping. 

## Pipeline
-Recquired Operating System : Unix-based operating systems such as Mac OS or Linux distributions.
-Recquired Python version : the python version should be compatible with anv'io (default=python3.6)(INDICATE THE PATH OF THE PYTHON INTERPRETER AT THE BEGINNING OF Data_Exploration.py)
-Recquired Python librairies : numpy, pandas, sys, time, multiprocess & os
-Recquired r librairies : ggplot2, ggdendro, nlme, grid, reshape2 & seqinr
-Recquired virtual environment : anvio5
-Workspace repertory architecture : 
	-Choose the repertory that will be the Main_Workspace from which the python3.6 interpreter will execute the script
	-In this Workspace :
		-put all the samples directory (see required content below)
		-put the 4 scripts file which are "Data_Exploration.py",                        "subsample_bam_v_df.r", "curve_Cov_dnds_Ne_.r" and "Samples_comparisons.r"
		-put the files “contigs.fa” (the file containing the reference sequence of      all the genes of interest in fasta format; Gene ids need to have the anvio      format which you can get with anvi-script-reformat-fasta),                      “conversion_fasta.id.txt” (the report file of anvi-script-reformat-fasta        indicating the correspondency between your gene_ids and simplified              anvio_gene_ids) and "contigs.db” (the anvio contigs database file which        can be generated from a fasta file of the genes reference sequence and           the anvio script anvi-gen-contigs-database)
		-put the Metadata file "Metadata_for_mobile_gene_dataset.csv" and the file      "HGTgenes_lookup.txt" for gene id conversion 
	-In each sample directory, put the non-sorted sample ".bam" files                 representing the mapping of the metagenomic reads of one individual on the     tested genes (the mapping file should have the extension                       ".hgt.mapped.merged.99.50.bam").

-To execute the pipeline:
	-cd $PATH_Main_Workspace
	-chmod u+x Data_Exploration.py subsample_bam_v_df.r curve_Cov_dnds_Ne_.r        Samples_comparisons.r
	-load modules AND anvio environment
		-module load nixpkgs/16.09 gcc/5.4.0 gsl prodigal/2.6.3 hmmer/3.1b2               samtools/1.3.1 hdf5/1.8.18 r/3.4.0 python/3.6
		-source $ANVIO_VIRTUAL_ENVIRONMENT_PATH/bin/activate
	-./Data_Exploration.py

## References
1. A language and environment for statistical computing (R Foundation for Statistical Computing, Vienna, Austria, 2019). 
2. Brito, I. L. et al. Mobile genes in the human microbiome are structured from global to individual scales. Nature 535, 435-439, doi:10.1038/nature18927 (2016).
3. Eren, A. M. et al. Anvi'o: an advanced analysis and visualization platform for 'omics data. PeerJ 3, e1319, doi:10.7717/peerj.1319 (2015).
