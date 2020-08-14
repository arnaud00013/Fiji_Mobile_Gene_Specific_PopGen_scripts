#@Author=Arnaud NG
#Script arguments :
#snv_matrix_filepath : the path of the SNVs matrix file that contains data to plot
#sample_name : the name of the metagenomic sample

#load libraries and read fasta files + SNV profile file
#ENLEVER le samtools bin /Users/arnaudng/miniconda3/bin/ ET les chemins de repertoire par defaut

library("ggplot2")
library("seqinr")
library("nlme")

snv_matrix_filepath <- as.character(commandArgs(TRUE)[1])   #Find the path of the SNV matrix file
#obtain repertory path from absolute path of the SNV matrix file
indexes <- gregexpr(pattern = "/",text = snv_matrix_filepath)[[1]]
rep_path <- substr(snv_matrix_filepath,start = 1,stop=indexes[length(indexes)])
wp_path <- substr(snv_matrix_filepath,start = 1,stop=indexes[length(indexes)-1])
root_workspace <- substr(snv_matrix_filepath,start = 1,stop=indexes[length(indexes)-2])
the_sample_name <- as.character(commandArgs(TRUE)[2])

#system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/Users/arnaudng/miniconda3/bin/samtools idxstats ", wp_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam > ",wp_path,"Read_count_per_gene.csv"),intern = FALSE,wait = TRUE)
#system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools idxstats ", wp_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam > ",wp_path,"Read_count_per_gene.csv"),intern = FALSE,wait = TRUE)
#df_the_sample_genes_length_and_nb_reads <- read.csv2(file=paste0(wp_path,"Read_count_per_gene.csv"),sep='\t',header = FALSE, stringsAsFactors=FALSE)


#Read the mobile genes fasta file and store data
fasta_mbg_seqs <- read.fasta(file = paste0(root_workspace,"contigs.fa"), as.string = TRUE,forceDNAtolower = FALSE)
df_SNV <- read.csv2(file=snv_matrix_filepath,sep='\t',stringsAsFactors=FALSE,header = TRUE) #load the SNVs matrix file
df_SNV <- df_SNV[,c(1:6,17,21:25,27:32)]
  
#make sure the data frame variables have the good type
df_SNV$departure_from_reference <- as.numeric(df_SNV$departure_from_reference)
df_SNV$departure_from_consensus <- as.numeric(df_SNV$departure_from_consensus)
df_SNV$n2n1ratio <- as.numeric(df_SNV$n2n1ratio)
df_SNV$entropy <- as.numeric(df_SNV$entropy)

#exclude false-positives SNVs in df_SNV
df_SNV <- subset(x = df_SNV, subset = (((A>0)+(C>0)+(G>0)+(T>0)) >= 2))

#create a dataframe containing the genes of interest having a significant coverage 
df_genes_of_interest <- read.csv(paste0(wp_path,"gene_IDs_of_interest_",the_sample_name,".txt"),header = FALSE,sep='\t',stringsAsFactors = FALSE)

#get the contigs having a significant coverage
df_depth <- read.csv(paste0(wp_path,"depth.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
df_depth_over10 <- read.csv(paste0(wp_path,"genes_with_depth_over10_",the_sample_name,".csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)

#Update df_genes_of_interest by selecting genes having a coverage > 10
df_genes_of_interest <- merge(x = df_genes_of_interest,y=df_depth_over10, by = "V1",sort = FALSE,na.rm=TRUE)

sorted_bam_file_name<-paste0("sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam")
#read the file containing the match between unique.gene.name (ex.: 10043_0 ) & Gene (ex: 643886116_contig45_186)
df_uniq_vs_cell_contig<-read.csv(paste0(root_workspace,"Metadata_for_mobile_gene_dataset.csv"),stringsAsFactors = FALSE)
#only keep the information for the gene name conversion
df_uniq_vs_cell_contig <- subset(df_uniq_vs_cell_contig,select = c(Unique.gene.name,Gene,Gene.length))
df_uniq_vs_cell_contig$Gene.length <- as.integer(gsub(x=(df_uniq_vs_cell_contig$Gene.length), pattern = ",", replacement = ""))
#read the file containing the match between clusternumber/Unique.gene.name (ex.: 10043_0 ) and bam ids (ex.: gnl|BL_ORD_ID|1007602)
df_uniq_bam_id <- read.csv(paste0(root_workspace,"HGTgenes_lookup.txt"),sep="\t",stringsAsFactors = FALSE)
#get header SQ lines in data frame
df_refname_in_header_bam <- read.csv(paste0(wp_path,"headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)

#get IDs in the body of the bam file
df_refname_in_body_bam <- read.csv(paste0(wp_path,"IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
vect_refname_in_body_bam <- df_refname_in_body_bam[[1]]
df_uniq_bam_id <- merge(x=df_uniq_bam_id,y=df_refname_in_body_bam,by.x="othername",by.y="V1")
df_uniq_bam_id$gene_len <- 0L
for (i in 1:nrow(df_uniq_bam_id)){
  df_uniq_bam_id$gene_len[i] <- (subset(x = df_refname_in_header_bam,subset= V2==df_uniq_bam_id$othername[i]))$V3
}

#Merge the 2 dataframes with the condition that Clusternumber should be the same as Unique.gene.name AND same Gene length
df_uniq_celletcontig_bamsn <- merge(df_uniq_bam_id,df_uniq_vs_cell_contig,by.x = c("Clusternumber","gene_len"), by.y=c("Unique.gene.name","Gene.length"))

#Get the list of genes name in the fasta file
df_anvi_ids_vs_original_ids <- read.csv(paste0(root_workspace,"conversion_fasta_id.txt"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
df_uniq_celletcontig_bamsn <- merge(df_uniq_celletcontig_bamsn,df_anvi_ids_vs_original_ids,by.x = "Gene", by.y="V2")
vect_seq_fasta <- read.csv(paste0(root_workspace,"conversion_fasta_id.txt"),sep="\t",stringsAsFactors = FALSE,header = FALSE)[[2]]
df_uniq_celletcontig_bamsn <- subset(df_uniq_celletcontig_bamsn,subset= Gene %in% vect_seq_fasta,select=names(df_uniq_celletcontig_bamsn))

df_uniq_celletcontig_bamsn$duplicated <- duplicated(df_uniq_celletcontig_bamsn$othername)
df_uniq_celletcontig_bamsn <- subset(df_uniq_celletcontig_bamsn, subset= duplicated == FALSE,select=names(df_uniq_celletcontig_bamsn))

#create a function that returns number of possible SINGLE-SITE synonymous mutations divided by 3 for a CODON
calculate_third_of_possible_ns_codon <- function(the_codon){
  the_codon <- toupper(the_codon)
  if (nchar(the_codon)!=3){
    stop("codon length should be 3!")
  }
  possible_single_site_mutated_codons <- rep("",9)
  num_mut_codon <-1
  for (pos_codon in 1:3){
    if (substr(the_codon,start = pos_codon,stop=pos_codon)=="A"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
      
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="T"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else if (substr(the_codon,start = pos_codon,stop=pos_codon)=="C"){
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "G"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }else{#G
      mut_codon_1 <- the_codon
      substr(mut_codon_1,start = pos_codon,stop=pos_codon) <- "A"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_1
      num_mut_codon=num_mut_codon+1
      mut_codon_2 <- the_codon
      substr(mut_codon_2,start = pos_codon,stop=pos_codon) <- "T"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_2
      num_mut_codon=num_mut_codon+1
      mut_codon_3 <- the_codon
      substr(mut_codon_3,start = pos_codon,stop=pos_codon) <- "C"
      possible_single_site_mutated_codons[num_mut_codon] <- mut_codon_3
      num_mut_codon=num_mut_codon+1
    }
  }
  #count the number of synonymous mutations based on the genetic code
  nb_unique_syn_mut_codons <-0 #default initialization
  if (the_codon == "TTT") {
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTC"])
    
  } else if (the_codon == "TTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons=="TTT"])
    
  } else if (the_codon == "TTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTG","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","CTT","CTC","CTA","CTG")])
    
  } else if (the_codon == "TCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCC","TCA","TCG","AGT","AGC")])
  } else if (the_codon == "TCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCC","TCG","AGT","AGC")])
    
  } else if (the_codon == "TCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TCT","TCA","TCC","AGT","AGC")])
    
  } else if (the_codon == "TAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAC")])
    
  } else if (the_codon == "TAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TAT")])
    
  } else if (the_codon == "TGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGC")])
    
  } else if (the_codon == "TGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TGT")])
    
  } else if (the_codon == "TGG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "CTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTC","CTA","CTG")])
    
  } else if (the_codon == "CTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTA","CTG")])
    
  } else if (the_codon == "CTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTG")])
    
  } else if (the_codon == "CTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("TTA","TTG","CTT","CTC","CTA")])
    
  } else if (the_codon == "CCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCC","CCA","CCG")])
    
  } else if (the_codon == "CCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCA","CCG")])
    
    
  } else if (the_codon == "CCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCG")])
    
  } else if (the_codon == "CCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CCT","CCC","CCA")])
    
  } else if (the_codon == "CAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAC")])
    
  } else if (the_codon == "CAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAT")])
    
  } else if (the_codon == "CAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAG")])
    
  } else if (the_codon == "CAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CAA")])
    
  } else if (the_codon == "CGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGC","CGA","CGG")])
    
  } else if (the_codon == "CGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGG")])
    
  } else if (the_codon == "CGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGC","CGG")])
    
  } else if (the_codon == "CGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("CGT","CGA","CGC")])
    
  } else if (the_codon == "ATT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATA")])
    
  } else if (the_codon == "ATC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATT","ATA")])
    
  } else if (the_codon == "ATA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ATC","ATT")])
    
  } else if (the_codon == "ATG"){
    nb_unique_syn_mut_codons <- 0
    
  } else if (the_codon == "ACT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACC","ACA","ACG")])
    
    
  } else if (the_codon == "ACC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACA","ACG")])
    
  } else if (the_codon == "ACA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACG")])
    
    
  } else if (the_codon == "ACG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("ACT","ACC","ACA")])
    
  } else if (the_codon == "AAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAC")])
    
  } else if (the_codon == "AAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAT")])
    
  } else if (the_codon == "AAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAG")])
    
  } else if (the_codon == "AAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AAA")])
    
  } else if (the_codon == "AGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGC","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGT","TCT","TCC","TCA","TCG")])
    
  } else if (the_codon == "AGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGG")])
    
  } else if (the_codon == "AGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("AGA")])
    
  } else if (the_codon == "GTT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTG")])
    
  } else if (the_codon == "GTC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTT","GTA","GTG")])
    
  } else if (the_codon == "GTA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTT","GTG")])
    
  } else if (the_codon == "GTG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GTC","GTA","GTT")])
    
  } else if (the_codon == "GCT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCG")])
    
  } else if (the_codon == "GCC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCT","GCA","GCG")])
    
  } else if (the_codon == "GCA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCT","GCG")])
    
  } else if (the_codon == "GCG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GCC","GCA","GCT")])
    
  } else if (the_codon == "GAT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAC")])
    
  } else if (the_codon == "GAC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAT")])
    
  } else if (the_codon == "GAA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAG")])
    
  } else if (the_codon == "GAG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GAA")])
    
  } else if (the_codon == "GGT"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGG")])
    
  } else if (the_codon == "GGC"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGT","GGA","GGG")])
    
  } else if (the_codon == "GGA"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGT","GGG")])
    
    
  } else if (the_codon == "GGG"){
    nb_unique_syn_mut_codons <- length(possible_single_site_mutated_codons[possible_single_site_mutated_codons %in% c("GGC","GGA","GGT")])
  }
  return((nb_unique_syn_mut_codons/3))
}

#create a function that does all possible pariwise COMPARISON (not alignment) between all the sequences in a vector(so if n is the number of sequences there will be n*(n-1) comparisons between reads pair codons. The algorithm execution time is in O(n^2) because of the if condition and the 2 for loop on the same vectorof n sequences.)
#Because at each SNV site, we're looking at a single nucleotide position, the maximum of pairwise difference to consider between 2 sequences will be one diffrent nucleotide at a single position for 2 reads)
#The output is the total number of different sequences. 
calculate_nb_pwd <- function(vector_read_codon_seqs){
  output <- 0 #initialize the variable that will contain the number of pairwise differences in single nucleotide sites
  for (i in 1:length(vector_read_codon_seqs)){
    for (j in i:length(vector_read_codon_seqs)){
      if (i != j){
        if (vector_read_codon_seqs[i]!=vector_read_codon_seqs[j]){
          output <- output + 1
        }
      }
    }
  }
  return(output)
}

#Create a function to Calculate dn_ds & number of segrating sites for a certain gene, using his id. 
calculate_dn_ds_ns_gene <- function(dataframe_SNVs,gene_id) {
  df_gene_SNVs <- subset(dataframe_SNVs, subset= contig_name==gene_id)
  nb_SNV_sites <- nrow(df_gene_SNVs) #number of SNV sites in the gene
  gene_seq <- fasta_mbg_seqs[[gene_id]][1]  
  gene_length <- nchar(gene_seq)
  gene_coverage <- (df_genes_of_interest$V2[df_genes_of_interest$V1==gene_id])[1]
  average_cov_seg_sites <- mean(x = df_gene_SNVs$coverage,na.rm = TRUE)
  Nb_sm_tot <- 0 #initialize the variable that will represent the number of synonymous mutations found for this gene
  Nb_nsm_tot <- 0 #initialize the variable that will represent the number of non-synonymous mutations found for this gene
  Nb_syn_sites <- 0 #initialize the variable that will represent the number of synonymous sites found in this gene
  Nb_nsyn_sites <- 0 #initialize the variable that will represent the number of non-synonymous sites found in this gene
  nb_stop_cods <- 0 #initialize the variable that will represent the number of stop codons for one gene
  Nb_pwdiff_gene <- 0 #initialize the variable that will represent the number of pairwise differences in the gene
  Nb_stop_codons_in_gene_reads <- 0 #initialize the variable that will represent the number of stop codons found in mobile genes mapped reads sequence considerint the reading frame of the gene
  
  #go through the gene sequence and calculate Nb_syn_sites
  for (pos_in_gene in seq(from = 1,to = gene_length,by = 3)){
    current_codon_gene <- substr(x = gene_seq,start = pos_in_gene,stop=pos_in_gene+2)
    Nb_syn_sites <- Nb_syn_sites + calculate_third_of_possible_ns_codon(current_codon_gene)
  }
  #Calculate Nb_nsyn_sites from Nb_syn_sites 
  Nb_nsyn_sites <- (gene_length - Nb_syn_sites)
  #print(paste0(gene_id,";",gene_seq,";",gene_length,";",Nb_syn_sites,";",Nb_nsyn_sites))
  
  #for each SNVs site in the gene, calculate the number of synonymous & non-synonymous mutations + determine if the site is a synonymous site and/or a non-synonymous site
  for (index_SNV_site in (1:nb_SNV_sites)) {
    df_site <- df_gene_SNVs[index_SNV_site,] #currently analysed SNV informations
    #find nucleotide at current site of the real reference gene 
    refseq_site_nucl <- ""
    site_nb_occ_A <- df_site$A
    site_nb_occ_C <- df_site$C
    site_nb_occ_G <- df_site$G
    site_nb_occ_T <- df_site$T
    max_nb_occ_max_site <- max(c(site_nb_occ_A,site_nb_occ_C,site_nb_occ_G,site_nb_occ_T),na.rm = TRUE)
    if (site_nb_occ_A == max_nb_occ_max_site){
      refseq_site_nucl <- "A"
    } else if (site_nb_occ_C == max_nb_occ_max_site){
      refseq_site_nucl <- "C"
    } else if (site_nb_occ_G == max_nb_occ_max_site){
      refseq_site_nucl <- "G"
    } else if (site_nb_occ_T == max_nb_occ_max_site){
      refseq_site_nucl <- "T"
    } else{
      stop("logical error in the script : No reference sequence nucleotide can be found based on maximum nucleotide frequency at the site. Contact the Author of the script!")
    }
    
    if ((((df_site$pos_in_contig+1)-1)%%3)==0){
      #the SNV is at the first position in its codon
      original_codon <- paste0(refseq_site_nucl,substr(x = gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1),substr(x=gene_seq,start = (df_site$pos_in_contig+1)+2,stop=(df_site$pos_in_contig+1)+2))
      mutated_codons <- c(rep(paste0("A",substr(x=gene_seq,start = (df_site$pos_in_contig+1)+1,stop = (df_site$pos_in_contig+1)+1),substr(x=gene_seq, start = (df_site$pos_in_contig+1)+2,stop=(df_site$pos_in_contig+1)+2)),df_site$A),rep(paste0("C",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1),substr(x=gene_seq, start=(df_site$pos_in_contig+1)+2, stop = (df_site$pos_in_contig+1)+2)),df_site$C),rep(paste0("G",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1),substr(x=gene_seq,start = (df_site$pos_in_contig+1)+2,stop = (df_site$pos_in_contig+1)+2)),df_site$G),rep(paste0("T",substr(x=gene_seq,start = (df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1),substr(x=gene_seq,start=(df_site$pos_in_contig+1)+2,stop = (df_site$pos_in_contig+1)+2)),df_site$T)) #initialization
      all_codons <- mutated_codons
      mutated_codons <- unique(mutated_codons[mutated_codons!=original_codon])#filter the codons to know which are different from the original_codon
      
    } else if ((((df_site$pos_in_contig+1)-2)%%3)==0){
      #the SNV is at the second position in its codon
      original_codon <- paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),refseq_site_nucl,substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1))
      mutated_codons <- c(rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"A",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1)),df_site$A),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"C",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1)),df_site$C),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"G",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1)),df_site$G),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"T",substr(x=gene_seq,start=(df_site$pos_in_contig+1)+1,stop=(df_site$pos_in_contig+1)+1)),df_site$T)) #initialization
      all_codons <- mutated_codons
      mutated_codons <- unique(mutated_codons[mutated_codons!=original_codon])#filter the codons to know which are different from the original_codon
      
    } else if ((((df_site$pos_in_contig+1)-3)%%3)==0){
      #the SNV is at the third position in its codon
      original_codon <- paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-2,stop=(df_site$pos_in_contig+1)-2),substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),refseq_site_nucl)
      mutated_codons <- c(rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-2,stop=(df_site$pos_in_contig+1)-2),substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"A"),df_site$A),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-2,stop=(df_site$pos_in_contig+1)-2),substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"C"),df_site$C),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-2,stop=(df_site$pos_in_contig+1)-2),substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"G"),df_site$G),rep(paste0(substr(x=gene_seq,start=(df_site$pos_in_contig+1)-2,stop=(df_site$pos_in_contig+1)-2),substr(x=gene_seq,start=(df_site$pos_in_contig+1)-1,stop=(df_site$pos_in_contig+1)-1),"T"),df_site$T)) #initialization
      all_codons <- mutated_codons
      mutated_codons <- unique(mutated_codons[mutated_codons!=original_codon])#filter the codons to know which are different from the original_codon
      
    }
    
    if (nchar(original_codon)!=3 || any(nchar(mutated_codons)!=3)){
      print(paste("problem with codon length for gene ",df_site$contig_name, " at 0-based position ",df_site$pos, ", original codon is ", original_codon, "mutated codons are ", mutated_codons," and gene sequence is ",gene_seq))
      stop()
    }
    nb_sm_site <- 0 #initialization (IMPORTANT WHEN THERE IS A STOP CODONS SO DON'T REMOVE THIS COMMAND)
    nb_nsm_site <- 0 #initialization (IMPORTANT WHEN THERE IS A STOP CODONS SO DON'T REMOVE THIS COMMAND)
    #count the number of synonymous and non-synonymous mutations based on the genetic code
    if (original_codon == "TTT") {
      nb_sm_site <- length(mutated_codons[mutated_codons=="TTC"])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TTC"){
      nb_sm_site <- length(mutated_codons[mutated_codons=="TTT"])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TTA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTG","CTT","CTC","CTA","CTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
     
    } else if (original_codon == "TTG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTA","CTT","CTC","CTA","CTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TCT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TCC","TCA","TCG","AGT","AGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
     
    } else if (original_codon == "TCC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TCT","TCA","TCG","AGT","AGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TCA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TCT","TCC","TCG","AGT","AGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
     
    } else if (original_codon == "TCG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TCT","TCA","TCC","AGT","AGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
     
    } else if (original_codon == "TAT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TAC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TAC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TAT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TAA"){
      nb_stop_cods <- nb_stop_cods + 1
    } else if (original_codon == "TAG"){
      nb_stop_cods <- nb_stop_cods + 1
    } else if (original_codon == "TGT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TGC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TGT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "TGA"){
      nb_stop_cods <- nb_stop_cods + 1
    } else if (original_codon == "TGG"){
      nb_sm_site <- 0
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CTT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTA","TTG","CTC","CTA","CTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CTC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTA","TTG","CTT","CTA","CTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CTA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTA","TTG","CTT","CTC","CTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CTG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("TTA","TTG","CTT","CTC","CTA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CCT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CCC","CCA","CCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CCC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CCT","CCA","CCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CCA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CCT","CCC","CCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CCG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CCT","CCC","CCA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CAT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CAC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
  
    } else if (original_codon == "CAC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CAT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CAA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CAG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CAG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CAA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CGT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CGC","CGA","CGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CGC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CGT","CGA","CGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CGA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CGT","CGC","CGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "CGG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("CGT","CGA","CGC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ATT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ATC","ATA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ATC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ATT","ATA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ATA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ATC","ATT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ATG"){
      nb_sm_site <- 0
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ACT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ACC","ACA","ACG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ACC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ACT","ACA","ACG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ACA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ACT","ACC","ACG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "ACG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("ACT","ACC","ACA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AAT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AAC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AAC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AAT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AAA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AAG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AAG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AAA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AGT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AGC","TCT","TCC","TCA","TCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AGC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AGT","TCT","TCC","TCA","TCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AGA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "AGG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("AGA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GTT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GTC","GTA","GTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GTC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GTT","GTA","GTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GTA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GTC","GTT","GTG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GTG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GTC","GTA","GTT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site

      
    } else if (original_codon == "GCT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GCC","GCA","GCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site

    } else if (original_codon == "GCC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GCT","GCA","GCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site

      
    } else if (original_codon == "GCA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GCC","GCT","GCG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GCG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GCC","GCA","GCT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GAT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GAC")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GAC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GAT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GAA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GAG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GAG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GAA")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
    
    } else if (original_codon == "GGT"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GGC","GGA","GGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GGC"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GGT","GGA","GGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GGA"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GGC","GGT","GGG")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    } else if (original_codon == "GGG"){
      nb_sm_site <- length(mutated_codons[mutated_codons %in% c("GGC","GGA","GGT")])
      nb_nsm_site <- length(mutated_codons) - nb_sm_site
      
    }

    nb_pwdiff_site <- calculate_nb_pwd(all_codons)
    Nb_pwdiff_gene <- Nb_pwdiff_gene + nb_pwdiff_site
    
    Nb_nsm_tot <- Nb_nsm_tot + nb_nsm_site
    Nb_sm_tot <- Nb_sm_tot + nb_sm_site
  
    Nb_stop_codons_in_gene_reads <- Nb_stop_codons_in_gene_reads + length(all_codons[all_codons %in% c("TAA","TAG","TGA")]) 
  }
  dn <- (Nb_nsm_tot)/(Nb_nsyn_sites)
  ds <- (Nb_sm_tot)/(Nb_syn_sites)
  dn_ds <- dn/ds
  if (is.infinite(dn_ds)){
    dn_ds <- NA
  }
  if (is.nan(dn_ds)){
    dn_ds <- NA
  }
  return(list(dn_ds, nb_SNV_sites, Nb_pwdiff_gene,nb_stop_cods,Nb_stop_codons_in_gene_reads,ds,dn,Nb_nsm_tot,Nb_nsyn_sites,Nb_sm_tot,Nb_syn_sites,gene_length,df_gene_SNVs$coverage,average_cov_seg_sites))
}

#initialize result data_frame
df_SNV <- merge(x = df_SNV,y=df_genes_of_interest,by.x = "contig_name", by.y = "V1")
  #Exclude sample if not enough genes to analyse after ids conversion
if (length(unique(df_SNV$contig_name)) < 500){
  print(paste0("Sample ",the_sample_name," is excluded from analysis because it only contains ",length(unique(df_SNV$contig_name))," genes."))
  opt <- options(show.error.messages=FALSE) 
  on.exit(options(opt)) 
  stop(paste0("END OF PIPELINE EXECUTION FOR Sample ",the_sample_name," because it only contains ",length(unique(df_SNV$contig_name))," genes."))
}
df_results <- data.frame(gene_id=unique(df_SNV$contig_name),dn_ds=c(0.0),gene_coverage=c(0.0),Ne_k_hat=c(0),Ne_S_Taj=c(0),gene_length=c(0L), D_Taj=c(0.0), ds= c(0.0), dn= c(0.0), nb_pol_sites =c(0L), nsm=c(0),nss=c(0),sm=c(0),ss=c(0),stringsAsFactors = FALSE)
#define dn_ds values and number of segrating sites for each gene in the appropriate dataframe.
vect_nb_seg_sites <- rep(0,length(df_results$gene_id))
vect_tot_pw_diffs <- rep(0,length(df_results$gene_id))
mu <- (10^(-10)) #average prokaryotic mutation rate from Drake et al. 1991
rep_path_no_filter <- paste0(rep_path,"Output_when_No_filter_on_genes_and_SNVs/")
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ",rep_path_no_filter),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ",rep_path_no_filter),intern = FALSE,wait = TRUE)
rep_path_single_filter <- paste0(rep_path,"Output_single_filter/")
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ",rep_path_single_filter),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ",rep_path_single_filter),intern = FALSE,wait = TRUE)
rep_path_10_or_more_ps_filter <- paste0(rep_path_single_filter,"Output_for_genes_with_10_or_more_ps/")
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ",rep_path_10_or_more_ps_filter),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ",rep_path_10_or_more_ps_filter),intern = FALSE,wait = TRUE)
rep_path_5_or_more_occ_mut <- paste0(rep_path_single_filter,"Output_for_genes_with_5_occ_each_SNVs_mutations/")
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ",rep_path_5_or_more_occ_mut),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ",rep_path_5_or_more_occ_mut),intern = FALSE,wait = TRUE)
rep_path_double_filter <- paste0(rep_path,"Output_with_Two_filters/")
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ",rep_path_double_filter),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ",rep_path_double_filter),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("echo '' > ",rep_path_no_filter,"Analysis_Report_",the_sample_name,".txt"),intern = FALSE,wait = TRUE)
for (i in 1:length(df_results$gene_id)){
  df_results$gene_coverage[i] <- floor(df_SNV$V2[df_SNV$contig_name==df_results$gene_id[i]][1]) #define gene_coverage values
  Nb_stop_codons <- 0 #initialize the variable that will report stop codons 
  results_current_gene <- calculate_dn_ds_ns_gene(df_SNV,df_results$gene_id[i])
  df_results$dn_ds[i] <- results_current_gene[[1]]
  df_results$ds[i] <- (results_current_gene[[6]])
  df_results$dn[i] <- (results_current_gene[[7]])
  df_results$nb_pol_sites[i] <- results_current_gene[[2]]
  df_results$nsm[i] <- results_current_gene[[8]]
  df_results$nss[i] <- results_current_gene[[9]]
  df_results$sm[i] <- results_current_gene[[10]]
  df_results$ss[i] <- results_current_gene[[11]]
  vect_nb_seg_sites[i] <- results_current_gene[[2]]
  vect_tot_pw_diffs[i] <- results_current_gene[[3]]
  Nb_stop_codons <- Nb_stop_codons +  results_current_gene[[4]]
  df_results$gene_length[i] <- results_current_gene[[12]]
  
  # #get number of reads of the current gene in the sample
  # nb_reads_current_gene <- (subset(x = df_the_sample_genes_length_and_nb_reads,subset = V1==df_results$gene_id[i]))$V3[1]
  # #parameters to calculate expected sqrt_variance
  # a1_current_gene <- sum((1:(nb_reads_current_gene-1))^-1) 
  # a_1 <- a1_current_gene
  # a2_current_gene <- sum((1:(nb_reads_current_gene-1))^-2) 
  # b1_current_gene <- (nb_reads_current_gene+1)/(3*(nb_reads_current_gene-1))
  # b2_current_gene <- (2*((nb_reads_current_gene^2)+nb_reads_current_gene+3))/((9*nb_reads_current_gene)*(nb_reads_current_gene-1))
  # c1_current_gene <- b1_current_gene - (1/a1_current_gene)
  # c2_current_gene <- b2_current_gene - ((nb_reads_current_gene+2)/(a1_current_gene*nb_reads_current_gene)) + (a2_current_gene/(a1_current_gene^2))
  # e1_current_gene <- c1_current_gene/a1_current_gene
  # e2_current_gene <- c2_current_gene/((a1_current_gene^2)+a2_current_gene)
  
  #get Average coverage at segregating sites for the current gene in the current sample
  mapping_sample_size <- results_current_gene[[14]] #Average coverage at segregating sites
  #parameters to calculate expected sqrt_variance
  a1_current_gene <- sum((1:(mapping_sample_size-1))^-1) 
  a_1 <- a1_current_gene
  a2_current_gene <- sum((1:(mapping_sample_size-1))^-2) 
  b1_current_gene <- (mapping_sample_size+1)/(3*(mapping_sample_size-1))
  b2_current_gene <- (2*((mapping_sample_size^2)+mapping_sample_size+3))/((9*mapping_sample_size)*(mapping_sample_size-1))
  c1_current_gene <- b1_current_gene - (1/a1_current_gene)
  c2_current_gene <- b2_current_gene - ((mapping_sample_size+2)/(a1_current_gene*mapping_sample_size)) + (a2_current_gene/(a1_current_gene^2))
  e1_current_gene <- c1_current_gene/a1_current_gene
  e2_current_gene <- c2_current_gene/((a1_current_gene^2)+a2_current_gene)
  
  #Find S , find a_1, calculate expected sqrt_variance with formula form Tajima (1989) and calculate Tajima's D and Ne_S_Taj according to it
  no_filter_S_current_gene <- vect_nb_seg_sites[i]
  sqrt_expected_variance_current_gene <- sqrt((e1_current_gene*no_filter_S_current_gene)+((e2_current_gene*no_filter_S_current_gene)*(no_filter_S_current_gene-1)))

  df_results$Ne_S_Taj[i] <- vect_nb_seg_sites[i]/(2*mu*a_1) #based on equation 5 of Tajima(1989) "Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism"
  df_results$Ne_k_hat[i] <- vect_tot_pw_diffs[i] / (2*mu*sum(choose(k=2,n=results_current_gene[[13]]),na.rm = TRUE))  #based on equation 10 of Tajima(1989) "Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism"
  df_results$D_Taj[i] <- ((df_results$Ne_k_hat[i] - df_results$Ne_S_Taj[i]) * (2*mu))/sqrt_expected_variance_current_gene
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo '-GENE",i,"(",df_results$gene_id[i],")"," : Gene length = ", df_results$gene_length[i],"; a1 = ", a_1, "; Average gene coverage = ",df_results$gene_coverage[i],"; number of segragating sites = ",vect_nb_seg_sites[i],"; Estimated Ne_S_Taj = ",df_results$Ne_S_Taj[i],"; Estimated Ne_k_hat = ",df_results$Ne_k_hat[i], "; Estimated dn/ds = ",df_results$dn_ds[i], "' >> ",rep_path_no_filter,"Analysis_Report_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo 'Reported number of stop codons in the gene reference sequence (if there is more than one, one of them may be a premature stop codon) = ",Nb_stop_codons, "' >> ", rep_path_no_filter,"Analysis_Report_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo 'Reported number of stop codons in the sequence of the gene mapped reads (if there is more than one, one of them may be a premature stop codon) = ",results_current_gene[[5]], "' >> ", rep_path_no_filter,"Analysis_Report_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
}

#******************NO FILTER************************
#plot Ne, dn/ds & average gene coverage density distribution for sample + save plots in the repertory of the SNV matrix file

#dn/ds distribution
ggplot(data=df_results) + geom_density(mapping=aes(x = dn_ds),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes dn/ds for sample ",the_sample_name," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 1, color = "black",lty=2)
#save graph
ggsave(filename = paste0("dn_ds_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#Average gene coverage distribution
ggplot(data=df_results) + geom_density(mapping=aes(x = gene_coverage),na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes coverage for sample ",the_sample_name," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("Average gene coverage") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("gene_coverage_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot Ne_k_hat distribution
ggplot(data=df_results) + geom_density(mapping=aes(x = Ne_k_hat), color = "blue",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Ne for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989)"," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("Estimated Ne") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot Ne_S_Taj distribution
ggplot(data=df_results) + geom_density(mapping=aes(x = Ne_S_Taj), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Ne for sample ",the_sample_name," using Ne estimator created from Equation 5 of Tajima (1989)"," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("Estimated Ne") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_S_Tajima_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot log_10_(Theta_s) distribution. Theta_s is also named Theta_Watterson
ggplot(data=df_results) + geom_density(mapping=aes(x = log10(2*Ne_S_Taj*mu)), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Theta_S in a logarithmic scale for sample ",the_sample_name," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("log10(Theta_s)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Theta_s_log10_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_k_hat vs gene_coverage
ggplot(data=df_results,aes(x=log10(gene_coverage), y=Ne_k_hat)) + geom_point(color = "orange",na.rm = TRUE) + ggtitle(label = paste0("Estimated mobile genes Ne in function of average gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (No Filter)"), subtitle = paste0("n = ", nrow(df_results)," genes ")) + xlab("log10(Average_gene_coverage)") + ylab("Estimated Ne") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_S_Taj vs gene_coverage
ggplot(data=df_results,aes(x=log10(gene_coverage), y=Ne_S_Taj)) + geom_point(color = "orange",na.rm = TRUE) + ggtitle(label = paste0("Estimated mobile genes Ne in function of average gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 5 of Tajima (1989) (No Filter)"), subtitle = paste0("n = ", nrow(df_results)," genes ")) + xlab("log10(Average_gene_coverage)") + ylab("Estimated Ne") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_S_Taj_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot Ne_k_hat in function of dn/ds
ggplot(data=df_results,aes(x=dn_ds, y=Ne_k_hat)) + geom_point(color = "blue",na.rm = TRUE)+ ggtitle(label = paste0("Correlation between mobile genes Ne and dn/ds for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (No Filter)"), subtitle = paste0("n = ", nrow(df_results)," genes ")) + xlab("dn/ds") + ylab("Estimated Ne") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_VS_dn_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot D_Taj distribution
ggplot(data=df_results) + geom_density(mapping=aes(x = D_Taj), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of Tajima's D for sample ",the_sample_name," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("Tajima's D") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2)
#save graph
ggsave(filename = paste0("D_Taj_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot dn_ds vs ds to identify noise and dn_ds cutoff
ggplot(data=df_results) + geom_line(mapping=aes(x = log10(ds), y = log10(dn_ds)), color = "black",na.rm = TRUE) + ggtitle(paste0("dn/ds in function of ds in a logarithmic scale for identification of ds noise in sample ",the_sample_name," (n = ", nrow(df_results)," genes)"),subtitle = "No Filter") + xlab("log10(ds)") + ylab("log10(dn/ds)")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("ds_noise_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot dn vs ds 
ggplot(data=df_results) + geom_point(mapping=aes(x = dn, y = ds), color = "black",na.rm = TRUE) + ggtitle(paste0("dn in function of ds in sample ",the_sample_name," (data with no filters; n = ", nrow(df_results)," genes)")) + xlab("ds") + ylab("dn")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("dn_VS_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_no_filter, width = 35, height = 20, units = "cm")

#plot a barplot of the number of genes in the 3 different categories of Tajima's D values
nb_gene_D_Taj_nega <- length(df_results$gene_id[df_results$D_Taj < 0])
nb_gene_D_Taj_zero <- length(df_results$gene_id[df_results$D_Taj == 0])
nb_gene_D_Taj_posi <- length(df_results$gene_id[df_results$D_Taj > 0])

png(filename = paste0(rep_path_no_filter,"Barplot_nb_Genes_and_Tajima_D_category_for_sample_",the_sample_name,".png"))
par(mar=c(15, 3, 3, 1)) #set appropriate margins for barplot
barplot(c(nb_gene_D_Taj_nega,nb_gene_D_Taj_zero,nb_gene_D_Taj_posi),names.arg=c("nb_gene_negative_D_Taj","nb_gene_D_Taj_zero","nb_gene_positive_D_Taj"),las=2,main=paste0("Number of Genes in the each Tajima's D category for sample ",the_sample_name),ylim=range(pretty(c(0,max(c(nb_gene_D_Taj_nega,nb_gene_D_Taj_zero,nb_gene_D_Taj_posi))))))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1)) #default margin

#******************SINGLE FILTER FOR KEEPING genes with 10+ polymorphic sites************************
df_results_filter_10_or_more_ps <- subset(df_results, subset= nb_pol_sites>=10)

#dn_ds distribution on filtered data
ggplot(data=df_results_filter_10_or_more_ps) + geom_density(mapping=aes(x = dn_ds),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes dn/ds for sample ",the_sample_name," (n = ", nrow(df_results_filter_10_or_more_ps)," genes); ","(data with single filter 10+ polymorphic sites)")) + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 1, color = "black",lty=2)
#save graph
ggsave(filename = paste0("dn_ds_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot log_10_(Theta_s) distribution. Theta_s is also named Theta_Watterson
ggplot(data=df_results_filter_10_or_more_ps) + geom_density(mapping=aes(x = log10(2*Ne_S_Taj*mu)), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Theta_S (Theta_Watterson) in a logarithmic scale for sample ",the_sample_name," (data with single filter 10+ polymorphic sites)"," (n = ", nrow(df_results_filter_10_or_more_ps)," genes)")) + xlab("log10(Theta_s)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Theta_s_log10_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot dn_ds vs ds to identify noise and dn_ds cutoff
ggplot(data=df_results_filter_10_or_more_ps) + geom_line(mapping=aes(x = log10(ds), y = log10(dn_ds)), color = "black",na.rm = TRUE) + ggtitle(paste0("dn/ds in function of ds in a logarithmic scale for identification of ds noise in sample ",the_sample_name," (data with single filter 10+ polymorphic sites)"," (n = ", nrow(df_results_filter_10_or_more_ps)," genes)")) + xlab("log10(ds)") + ylab("log10(dn/ds)")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("ds_noise_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot dn vs ds 
ggplot(data=df_results_filter_10_or_more_ps) + geom_point(mapping=aes(x = dn, y = ds), color = "black",na.rm = TRUE) + ggtitle(paste0("dn in function of ds in sample ",the_sample_name," (data with single filter 10+ polymorphic sites; n = ", nrow(df_results_filter_10_or_more_ps)," genes)")) + xlab("ds") + ylab("dn")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("dn_VS_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot Ne_k_hat in function of dn/ds
ggplot(data=df_results_filter_10_or_more_ps,aes(x=log10(dn_ds), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE)+ ggtitle(label = paste0("Correlation between mobile genes Ne and dn/ds in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989)"), subtitle = paste0("n = ", nrow(df_results_filter_10_or_more_ps)," genes","; (Single Filter 10+ polymorphic sites)")) + xlab("log10(dn/ds)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_VS_dn_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_k_hat vs gene_coverage
ggplot(data=df_results_filter_10_or_more_ps,aes(x=log10(gene_coverage), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989)"), subtitle = paste0("n = ", nrow(df_results_filter_10_or_more_ps)," genes; (Single Filter 10+ polymorphic sites)")) + xlab("log10(Average gene coverage)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_S_Taj vs gene_coverage
ggplot(data=df_results_filter_10_or_more_ps,aes(x=log10(gene_coverage), y=log10(Ne_S_Taj))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 5 of Tajima (1989)"), subtitle = paste0("n = ", nrow(df_results_filter_10_or_more_ps)," genes; (Single Filter 10+ polymorphic sites)")) + xlab("log10(Average gene coverage)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_S_Taj_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#dn_ds distribution on filtered data
ggplot(data=df_results_filter_10_or_more_ps) + geom_density(mapping=aes(x = dn_ds),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes dn/ds for sample ",the_sample_name," (n = ", nrow(df_results_filter_10_or_more_ps)," genes); ","(data with single filter 10+ polymorphic sites)")) + xlab("dn/ds)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 1, color = "black",lty=2)
#save graph
ggsave(filename = paste0("dn_ds_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#Ne_k_hat distribution on filtered data
ggplot(data=df_results_filter_10_or_more_ps) + geom_density(mapping=aes(x = Ne_k_hat),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Ne for sample ",the_sample_name, " using Eq. 10 of Tajima (1989)"," (n = ", nrow(df_results_filter_10_or_more_ps)," genes; ","(data with single filter 10+ polymorphic sites)")) + xlab("Estimated Ne") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))

#save graph
ggsave(filename = paste0("Ne_k_hat_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")

#plot D_Taj distribution
ggplot(data=df_results_filter_10_or_more_ps) + geom_density(mapping=aes(x = D_Taj), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of Tajima's D for sample ",the_sample_name," (n = ", nrow(df_results_filter_10_or_more_ps)," genes; ","(data with single filter 10+ polymorphic sites)")) + xlab("Tajima's D") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2)
#save graph
ggsave(filename = paste0("D_Taj_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_10_or_more_ps_filter, width = 35, height = 20, units = "cm")


#plot a barplot of the number of genes likely evolving under specific evolution force
nb_gene_on_pos_sel <- length(df_results_filter_10_or_more_ps$gene_id[df_results_filter_10_or_more_ps$dn_ds > 1])
nb_gene_neutral <- length(df_results_filter_10_or_more_ps$gene_id[df_results_filter_10_or_more_ps$dn_ds == 1])
nb_gene_on_pur_sel <- length(df_results_filter_10_or_more_ps$gene_id[df_results_filter_10_or_more_ps$dn_ds < 1])

png(filename = paste0(rep_path_10_or_more_ps_filter,"nonNull_nb_sm_or_no_nsm_data_Barplot_Genes_and_EvolutionForce_based_on_dn_ds_sample_",the_sample_name,".png"))
par(mar=c(15, 3, 3, 1)) #set appropriate margins for barplots
barplot(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel),names.arg=c("nb_gene_on_pur_sel","nb_gene_neutral","nb_gene_on_pos_sel"),las=2,main=paste0("Number of genes under different Evolution Forces based on dn/ds interpretation of sample ",the_sample_name,"; (data with single filter 10+ polymorphic sites)"),ylim = range(pretty(c(0,max(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel))))))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#default margins

#plot a barplot of the number of genes likely evolving under specific evolution force
nb_gene_on_pos_sel <- length(df_results$gene_id[df_results$dn_ds > 1])
nb_gene_neutral <- length(df_results$gene_id[df_results$dn_ds == 1])
nb_gene_on_pur_sel <- length(df_results$gene_id[df_results$dn_ds < 1])

png(filename = paste0(rep_path_no_filter,"No_filter_data_Barplot_Genes_and_EvolutionForce_based_on_dn_ds_sample_",the_sample_name,".png"))
par(mar=c(15, 3, 3, 1)) #set appropriate margins for barplots
barplot(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel),names.arg=c("nb_gene_on_pur_sel","nb_gene_neutral","nb_gene_on_pos_sel"),las=2,main=paste0("Number of genes evolving under different Evolution Forces based on dn/ds interpretation of sample ",the_sample_name,"; (data with no filter)"),ylim = range(pretty(c(0,max(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel))))))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#default margins

#******************SINGLE FILTER FOR KEEPING gene sites with at least 5 occurences of each reads containing a mutation************************
#initialize result data_frame
df_SNV_filter_5_for_all_mut <- subset(x = df_SNV, subset = ((A>=5) + (C>=5) + (G>=5) + (T>=5) >= ((A>0) + (C>0) + (G>0) + (T>0))))
df_results_5_or_more_occ_all_mut <- data.frame(gene_id=unique(df_SNV_filter_5_for_all_mut$contig_name),dn_ds=c(0.0),gene_coverage=c(0.0),Ne_k_hat=c(0),Ne_S_Taj=c(0),gene_length=c(0L), D_Taj=c(0.0), ds=c(0.0), dn=c(0.0), nb_pol_sites =c(0L), nsm=c(0),nss=c(0),sm=c(0),ss=c(0), stringsAsFactors = FALSE)
#define dn_ds values and number of segrating sites for each gene in the appropriate dataframe.
vect_nb_seg_sites <- rep(0,length(df_results_5_or_more_occ_all_mut$gene_id))
vect_tot_pw_diffs <- rep(0,length(df_results_5_or_more_occ_all_mut$gene_id))
mu <- (10^(-10)) #average prokaryotic mutation rate from Drake et al. 1991
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("echo '' > ",rep_path_5_or_more_occ_mut,"filtered_genes_Analysis_Report_with_single_filter_10_or_more_ps_",the_sample_name,".txt"),intern = FALSE,wait = TRUE)
for (i in 1:length(df_results_5_or_more_occ_all_mut$gene_id)){
  df_results_5_or_more_occ_all_mut$gene_coverage[i] <- floor(df_SNV_filter_5_for_all_mut$V2[df_SNV_filter_5_for_all_mut$contig_name==df_results_5_or_more_occ_all_mut$gene_id[i]][1]) #define gene_coverage values
  Nb_stop_codons <- 0 #initialize the variable that will report stop codons 
  results_current_gene <- calculate_dn_ds_ns_gene(df_SNV_filter_5_for_all_mut,df_results_5_or_more_occ_all_mut$gene_id[i])
  df_results_5_or_more_occ_all_mut$dn_ds[i] <- results_current_gene[[1]]
  df_results_5_or_more_occ_all_mut$ds[i] <- (results_current_gene[[6]])
  df_results_5_or_more_occ_all_mut$dn[i] <- (results_current_gene[[7]])
  df_results_5_or_more_occ_all_mut$nb_pol_sites[i] = results_current_gene[[2]]
  df_results_5_or_more_occ_all_mut$nsm[i] <- results_current_gene[[8]]
  df_results_5_or_more_occ_all_mut$nss[i] <- results_current_gene[[9]]
  df_results_5_or_more_occ_all_mut$sm[i] <- results_current_gene[[10]]
  df_results_5_or_more_occ_all_mut$ss[i] <- results_current_gene[[11]]
  vect_nb_seg_sites[i] <- results_current_gene[[2]]
  vect_tot_pw_diffs[i] <- results_current_gene[[3]]
  Nb_stop_codons <- Nb_stop_codons +  results_current_gene[[4]]
  df_results_5_or_more_occ_all_mut$gene_length[i] <- results_current_gene[[12]]
  
  # #get number of reads of the current gene in the sample
  # nb_reads_current_gene <- (subset(x = df_the_sample_genes_length_and_nb_reads,subset = V1==df_results$gene_id[i]))$V3[1]
  # #parameters to calculate expected sqrt_variance
  # a1_current_gene <- sum((1:(nb_reads_current_gene-1))^-1) 
  # a_1 <- a1_current_gene
  # a2_current_gene <- sum((1:(nb_reads_current_gene-1))^-2) 
  # b1_current_gene <- (nb_reads_current_gene+1)/(3*(nb_reads_current_gene-1))
  # b2_current_gene <- (2*((nb_reads_current_gene^2)+nb_reads_current_gene+3))/((9*nb_reads_current_gene)*(nb_reads_current_gene-1))
  # c1_current_gene <- b1_current_gene - (1/a1_current_gene)
  # c2_current_gene <- b2_current_gene - ((nb_reads_current_gene+2)/(a1_current_gene*nb_reads_current_gene)) + (a2_current_gene/(a1_current_gene^2))
  # e1_current_gene <- c1_current_gene/a1_current_gene
  # e2_current_gene <- c2_current_gene/((a1_current_gene^2)+a2_current_gene)
  
  #get Average coverage at segregating sites for the current gene in the current sample
  mapping_sample_size <- results_current_gene[[14]] #Average coverage at segregating sites
  #parameters to calculate expected sqrt_variance
  a1_current_gene <- sum((1:(mapping_sample_size-1))^-1) 
  a_1 <- a1_current_gene
  a2_current_gene <- sum((1:(mapping_sample_size-1))^-2) 
  b1_current_gene <- (mapping_sample_size+1)/(3*(mapping_sample_size-1))
  b2_current_gene <- (2*((mapping_sample_size^2)+mapping_sample_size+3))/((9*mapping_sample_size)*(mapping_sample_size-1))
  c1_current_gene <- b1_current_gene - (1/a1_current_gene)
  c2_current_gene <- b2_current_gene - ((mapping_sample_size+2)/(a1_current_gene*mapping_sample_size)) + (a2_current_gene/(a1_current_gene^2))
  e1_current_gene <- c1_current_gene/a1_current_gene
  e2_current_gene <- c2_current_gene/((a1_current_gene^2)+a2_current_gene)
  
  #Find S , find a_1, calculate expected sqrt_variance with formula form Tajima (1989) and calculate Tajima's D and Ne_S_Taj according to it
  single_filter_S_current_gene <- vect_nb_seg_sites[i]
  sqrt_expected_variance_current_gene <- sqrt((e1_current_gene*single_filter_S_current_gene)+((e2_current_gene*single_filter_S_current_gene)*(single_filter_S_current_gene-1)))
  
  df_results_5_or_more_occ_all_mut$Ne_S_Taj[i] <- vect_nb_seg_sites[i]/(2*mu*a_1) #based on equation 5 of Tajima(1989) "Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism"
  df_results_5_or_more_occ_all_mut$Ne_k_hat[i] <- vect_tot_pw_diffs[i] / (2*mu*sum(choose(k=2,n=results_current_gene[[13]]),na.rm = TRUE))  #based on equation 10 of Tajima(1989) "Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism"
  df_results_5_or_more_occ_all_mut$D_Taj[i] <- ((df_results_5_or_more_occ_all_mut$Ne_k_hat[i] - df_results_5_or_more_occ_all_mut$Ne_S_Taj[i]) * (2*mu))/sqrt_expected_variance_current_gene
  
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo '-GENE",i,"(",df_results_5_or_more_occ_all_mut$gene_id[i],")"," : Gene length = ", df_results_5_or_more_occ_all_mut$gene_length[i],"; a1 = ", a_1, "; Average gene coverage = ",df_results_5_or_more_occ_all_mut$gene_coverage[i],"; number of segragating sites = ",vect_nb_seg_sites[i],"; Estimated Ne_S_Taj = ",df_results_5_or_more_occ_all_mut$Ne_S_Taj[i],"; Estimated Ne_k_hat = ",df_results_5_or_more_occ_all_mut$Ne_k_hat[i], "; Estimated dn/ds = ",df_results_5_or_more_occ_all_mut$dn_ds[i], "' >> ",rep_path_5_or_more_occ_mut,"filtered_genes_Analysis_Report_with_single_filter_10_or_more_ps_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo 'Reported number of stop codons in the gene reference sequence (if there is more than one, one of them may be a premature stop codon) = ",Nb_stop_codons, "' >> ", rep_path_5_or_more_occ_mut,"filtered_genes_Analysis_Report_with_single_filter_10_or_more_ps_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE, command = paste0("echo 'Reported number of stop codons in the sequence of the gene mapped reads (if there is more than one, one of them may be a premature stop codon) = ",results_current_gene[[5]], "' >> ", rep_path_5_or_more_occ_mut,"filtered_genes_Analysis_Report_with_single_filter_10_or_more_ps_",the_sample_name,".txt"),intern = FALSE, wait = TRUE)
}

#plot log_10_(Theta_s) distribution. Theta_s is also named Theta_Watterson
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_density(mapping=aes(x = log10(2*Ne_S_Taj*mu)), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Theta_S (Theta_Watterson) in a logarithmic scale for sample ",the_sample_name," (data with single filter SNVs mutations frequency >= 5)"," (n = ", nrow(df_results_5_or_more_occ_all_mut)," genes)")) + xlab("log10(Theta_s)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Theta_s_log10_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#plot dn_ds vs ds to identify noise and dn_ds cutoff
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_line(mapping=aes(x = log10(ds), y = log10(dn_ds)), color = "black",na.rm = TRUE) + ggtitle(paste0("dn/ds in function of ds in a logarithmic scale for identification of ds noise in sample ",the_sample_name," (data with single filter SNVs mutations frequency >= 5)"," (n = ", nrow(df_results_5_or_more_occ_all_mut)," genes)")) + xlab("log10(ds)") + ylab("log10(dn/ds)")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("ds_noise_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#plot dn vs ds 
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_point(mapping=aes(x = dn, y = ds), color = "black",na.rm = TRUE) + ggtitle(paste0("dn in function of ds in sample ",the_sample_name," (data with single filter SNVs mutations frequency >= 5; n = ", nrow(df_results_5_or_more_occ_all_mut)," genes)")) + xlab("ds") + ylab("dn")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("dn_VS_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#plot Ne_k_hat in function of dn/ds
ggplot(data=df_results_5_or_more_occ_all_mut,aes(x=log10(dn_ds), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and dn/ds in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (data with single filter SNVs mutations frequency >= 5)"), subtitle = paste0("n = ", nrow(df_results_5_or_more_occ_all_mut)," genes")) + xlab("log10(dn/ds)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_VS_dn_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#plot Estimation Ne_k_hat vs gene_coverage
ggplot(data=df_results_5_or_more_occ_all_mut,aes(x=log10(gene_coverage), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (data with single filter SNVs mutations frequency >= 5)"), subtitle = paste0("n = ", nrow(df_results_5_or_more_occ_all_mut)," genes")) + xlab("log10(Average gene coverage)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#dn_ds distribution on filtered data
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_density(mapping=aes(x = dn_ds),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes dn/ds for sample ",the_sample_name," (n = ", nrow(df_results_5_or_more_occ_all_mut)," genes); ","(data with single filter SNVs mutations frequency >= 5)")) + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 1, color = "black",lty=2)
#save graph
ggsave(filename = paste0("dn_ds_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#Ne_k_hat distribution on filtered data
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_density(mapping=aes(x = Ne_k_hat),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Ne for sample ",the_sample_name, " using Eq. 10 of Tajima (1989)"," (n = ", nrow(df_results_5_or_more_occ_all_mut)," genes); ","(data with single filter SNVs mutations frequency >= 5)")) + xlab("Estimated Ne") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))

#save graph
ggsave(filename = paste0("Ne_k_hat_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")

#plot D_Taj distribution
ggplot(data=df_results_5_or_more_occ_all_mut) + geom_density(mapping=aes(x = D_Taj), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of Tajima's D for sample ",the_sample_name," (n = ", nrow(df_results_5_or_more_occ_all_mut)," genes); ","(data with single filter SNVs mutations frequency >= 5)")) + xlab("Tajima's D") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2)
#save graph
ggsave(filename = paste0("D_Taj_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_5_or_more_occ_mut, width = 35, height = 20, units = "cm")


#plot a barplot of the number of genes likely evolving under specific evolution force
nb_gene_on_pos_sel <- length(df_results_5_or_more_occ_all_mut$gene_id[df_results_5_or_more_occ_all_mut$dn_ds > 1])
nb_gene_neutral <- length(df_results_5_or_more_occ_all_mut$gene_id[df_results_5_or_more_occ_all_mut$dn_ds == 1])
nb_gene_on_pur_sel <- length(df_results_5_or_more_occ_all_mut$gene_id[df_results_5_or_more_occ_all_mut$dn_ds < 1])

png(filename = paste0(rep_path_5_or_more_occ_mut,"nonNull_nb_sm_or_no_nsm_data_Barplot_Genes_and_EvolutionForce_based_on_dn_ds_sample_",the_sample_name,".png"))
par(mar=c(15, 3, 3, 1)) #set appropriate margins for barplots
barplot(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel),names.arg=c("nb_gene_on_pur_sel","nb_gene_neutral","nb_gene_on_pos_sel"),las=2,main=paste0("Number of genes under different Evolution Forces based on dn/ds interpretation of sample ",the_sample_name,"; (data with single filter SNVs mutations frequency >= 5)"),ylim = range(pretty(c(0,max(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel))))))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#default margins


#****************************************TWO FILTERS*********************************************
#Combine the 2 previous filters

df_results_two_filters <- subset(df_results_5_or_more_occ_all_mut, subset=nb_pol_sites>=10)

#plot log_10_(Theta_s) distribution. Theta_s is also named Theta_Watterson
ggplot(data=df_results_two_filters) + geom_density(mapping=aes(x = log10(2*Ne_S_Taj*mu)), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes Theta_S (Theta_Watterson) in a logarithmic scale for sample ",the_sample_name," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("log10(Theta_s)") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Theta_s_log10_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot dn_ds vs ds to identify noise and dn_ds cutoff
ggplot(data=df_results_two_filters) + geom_line(mapping=aes(x = log10(ds), y = log10(dn_ds)), color = "black",na.rm = TRUE) + ggtitle(paste0("dn/ds in function of ds in a logarithmic scale for identification of ds noise in sample ",the_sample_name," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("log10(ds)") + ylab("log10(dn/ds)")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("ds_noise_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot dn vs ds 
ggplot(data=df_results_two_filters) + geom_point(mapping=aes(x = dn, y = ds), color = "black",na.rm = TRUE) + ggtitle(paste0("dn in function of ds in sample ",the_sample_name," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("ds") + ylab("dn")  + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("dn_VS_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot Ne_k_hat in function of dn/ds
ggplot(data=df_results_two_filters,aes(x=log10(dn_ds), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE)+ ggtitle(label = paste0("Correlation between mobile genes Ne and dn/ds in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (Two Filters)"), subtitle = paste0("n = ", nrow(df_results_two_filters)," genes")) + xlab("log10(dn/ds)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_VS_dn_ds_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_k_hat vs gene_coverage
ggplot(data=df_results_two_filters,aes(x=log10(gene_coverage), y=log10(Ne_k_hat))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (Two Filters)"), subtitle = paste0("n = ", nrow(df_results_two_filters)," genes")) + xlab("log10(Average gene coverage)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_k_hat_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot Estimation Ne_S_Taj vs gene_coverage
ggplot(data=df_results_two_filters,aes(x=log10(gene_coverage), y=log10(Ne_S_Taj))) + geom_point(color = "blue",na.rm = TRUE) + ggtitle(label = paste0("Correlation between mobile genes Ne and gene coverage in a logarithmic scale for sample ",the_sample_name," using Ne estimator created from Equation 10 of Tajima (1989) (Two Filters)"), subtitle = paste0("n = ", nrow(df_results_two_filters)," genes")) + xlab("log10(Average gene coverage)") + ylab("log10(Estimated Ne)") + theme(plot.title=element_text(hjust=0,size=8))
#save graph
ggsave(filename = paste0("Ne_S_Taj_vs_gene_coverage_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#dn_ds distribution on filtered data
ggplot(data=df_results_two_filters) + geom_density(mapping=aes(x = dn_ds),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of mobile genes dn/ds for sample ",the_sample_name," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("dn/ds") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 1, color = "black",lty=2)
#save graph
ggsave(filename = paste0("dn_ds_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#Ne_S_Taj distribution on filtered data
ggplot(data=df_results_two_filters) + geom_density(mapping=aes(x = Ne_S_Taj),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of estimated mobile genes Ne for sample ",the_sample_name, " using eq. 5 of Tajima (1989)"," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("Estimated Ne") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8)) 
#save graph
ggsave(filename = paste0("Ne_S_Taj_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#Ne_k_hat distribution on filtered data
ggplot(data=df_results_two_filters) + geom_density(mapping=aes(x = Ne_k_hat),color="black",na.rm = TRUE) + ggtitle(paste0("Distribution of estimated mobile genes Ne for sample ",the_sample_name, " using eq. 10 of Tajima (1989)"," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("Estimated Ne") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))

#save graph
ggsave(filename = paste0("Ne_K_hat_estimation_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot D_Taj distribution
ggplot(data=df_results_two_filters) + geom_density(mapping=aes(x = D_Taj), color = "red",na.rm = TRUE) + ggtitle(paste0("Distribution of Tajima's D for sample ",the_sample_name," (data with two filters; n = ", nrow(df_results_two_filters)," genes)")) + xlab("Tajima's D") + ylab("Density")+ theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2)
#save graph
ggsave(filename = paste0("D_Taj_distribution_SAMPLE_",the_sample_name,".png"), path=rep_path_double_filter, width = 35, height = 20, units = "cm")

#plot a barplot of the number of genes likely evolving under specific evolution force
nb_gene_on_pos_sel <- length(df_results_two_filters$gene_id[df_results_two_filters$dn_ds > 1])
nb_gene_neutral <- length(df_results_two_filters$gene_id[df_results_two_filters$dn_ds == 1])
nb_gene_on_pur_sel <- length(df_results_two_filters$gene_id[df_results_two_filters$dn_ds < 1])

png(filename = paste0(rep_path_double_filter,"Double_filtered_data_Barplot_Genes_and_EvolutionForce_based_on_dn_ds_sample_",the_sample_name,".png"))
par(mar=c(15, 3, 3, 1)) #set appropriate margins for barplots
barplot(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel),names.arg=c("nb_gene_negative_sel","nb_gene_neutral","nb_gene_positive_sel"),las=2,main=paste0("Number of genes under different Evolution Forces based on dn/ds interpretation for double-filtered data of sample ",the_sample_name),ylim=range(pretty(c(0,max(c(nb_gene_on_pur_sel,nb_gene_neutral,nb_gene_on_pos_sel))))))
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#default margins

#add the sample name in the dataframe before saving each table
df_results$Sample <- the_sample_name
df_results_filter_10_or_more_ps$Sample <- the_sample_name
df_results_5_or_more_occ_all_mut$Sample <- the_sample_name
df_results_two_filters$Sample <- the_sample_name

#save the information about the number of genes analysed (depth over 10) & the number of genes kept after the filters
png(filename = paste0(rep_path,"Number_of_genes_kept_for_analysis_for_sample_",the_sample_name,".png"))
par(mar=c(20, 3, 3, 1)) #set appropriate margins for barplots
barplot(c(nrow(df_results),nrow(df_results_filter_10_or_more_ps),nrow(df_results_5_or_more_occ_all_mut),nrow(df_results_two_filters)),names.arg=c("total number of genes(depth over 10)","number of genes after filter 10+ ps","number of genes after filter 5+ mutation occurences","number of genes after both filters"),las=2,main = paste0("Number of genes kept for analysis depending on filters for sample ",the_sample_name), ylim=range(pretty(c(0, nrow(df_results)))),xpd = FALSE)
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))#default margins

#save df_results in a file named Data_Exploration_results_sample_#SAMPLENAME.csv and in the appropriate repertory
write.table(x = df_results, file = paste0(rep_path_no_filter,"Data_Exploration_results_sample_",the_sample_name,".csv"),row.names=FALSE, na="",col.names=TRUE, sep="\t", quote = FALSE, append = FALSE)

#save df_results_filter_10_or_more_ps in a file of the appropriate repertory
write.table(x = df_results_filter_10_or_more_ps, file = paste0(rep_path_10_or_more_ps_filter,"Data_Exploration_results_sample_",the_sample_name,"_filter_10_or_more_ps.csv"),row.names=FALSE, na="",col.names=TRUE, sep="\t", quote = FALSE, append = FALSE)

#save df_results_5_or_more_occ_all_mut in a file of the appropriate repertory
write.table(x = df_results_5_or_more_occ_all_mut, file = paste0(rep_path_5_or_more_occ_mut,"Data_Exploration_results_sample_",the_sample_name,"_filter_5_or_more_occur_mut.csv"),row.names=FALSE, na="",col.names=TRUE, sep="\t", quote = FALSE, append = FALSE)

#save df_results_two_filters in a file of the appropriate repertory
write.table(x = df_results_two_filters, file = paste0(rep_path_double_filter,"Data_Exploration_results_sample_",the_sample_name,"_two_filters_for_nb_pol_sites_and_mutation_read_freq.csv"),row.names=FALSE, na="",col.names=TRUE, sep="\t", quote = FALSE, append = FALSE)


#print short reports
print(paste0("----------Short report for sample ", the_sample_name,"----------"))
print(paste0("plots are in 'Output' directories in ", rep_path))
print(paste0("(No filter) Results file path is : '",rep_path_no_filter,"Data_Exploration_results_sample_",the_sample_name,".csv"))
print(paste0("(first Single Filter) Results file path is : '",rep_path_10_or_more_ps_filter,"Data_Exploration_results_sample_",the_sample_name,"_filter_10_or_more_ps.csv"))
print(paste0("(second Single Filter) Results file path is : '",rep_path_5_or_more_occ_mut,"Data_Exploration_results_sample_",the_sample_name,"_filter_5_or_more_occur_mut.csv"))
print(paste0("(Two filters) Results file path is : '",rep_path_double_filter,"Data_Exploration_results_sample_",the_sample_name,"_two_filters_for_nb_pol_sites_and_mutation_read_freq.csv"))
print("********************************")

system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("echo '",the_sample_name,"' >> ",root_workspace,"Appropriate_samples_list.txt"),intern = FALSE,wait = TRUE)
