#@Author=Arnaud NG
#This script converts gene IDs in a sample .bam file into anvio-simplified fasta gene IDs, calls mobile genes SNVs from sample metagenomic reads alignment and execute curve_Cov_dnds_Ne.r, which computes popGen metrics for each mobile gene for which SNVs were called
#SIMPLIFY GENE ID IN FASTA WITH anvi-script-reformat-fasta    ###to do before or with system(...)
#Script arguments :
#sample_workspace_path : the absolute path of the workspace where there are the bam file, the fasta file and the files for gene IDs conversion
  #SAMPLE worspace name should be the name of the sample (ex: G30479)
#num_threads : Number of threads to allow to anvio
#anvio_bin : the absolute path of the directory where there are anvio execution files

#library for curve_Cov_dnds_Ne.r calling
library("ggplot2")
library("seqinr")
library("nlme")

#script arguments
workspace_path <- as.character(commandArgs(TRUE)[1]) #the SAMPLE workspace path is given in argument
num_threads <- as.character(commandArgs(TRUE)[2])
anvio_bin <- as.character(commandArgs(TRUE)[3])

#obtain sample name from absolute path of the workspace
indexes <- gregexpr(pattern = "/",text = workspace_path)[[1]]
the_sample_name <- substr(x=workspace_path,start=indexes[length(indexes)-1] + 1, stop=nchar(workspace_path) - 1)
parent_workspace <- substr(x=workspace_path,start=1, stop=indexes[length(indexes) - 1])
sorted_bam_file_name<-paste0("sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam")

#sort and index sample bam file
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools sort -o ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam ",workspace_path,the_sample_name,".hgt.mapped.merged.99.50.bam"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam ", workspace_path,sorted_bam_file_name),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools index -b ",workspace_path,sorted_bam_file_name," ", workspace_path,sorted_bam_file_name,".bai"),intern = FALSE,wait = TRUE)

#read the file containing the match between unique.gene.name (ex.: 10043_0 ) & Gene (ex: 643886116_contig45_186)
df_uniq_vs_cell_contig<-read.csv(paste0(parent_workspace,"Metadata_for_mobile_gene_dataset.csv"),stringsAsFactors = FALSE)
#only keep the information for the gene name conversion
df_uniq_vs_cell_contig <- subset(df_uniq_vs_cell_contig,select = c(Unique.gene.name,Gene,Gene.length))
df_uniq_vs_cell_contig$Gene.length <- as.integer(gsub(x=(df_uniq_vs_cell_contig$Gene.length), pattern = ",", replacement = ""))
#read the file containing the match between clusternumber/Unique.gene.name (ex.: 10043_0 ) and bam ids (ex.: gnl|BL_ORD_ID|1007602)
df_uniq_bam_id <- read.csv(paste0(parent_workspace,"HGTgenes_lookup.txt"),sep="\t",stringsAsFactors = FALSE)

#save header SQ lines in data frame
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools view -H ",workspace_path, sorted_bam_file_name," | grep '@SQ' | sed 's/SN://g'| sed 's/LN://g' > ", workspace_path,"headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("sort -u ",workspace_path,"headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"," > ",workspace_path,"tmp_headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("mv -f ",workspace_path,"tmp_headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv ", workspace_path,"headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
df_refname_in_header_bam <- read.csv(paste0(workspace_path,"headers_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)

#get IDs in the body of the bam file
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools view ",workspace_path, sorted_bam_file_name," | awk -F '\t' '{print $3}' > ",workspace_path,"IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("sort -u ",workspace_path,"IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv > ", workspace_path,"tmp_IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("mv -f ",workspace_path,"tmp_IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv ",workspace_path,"IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),intern = FALSE,wait = TRUE)
df_refname_in_body_bam <- read.csv(paste0(workspace_path,"IDs_body_",the_sample_name,".hgt.mapped.merged.99.50.bam.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
vect_refname_in_body_bam <- df_refname_in_body_bam[[1]]
df_uniq_bam_id <- merge(x=df_uniq_bam_id,y=df_refname_in_body_bam,by.x="othername",by.y="V1")
df_uniq_bam_id$gene_len <- 0L
for (i in 1:nrow(df_uniq_bam_id)){
  df_uniq_bam_id$gene_len[i] <- (subset(x = df_refname_in_header_bam,subset= V2==df_uniq_bam_id$othername[i]))$V3
}


#Merge the 2 dataframes with the condition that Clusternumber should be the same as Unique.gene.name AND same Gene length
df_uniq_celletcontig_bamsn <- merge(df_uniq_bam_id,df_uniq_vs_cell_contig,by.x = c("Clusternumber","gene_len"), by.y=c("Unique.gene.name","Gene.length"))

#Get the list of genes name in the fasta file
df_anvi_ids_vs_original_ids <- read.csv(paste0(parent_workspace,"conversion_fasta_id.txt"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
df_uniq_celletcontig_bamsn <- merge(df_uniq_celletcontig_bamsn,df_anvi_ids_vs_original_ids,by.x = "Gene", by.y="V2")
vect_seq_fasta <- read.csv(paste0(parent_workspace,"conversion_fasta_id.txt"),sep="\t",stringsAsFactors = FALSE,header = FALSE)[[2]]
df_uniq_celletcontig_bamsn <- subset(df_uniq_celletcontig_bamsn,subset= Gene %in% vect_seq_fasta,select=names(df_uniq_celletcontig_bamsn))

df_uniq_celletcontig_bamsn$duplicated <- duplicated(df_uniq_celletcontig_bamsn$othername)
df_uniq_celletcontig_bamsn <- subset(df_uniq_celletcontig_bamsn, subset= duplicated == FALSE,select=names(df_uniq_celletcontig_bamsn))

#import sorted bam body in dataframe AND filter the genes that can be converted into fasta file IDs
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools view ",workspace_path,sorted_bam_file_name," > ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
df_bam <- read.csv(paste0(workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
df_bam <-subset(x = df_bam,select = c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13))
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
for (i in 1:nrow(df_uniq_celletcontig_bamsn)){
  df_bam$V3<-replace(df_bam$V3,df_bam$V3==df_uniq_celletcontig_bamsn$othername[i],df_uniq_celletcontig_bamsn$V1[i])
}

df_bam <- subset(x = df_bam, subset=V10!="")
df_bam$V3 <- gsub(pattern = "gnl|BL_ORD_ID|",replacement = "gnl_BL_ORD_ID_",x=df_bam$V3,fixed=TRUE)
df_bam$V7 <- gsub(pattern = "gnl|BL_ORD_ID|",replacement = "gnl_BL_ORD_ID_",x=df_bam$V7,fixed=TRUE)

df_genes_of_interest <- data.frame(unique(df_bam$V3[grepl(x = df_bam$V3, pattern = "c_", fixed=TRUE)  & df_bam$V3 %in% df_anvi_ids_vs_original_ids$V1]))


#remove from the header file the line containing gene IDs that are not in the bam file
vect_refname_in_header_bam <- df_refname_in_header_bam$V2
for (i in 1:nrow(df_uniq_celletcontig_bamsn)){
  vect_refname_in_header_bam<-replace(vect_refname_in_header_bam,vect_refname_in_header_bam==df_uniq_celletcontig_bamsn$othername[i],df_uniq_celletcontig_bamsn$V1[i])
}
vect_refname_in_header_bam <- gsub(pattern = "gnl|BL_ORD_ID|",replacement = "gnl_BL_ORD_ID_",x=vect_refname_in_header_bam,fixed=TRUE)
vect_absent_in_body_bam <- unique(vect_refname_in_header_bam[!vect_refname_in_header_bam %in% df_bam$V3])
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools view -H ",workspace_path,sorted_bam_file_name," | sed 's/gnl|BL_ORD_ID|/gnl_BL_ORD_ID_/g' > ",workspace_path,"final_headers_",sorted_bam_file_name,".csv"),intern = FALSE,wait = TRUE)
for (i in 1:length(vect_absent_in_body_bam)){
  system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("cat ",workspace_path,"final_headers_",sorted_bam_file_name,".csv"," | grep -v '",vect_absent_in_body_bam[i],"' > ",workspace_path,"tmp_final_headers_",sorted_bam_file_name,".csv"),intern = FALSE,wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"tmp_final_headers_",sorted_bam_file_name,".csv ",workspace_path,"final_headers_",sorted_bam_file_name,".csv"),intern = FALSE,wait = TRUE)
}



# Write files 
write.table(df_genes_of_interest, file = paste0(workspace_path,"gene_IDs_of_interest_",the_sample_name,".txt"),row.names=FALSE, na="",col.names=FALSE, sep="\t")
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("sed 's/\"//g' ",workspace_path,"gene_IDs_of_interest_",the_sample_name,".txt > ",workspace_path,"tmp_gene_IDs_of_interest_",the_sample_name,".txt"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"tmp_gene_IDs_of_interest_",the_sample_name,".txt ",workspace_path,"gene_IDs_of_interest_",the_sample_name,".txt"),intern = FALSE,wait = TRUE)

system(ignore.stdout = FALSE, ignore.stderr = TRUE, command=paste0("/bin/rm ",workspace_path,sorted_bam_file_name),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE, command=paste0("/bin/rm ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/cp ",workspace_path,"final_headers_",sorted_bam_file_name,".csv ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
#replace gene IDs in the header of the file to fasta IDs 
for (i in 1:nrow(df_uniq_celletcontig_bamsn)){
  system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("cat ", workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam | sed 's/",gsub(pattern = "gnl|BL_ORD_ID|",replacement = "gnl_BL_ORD_ID_",x=df_uniq_celletcontig_bamsn$othername[i],fixed=TRUE),"/",df_uniq_celletcontig_bamsn$V1[i],"/g' > ",workspace_path,"tmp_sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
  system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"tmp_sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
}

write.table(df_bam, file = paste0(workspace_path,"tmp_body.csv"),row.names=FALSE, na="",col.names=FALSE, sep="\t")
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("cat ",workspace_path,"tmp_body.csv | sed 's/\"//g' >> ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm ",workspace_path,"tmp_body.csv"),intern = FALSE,wait = TRUE)

system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools view -bS ", workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam > ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam ", workspace_path,sorted_bam_file_name),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE, command=paste0("/bin/rm ",workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.sam"),intern = FALSE,wait = TRUE)


####Re-sort and re-index bam file avec samtools and execute SNV profile
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools sort -o ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam ",workspace_path,sorted_bam_file_name),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"p_sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam ", workspace_path,sorted_bam_file_name),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools index -b ",workspace_path,sorted_bam_file_name," ", workspace_path,sorted_bam_file_name,".bai"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("source ", anvio_bin,"activate; ",anvio_bin, "anvi-profile -i ",workspace_path,sorted_bam_file_name," -c ",parent_workspace,"contigs.db -o ",workspace_path,"profile_subset_",the_sample_name," --sample-name subset_",the_sample_name," --min-coverage-for-variability 10 --min-contig-length 50 --report-variability-full --skip-hierarchical-clustering --contigs-of-interest ",workspace_path,"gene_IDs_of_interest_",the_sample_name,".txt -T ",num_threads," ;"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("source ", anvio_bin,"activate; ",anvio_bin, "anvi-script-add-default-collection -p ",workspace_path,"profile_subset_",the_sample_name,"/PROFILE.db ;"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("source ", anvio_bin,"activate; ",anvio_bin, "anvi-gen-variability-profile -c ",parent_workspace,"contigs.db -p ",workspace_path,"profile_subset_",the_sample_name,"/PROFILE.db -C DEFAULT -b EVERYTHING --include-contig-names -o ",workspace_path,"profile_subset_",the_sample_name,"/SNVs_matrix_subset_",the_sample_name,".csv --compute-gene-coverage-stats --engine NT ;"),intern = FALSE,wait = TRUE)

#execute SAMTOOLS DEPTH ANALYSIS AND PRODUCE FILE depth.csv in workspace!!!!
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("samtools depth -a ", workspace_path,"sorted_",the_sample_name,".hgt.mapped.merged.99.50.bam > ",workspace_path,"depth.csv"),intern = FALSE,wait = TRUE)

#get the contigs having a significant coverage
df_depth <- read.csv(paste0(workspace_path,"depth.csv"),header = FALSE,sep='\t',stringsAsFactors = FALSE)
df_depth_over10 <- data.frame(V1=levels(as.factor(df_depth$V1)))
df_depth_over10$V1 <- as.character(df_depth_over10$V1)
library("nlme")
df_depth_over10$average_cov <- gapply(object = df_depth,which = "V3",FUN = colMeans,groups = as.factor(df_depth$V1))
df_depth_over10 <- subset(x = df_depth_over10, subset = average_cov > 10.0)
write.table(df_depth_over10, file = paste0(workspace_path,"genes_with_depth_over10_",the_sample_name,".csv"),row.names=FALSE, na="",col.names=FALSE, sep="\t")
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("sed 's/\"//g' ",workspace_path,"genes_with_depth_over10_",the_sample_name,".csv > ",workspace_path,"tmp_genes_with_depth_over10_",the_sample_name,".csv"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mv ",workspace_path,"tmp_genes_with_depth_over10_",the_sample_name,".csv ",workspace_path,"genes_with_depth_over10_",the_sample_name,".csv"),intern = FALSE,wait = TRUE)

system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("Rscript ",parent_workspace,"curve_Cov_dnds_Ne_.r ",workspace_path,"profile_subset_",the_sample_name,"/SNVs_matrix_subset_",the_sample_name,".csv ", the_sample_name),intern = FALSE,wait = TRUE)
