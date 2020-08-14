#@Author=Arnaud NG
#Script arguments :
#workspace_path : the absolute path of the workspace where there are the samples directories
#sdn_list_filepath : the path of file containing the list of the sample directories

#load libraries
library("ggplot2")
library("nlme")
library("ggdendro")
library("reshape2")
library("grid")
library("RColorBrewer")
library("cluster")
library("gplots")
library("FD")

workspace_path <- as.character(commandArgs(TRUE)[1])
sdn_list_filepath <- as.character(commandArgs(TRUE)[2])

#import sample names in vector
v_samples_directories <- (read.csv2(file=sdn_list_filepath,sep='\t',header = FALSE, stringsAsFactors=FALSE))$V1 #load the list of the samples directories
#delete Output folder if it already exists
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/rm -rf ", workspace_path,"Output_Samples_comparisons"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Samples_size"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Tajima_D"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Tajima_D/No_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Tajima_D/Single_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Tajima_D/Double_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/dn_ds"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Ne"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Ne/No_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Ne/Single_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Ne/Double_filter"),intern = FALSE,wait = TRUE)
system(ignore.stdout = FALSE, ignore.stderr = TRUE,command=paste0("/bin/mkdir ", workspace_path,"Output_Samples_comparisons/Comparison_Matrixes"),intern = FALSE,wait = TRUE)


#Build the 3 big dataframes containing all samples results (no filter + single filter for excluding genes with undefined dn_ds due to zero synonymous mutation and/or zero non-synonymous mutation + double filter which is the previous filter and a filter for having results calculated without considering SNVs with mutations observed less than 5 times)
no_filter_df_results <- read.csv(paste0(workspace_path,v_samples_directories[1],"/profile_subset_",v_samples_directories[1],"/Output_when_No_filter_on_genes_and_SNVs/Data_Exploration_results_sample_",v_samples_directories[1],".csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
no_filter_df_results$gene_id <- as.character(no_filter_df_results$gene_id)
no_filter_df_results$Sample <- as.character(no_filter_df_results$Sample)
single_filter_df_results <- read.csv(paste0(workspace_path,v_samples_directories[1],"/profile_subset_",v_samples_directories[1],"/Output_single_filter/Output_for_genes_with_5_occ_each_SNVs_mutations/Data_Exploration_results_sample_",v_samples_directories[1],"_filter_5_or_more_occur_mut.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
single_filter_df_results$gene_id <- as.character(single_filter_df_results$gene_id)
single_filter_df_results$Sample <- as.character(single_filter_df_results$Sample)
double_filter_df_results <- read.csv(paste0(workspace_path,v_samples_directories[1],"/profile_subset_",v_samples_directories[1],"/Output_with_Two_filters/Data_Exploration_results_sample_",v_samples_directories[1],"_two_filters_for_nb_pol_sites_and_mutation_read_freq.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
double_filter_df_results$gene_id <- as.character(double_filter_df_results$gene_id)
double_filter_df_results$Sample <- as.character(double_filter_df_results$Sample)
for (i in 2:length(v_samples_directories)){
  nf_current_sample_df_results <- read.csv(paste0(workspace_path,v_samples_directories[i],"/profile_subset_",v_samples_directories[i],"/Output_when_No_filter_on_genes_and_SNVs/Data_Exploration_results_sample_",v_samples_directories[i],".csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
  sf_current_sample_df_results <- read.csv(paste0(workspace_path,v_samples_directories[i],"/profile_subset_",v_samples_directories[i],"/Output_single_filter/Output_for_genes_with_5_occ_each_SNVs_mutations/Data_Exploration_results_sample_",v_samples_directories[i],"_filter_5_or_more_occur_mut.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
  Df_current_sample_df_results <- read.csv(paste0(workspace_path,v_samples_directories[i],"/profile_subset_",v_samples_directories[i],"/Output_with_Two_filters/Data_Exploration_results_sample_",v_samples_directories[i],"_two_filters_for_nb_pol_sites_and_mutation_read_freq.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE)
  no_filter_df_results <- rbind(no_filter_df_results,nf_current_sample_df_results)
  no_filter_df_results$gene_id <- as.character(no_filter_df_results$gene_id)
  no_filter_df_results$Sample <- as.character(no_filter_df_results$Sample)
  single_filter_df_results <- rbind(single_filter_df_results,sf_current_sample_df_results)
  single_filter_df_results$gene_id <- as.character(single_filter_df_results$gene_id)
  single_filter_df_results$Sample <- as.character(single_filter_df_results$Sample)
  double_filter_df_results <- rbind(double_filter_df_results,Df_current_sample_df_results)
  double_filter_df_results$gene_id <- as.character(double_filter_df_results$gene_id)
  double_filter_df_results$Sample <- as.character(double_filter_df_results$Sample)
}
  
  #No Filter Data comparisons

#plot D_Taj comparisons
ggplot(data=no_filter_df_results,aes(x = D_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Tajima's D distribution between samples (No Filter)")) + xlab("Tajima's D") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_D_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Tajima_D/No_filter"), width = 15, height = 10, units = "cm")

#plot Ne_k_hat comparisons
ggplot(data=no_filter_df_results,aes(x = Ne_k_hat,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (No Filter)")) + xlab("Ne_S_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_Ne_k_hat_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/No_filter"), width = 15, height = 10, units = "cm")

#plot Ne_k_hat comparisons (log10 SCALE)
ggplot(data=no_filter_df_results,aes(x = log10(Ne_k_hat+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (No Filter;log10 scale)")) + xlab("log10(Ne_S_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_Ne_k_hat_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/No_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons
ggplot(data=no_filter_df_results,aes(x = Ne_S_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (No Filter)")) + xlab("Ne_k_hat_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_Ne_S_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/No_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons (log10 SCALE)
ggplot(data=no_filter_df_results,aes(x = log10(Ne_S_Taj+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (No Filter; log10 scale)")) + xlab("log10(Ne_k_hat_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_Ne_S_Taj_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/No_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons 
ggplot(data=no_filter_df_results,aes(x = dn_ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (No filter) ")) + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 1, color = "black",lty=2) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_dn_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons (log10 SCALE)
ggplot(data=no_filter_df_results,aes(x = log10(dn_ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (No filter;log10 scale) ")) + xlab("log10(dn/ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_dn_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons 
ggplot(data=no_filter_df_results,aes(x = dn,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (No filter) ")) + xlab("dn") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))  + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_dn_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons (log10 SCALE)
ggplot(data=no_filter_df_results,aes(x = log10(dn+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (No filter;log10 scale) ")) + xlab("log10(dn + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_dn_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons 
ggplot(data=no_filter_df_results,aes(x = ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (No filter) ")) + xlab("ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none") 
#save graph
ggsave(filename = paste0("no_filter_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons (log10 SCALE)
ggplot(data=no_filter_df_results,aes(x = log10(ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (No filter;log10 scale) ")) + xlab("log10(ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("no_filter_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/No_filter"), width = 15, height = 10, units = "cm")


  #Single Filter Data comparisons

#plot D_Taj comparisons
ggplot(data=single_filter_df_results,aes(x = D_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of Tajima's D distribution between samples (filter on mutation frequency 5+)")) + xlab("Tajima's D") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_D_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Tajima_D/Single_filter"), width = 15, height = 10, units = "cm")

#plot Ne_k_hat comparisons
ggplot(data=single_filter_df_results,aes(x = Ne_k_hat,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (filter on mutation frequency 5+)")) + xlab("Ne_S_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_Ne_k_hat_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Single_filter"), width = 15, height = 10, units = "cm")

#plot Ne_k_hat comparisons (log10 SCALE)
ggplot(data=single_filter_df_results,aes(x = log10(Ne_k_hat+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (filter on mutation frequency 5+;log10 scale)")) + xlab("log10(Ne_S_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_Ne_k_hat_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Single_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons
ggplot(data=single_filter_df_results,aes(x = Ne_S_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (filter on mutation frequency 5+)")) + xlab("Ne_k_hat_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_Ne_S_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Single_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons (log10 SCALE)
ggplot(data=single_filter_df_results,aes(x = log10(Ne_S_Taj+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (filter on mutation frequency 5+; log10 scale)")) + xlab("log10(Ne_k_hat_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_Ne_S_Taj_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Single_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons 
ggplot(data=single_filter_df_results,aes(x = dn_ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (filter on mutation frequency 5+) ")) + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 1, color = "black",lty=2) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_dn_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons (log10 SCALE)
ggplot(data=single_filter_df_results,aes(x = log10(dn_ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (filter on mutation frequency 5+;log10 scale) ")) + xlab("log10(dn/ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_dn_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons 
ggplot(data=single_filter_df_results,aes(x = dn,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (filter on mutation frequency 5+) ")) + xlab("dn") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))  + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_dn_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons (log10 SCALE)
ggplot(data=single_filter_df_results,aes(x = log10(dn+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (filter on mutation frequency 5+;log10 scale) ")) + xlab("log10(dn + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_dn_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons 
ggplot(data=single_filter_df_results,aes(x = ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (filter on mutation frequency 5+) ")) + xlab("ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))  + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons (log10 SCALE)
ggplot(data=single_filter_df_results,aes(x = log10(ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (filter on mutation frequency 5+;log10 scale) ")) + xlab("log10(ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("single_filter_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Single_filter"), width = 15, height = 10, units = "cm")

  #Double Filter Data comparisons

#plot D_Taj comparisons
ggplot(data=double_filter_df_results,aes(x = D_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of Tajima's D distribution between samples (Two Filters)")) + xlab("Tajima's D") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))+ geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_D_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Tajima_D/Double_filter"), width = 15, height = 10, units = "cm")


#plot Ne_k_hat comparisons
ggplot(data=double_filter_df_results,aes(x = Ne_k_hat,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (Two filters)")) + xlab("Ne_S_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_Ne_k_hat_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Double_filter"), width = 15, height = 10, units = "cm")

#plot Ne_k_hat comparisons (log10 SCALE)
ggplot(data=double_filter_df_results,aes(x = log10(Ne_k_hat+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 10 from Tajima 1989) distribution between samples (Two filters;log10 scale)")) + xlab("log10(Ne_S_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_Ne_k_hat_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Double_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons
ggplot(data=double_filter_df_results,aes(x = Ne_S_Taj,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (Two filters)")) + xlab("Ne_k_hat_Tajima") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_Ne_S_Taj_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Double_filter"), width = 15, height = 10, units = "cm")

#plot Ne_S_Taj comparisons (log10 SCALE)
ggplot(data=double_filter_df_results,aes(x = log10(Ne_S_Taj+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison of mobile genes Ne (Eq. 5 from Tajima 1989) distribution between samples (Two filters; log10 scale)")) + xlab("log10(Ne_k_hat_Tajima + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(3,13)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_Ne_S_Taj_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/Ne/Double_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons 
ggplot(data=double_filter_df_results,aes(x = dn_ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (Two filters) ")) + xlab("dn/ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 1, color = "black",lty=2) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_dn_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")

#plot dn_ds distribution comparisons (log10 SCALE)
ggplot(data=double_filter_df_results,aes(x = log10(dn_ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn/ds distribution of each sample (Two filters;log10 scale) ")) + xlab("log10(dn/ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + geom_vline(xintercept = 0, color = "black",lty=2) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_dn_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons 
ggplot(data=double_filter_df_results,aes(x = dn,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (Two filters) ")) + xlab("dn") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))  + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_dn_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")

#plot dn distribution comparisons (log10 SCALE)
ggplot(data=double_filter_df_results,aes(x = log10(dn+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes dn distribution of each sample (Two filters;log10 scale) ")) + xlab("log10(dn + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_dn_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons 
ggplot(data=double_filter_df_results,aes(x = ds,fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (Two filters) ")) + xlab("ds") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8))  + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_ds_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")

#plot ds distribution comparisons (log10 SCALE)
ggplot(data=double_filter_df_results,aes(x = log10(ds+(1E-16)),fill=Sample,colour=Sample)) + geom_density(alpha=0.1,na.rm = TRUE) + ggtitle(paste0("Comparison between mobile genes ds distribution of each sample (Two filters;log10 scale) ")) + xlab("log10(ds + 1E-16)") + ylab("Density") + theme(plot.title=element_text(hjust=0,size=8)) + xlim(c(-5,5)) + theme(legend.position="none")
#save graph
ggsave(filename = paste0("double_filter_ds_log10_comparisons",".png"), path=paste0(workspace_path,"Output_Samples_comparisons/dn_ds/Double_filter"), width = 15, height = 10, units = "cm")


#Plot a barplot of the number of genes found in each sample with and without filter(s)
  #no filter
df_nb_genes_no_filter <- data.frame(Sample=no_filter_df_results$Sample,count=1)
df_nb_genes_no_filter <- data.frame(Sample=unique(df_nb_genes_no_filter$Sample), count=gapply(object = df_nb_genes_no_filter,which = "count",FUN = colSums,groups = as.factor(df_nb_genes_no_filter$Sample)))
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Samples_size/Nb_genes_in_samples_no_filter.png"))
barplot(df_nb_genes_no_filter$count,names.arg=df_nb_genes_no_filter$Sample,main="Number of Genes from each sample considered for the analysis without filter",las=2,ylim=range(pretty(c(0,max(df_nb_genes_no_filter$count)))))
dev.off()
  #single filter (exclude genes for which each mutations at the different site are not seen at least 5 times in their site)
df_nb_genes_single_filter <- data.frame(Sample=single_filter_df_results$Sample,count=1)
df_nb_genes_single_filter <- data.frame(Sample=unique(df_nb_genes_single_filter$Sample), count=gapply(object = df_nb_genes_single_filter,which = "count",FUN = colSums,groups = as.factor(df_nb_genes_single_filter$Sample)))
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Samples_size/Nb_genes_in_samples_after_single_filter.png"))
barplot(df_nb_genes_single_filter$count,names.arg=df_nb_genes_single_filter$Sample,main="Number of Genes from each sample considered for the analysis with first filter",las=2,ylim=range(pretty(c(0,max(df_nb_genes_single_filter$count)))))
dev.off()
  #double filter (exclude genes for which each mutations at the different site are not seen at least 5 times in their site + exclude SNVs for which mutations are observed less than 5 times)
df_nb_genes_double_filter <- data.frame(Sample=double_filter_df_results$Sample,count=1)
df_nb_genes_double_filter <- data.frame(Sample=unique(df_nb_genes_double_filter$Sample), count=gapply(object = df_nb_genes_double_filter,which = "count",FUN = colSums,groups = as.factor(df_nb_genes_double_filter$Sample)))
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Samples_size/Nb_genes_in_samples_after_two_filters.png"))
barplot(df_nb_genes_double_filter$count,names.arg=df_nb_genes_double_filter$Sample,main="Number of Genes from each sample considered for the analysis with two filters",las=2,ylim=range(pretty(c(0,max(df_nb_genes_double_filter$count)))))
dev.off()

#Create AND SAVE the 18 comparison matrixes 

mu <- (10^(-10)) #average prokaryotic mutation rate from Drake et al. 1991

#Tajima's D, dn/ds, coverage, Theta_S_Taj, Theta_pi, ns, dn and ds no filter Comparison Matrixes
Row_names <- unique(no_filter_df_results$gene_id)
Col_names <- unique(no_filter_df_results$Sample)
no_filter_CM_D_Taj <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_Theta_S_Taj <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_Theta_pi <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_dn_ds <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_coverage <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names)) #average gene coverage
no_filter_CM_ns <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_dn <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
no_filter_CM_ds <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
for (i in 1:length(Row_names)){
  for (j in 1:length(Col_names)){
    subset_current_gene_sample_result <- subset(x=no_filter_df_results,subset= (gene_id==Row_names[i]) & (Sample==Col_names[j]))
    if (nrow(subset_current_gene_sample_result) == 0){
      no_filter_CM_D_Taj[i,j] <- NA
      no_filter_CM_Theta_S_Taj[i,j] <- NA
      no_filter_CM_Theta_pi[i,j] <- NA
      no_filter_CM_dn_ds[i,j] <- NA
      no_filter_CM_coverage[i,j] <- NA
      no_filter_CM_ns[i,j] <- NA
      no_filter_CM_dn[i,j] <- NA
      no_filter_CM_ds[i,j] <- NA
    } else {
      no_filter_CM_D_Taj[i,j] <- subset_current_gene_sample_result$D_Taj[1]
      no_filter_CM_Theta_S_Taj[i,j] <- 2* (subset_current_gene_sample_result$Ne_S_Taj[1]) * mu
      no_filter_CM_Theta_pi[i,j] <- 2 * (subset_current_gene_sample_result$Ne_k_hat[1]) * mu
      no_filter_CM_dn_ds[i,j] <- subset_current_gene_sample_result$dn_ds[1]
      no_filter_CM_coverage[i,j] <- subset_current_gene_sample_result$gene_coverage[1]
      no_filter_CM_ns[i,j] <- subset_current_gene_sample_result$nsm[1]
      no_filter_CM_dn[i,j] <- subset_current_gene_sample_result$dn[1]
      no_filter_CM_ds[i,j] <- subset_current_gene_sample_result$ds[1]
    }
  }
}
colnames(no_filter_CM_D_Taj) <- Col_names
rownames(no_filter_CM_D_Taj) <- Row_names
write.table(no_filter_CM_D_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_D_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_Theta_S_Taj) <- Col_names
rownames(no_filter_CM_Theta_S_Taj) <- Row_names
write.table(no_filter_CM_Theta_S_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_S_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_Theta_pi) <- Col_names
rownames(no_filter_CM_Theta_pi) <- Row_names
write.table(no_filter_CM_Theta_pi, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_pi.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_dn_ds) <- Col_names
rownames(no_filter_CM_dn_ds) <- Row_names
write.table(no_filter_CM_dn_ds, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_dn_ds.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_coverage) <- Col_names
rownames(no_filter_CM_coverage) <- Row_names
write.table(no_filter_CM_coverage, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_coverage.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_ns) <- Col_names
rownames(no_filter_CM_ns) <- Row_names
write.table(no_filter_CM_ns, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_ns.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_dn) <- Col_names
rownames(no_filter_CM_dn) <- Row_names
write.table(no_filter_CM_dn, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_dn.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(no_filter_CM_ds) <- Col_names
rownames(no_filter_CM_ds) <- Row_names
write.table(no_filter_CM_ds, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_ds.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")


#Tajima's D, dn/ds, Ne_k_hat, Theta_S_Taj and Theta_pi (number of pairwise differences) SINGLE filter Comparison Matrixes
single_filter_CM_D_Taj <- no_filter_CM_D_Taj[unique(single_filter_df_results$gene_id),]
write.table(single_filter_CM_D_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_D_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
single_filter_CM_Theta_S_Taj <- no_filter_CM_Theta_S_Taj[unique(single_filter_df_results$gene_id),]
write.table(single_filter_CM_Theta_S_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_S_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
single_filter_CM_Theta_pi <- no_filter_CM_Theta_pi[unique(single_filter_df_results$gene_id),]
write.table(single_filter_CM_Theta_pi, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_pi.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
single_filter_CM_Ne_k_hat <- single_filter_CM_Theta_pi/(2*mu)
write.table(single_filter_CM_Ne_k_hat, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Ne_k_hat.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
Row_names <- unique(single_filter_df_results$gene_id)
Col_names <- unique(single_filter_df_results$Sample)
single_filter_CM_dn_ds <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
for (i in 1:length(Row_names)){
  for (j in 1:length(Col_names)){
    subset_current_gene_sample_result <- subset(x=single_filter_df_results,subset= (gene_id==Row_names[i]) & (Sample==Col_names[j]))
    if (nrow(subset_current_gene_sample_result) == 0){
      single_filter_CM_dn_ds[i,j] <- NA
    } else {
      single_filter_CM_dn_ds[i,j] <- subset_current_gene_sample_result$dn_ds[1]
    }
  }
}
colnames(single_filter_CM_dn_ds) <- Col_names
rownames(single_filter_CM_dn_ds) <- Row_names
write.table(single_filter_CM_dn_ds, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_dn_ds.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

#Tajima's D, dn/ds, Ne_k_hat, Theta_S_Taj and Theta_pi (number of pairwise differences) DOUBLE filter Comparison Matrixes
Row_names <- unique(double_filter_df_results$gene_id)
Col_names <- unique(double_filter_df_results$Sample)
double_filter_CM_D_Taj <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
double_filter_CM_Theta_S_Taj <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
double_filter_CM_Theta_pi <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
double_filter_CM_dn_ds <- matrix(c(0.0),nrow = length(Row_names), ncol=length(Col_names))
for (i in 1:length(Row_names)){
  for (j in 1:length(Col_names)){
    subset_current_gene_sample_result <- subset(x=double_filter_df_results,subset= (gene_id==Row_names[i]) & (Sample==Col_names[j]))
    if (nrow(subset_current_gene_sample_result) == 0){
      double_filter_CM_D_Taj[i,j] <- NA
      double_filter_CM_Theta_S_Taj[i,j] <- NA
      double_filter_CM_Theta_pi[i,j] <- NA
      double_filter_CM_dn_ds[i,j] <- NA
    } else {
      double_filter_CM_D_Taj[i,j] <- subset_current_gene_sample_result$D_Taj[1]
      double_filter_CM_Theta_S_Taj[i,j] <- 2* (subset_current_gene_sample_result$Ne_S_Taj[1]) * mu
      double_filter_CM_Theta_pi[i,j] <- 2 * (subset_current_gene_sample_result$Ne_k_hat[1]) * mu
      double_filter_CM_dn_ds[i,j] <- subset_current_gene_sample_result$dn_ds[1]
    }
  }
}

colnames(double_filter_CM_D_Taj) <- Col_names
rownames(double_filter_CM_D_Taj) <- Row_names
write.table(double_filter_CM_D_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_D_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(double_filter_CM_Theta_S_Taj) <- Col_names
rownames(double_filter_CM_Theta_S_Taj) <- Row_names
write.table(double_filter_CM_Theta_S_Taj, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_S_Taj.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(double_filter_CM_Theta_pi) <- Col_names
rownames(double_filter_CM_Theta_pi) <- Row_names
write.table(double_filter_CM_Theta_pi, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_pi.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
colnames(double_filter_CM_dn_ds) <- Col_names
rownames(double_filter_CM_dn_ds) <- Row_names
write.table(double_filter_CM_dn_ds, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_dn_ds.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")
double_filter_CM_Ne_k_hat <- double_filter_CM_Theta_pi/(2*mu)
write.table(double_filter_CM_Ne_k_hat, file = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Ne_k_hat.csv"),row.names=TRUE, na="NA",col.names=TRUE, sep="\t")

#print short reports
print(paste0("----------Samples Comparison report  ----------"))
print(paste0("Compared Samples are : ",unique(no_filter_df_results$Sample)))
print(paste0("plots are saved in the directories of ", workspace_path,"Output_Samples_comparisons/"))
print(paste0("Comparison matrixes are saved at ", workspace_path,"Output_Samples_comparisons/Comparison_Matrixes"))
print("********************************")

#Create Heatmaps for the appropriate comparison matrixes

#Heatmap no filter Tajima's D
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_D_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_D_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*min(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap for Tajima's D (no filter)",trace="none",scale="none", col = colorRampPalette(c("red","darkmagenta","blue","black"), space = "rgb")(24))
dev.off()

#Heatmap no filter dn/ds
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_dn_ds.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_dn_ds_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap dn/ds (no filter)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap no filter coverage
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_coverage.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
current_cmp_mtrx <- log10(current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_coverage_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "log(Coverage) (no filter)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap no filter Theta Watterson
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_S_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_S_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_S (no filter)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap no filter Theta pi
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_pi.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/no_filter_CM_Theta_pi_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_pi (no filter)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap single filter Tajima's D
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_D_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_D_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*min(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap Tajima's D (f : 5+ occurences for each mutation)",trace="none",scale="none", col = colorRampPalette(c("red","darkmagenta","blue","black"), space = "rgb")(24))
dev.off()

#Heatmap single filter Theta Watterson
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_S_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_S_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_S (f : 5+ occurences for each mutation)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap single filter Theta Pi
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_pi.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Theta_pi_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_pi (f : 5+ occurences for each mutation)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap single filter dn/ds
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_dn_ds.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_dn_ds_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap dn/ds (f : 5+ occurences for each mutation)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap single filter Ne_k_hat
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Ne_k_hat.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
current_cmp_mtrx <- log10(current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/single_filter_CM_Ne_k_hat_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "log(Ne) (f : 5+ occurences for each mutation)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(12))
dev.off()

#Heatmap two filters Tajima's D
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_D_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_D_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -10*min(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap Tajima's D (2 filters)",trace="none",scale="none", col = colorRampPalette(c("red","darkmagenta","blue","black"), space = "rgb")(12))
dev.off()

#Heatmap two filters Theta Watterson
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_S_Taj.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_S_Taj_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_S (2 filters)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap two filters Theta Pi
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_pi.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Theta_pi_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Theta_pi (2 filters)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap two filters dn/ds
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_dn_ds.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_dn_ds_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "Heatmap dn/ds (2 filters)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(24))
dev.off()

#Heatmap two filters Ne_k_hat
df_current_cmp_mtrx <- read.csv(paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Ne_k_hat.csv"),header = TRUE,sep='\t',stringsAsFactors = FALSE,na.strings = "NA")
current_cmp_mtrx <- as.matrix(x = df_current_cmp_mtrx)
current_cmp_mtrx <- log10(current_cmp_mtrx)
png(filename = paste0(workspace_path,"Output_Samples_comparisons/Comparison_Matrixes/double_filter_CM_Ne_k_hat_Heatmap.png"))
current_cmp_mtrx[is.na(current_cmp_mtrx)] <- -1*max(current_cmp_mtrx,na.rm = TRUE)
heatmap.2(x = current_cmp_mtrx, main = "log(Ne) (2 filters)",trace="none",scale="none", col = colorRampPalette(c("black","white","blue","red"), space = "rgb")(12))
dev.off()

#To zoom on heatmap, identify zone of interest in complete heatmap. Then, assign heatmap to a variable. Extract labels orders from heatmaps (column indexes are in variable$colInd and row indexes are in variable$rowInd). Read the matrix and save it in a variable. Reorder colums and lines based on heatmap labels and finally build a new heatmap with the region of interest (guess the index) 
