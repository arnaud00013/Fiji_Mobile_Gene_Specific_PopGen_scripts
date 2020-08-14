#!PYTHON_PATH

#BEFORE EXECUTING THIS SCRIPT, make sure to activate anvio environment and the following dependencies 
    #nixpkgs/16.09 gcc/5.4.0 gsl prodigal/2.6.3 hmmer/3.1b2 samtools/1.3.1 hdf5/1.8.18 r/3.4.0 python/3.6

import sys
import time
import multiprocessing as mp
import os


#time feedback
start_time = time.time()

#A function that execute the anvio profiling step of a certain sample for which data are stored in wp_path/current_sample/
def lancer_Rscript_Curve(wp_path,current_sample, nb_t_profiling,anv_bin_path):
    current_sample_rep = "{0}{1}/".format(wp_path,current_sample)
    os.system("Rscript {0}subsample_bam_v_df.r {1} {2} {3}".format(wp_path,current_sample_rep, nb_t_profiling,anv_bin_path))

#create a function that test if a string is numeric
def is_numeric(the_str):
    try:
        float_conv = float(the_str)
        return True
    except (ValueError, TypeError):
        return False

if __name__ == "__main__":
    #Initializing global variables
    #Path of the Workspace where there are all samples directories
    workspace_path ="/home/p1211536/scratch/Cedar_FijiComp_Workspace/"
    #making sure to have "/" at the end of the workspace absolute path
    if ((workspace_path[len(workspace_path)-1]) != "/"):
        workspace_path = "{0}/".format(workspace_path)
    print("Chosen Workspace is : '{0}'".format(workspace_path))

    #Path of the anvio bin directory
    anvio_bin_path = "/home/p1211536/projects/def-shapiro/p1211536/virtual-envs/anvio5/bin/"
    #making sure to have "/" at the end of the anvio_bin absolute path
    if ((anvio_bin_path[len(anvio_bin_path)-1]) != "/"):
        anvio_bin_path = "{0}/".format(anvio_bin_path)
    print("Chosen anvio_bin is : '{0}'".format(anvio_bin_path))

    #Verify that workspace architecture is respected
    #create the variable worskpace_files_list that will contain the list of workspace files
    os.system("ls -1 {0} > {0}Workspace_file_list.txt".format(workspace_path))
    worskpace_files_list=""
    with open("{0}Workspace_file_list.txt".format(workspace_path)) as f:
        worskpace_files_list = f.readlines()


    index_line = 0
    for line in  worskpace_files_list:
        worskpace_files_list[index_line] = line.replace("\n", "")
        index_line = index_line + 1

    if (not (("HGTgenes_lookup.txt" in worskpace_files_list) and ("subsample_bam_v_df.r" in worskpace_files_list) and ("Metadata_for_mobile_gene_dataset.csv" in worskpace_files_list) and ("conversion_fasta_id.txt" in worskpace_files_list) and ("contigs.db" in worskpace_files_list) and ("Samples_comparisons.r" in worskpace_files_list) and ("curve_Cov_dnds_Ne_.r" in worskpace_files_list) and ("Data_Exploration.py" in worskpace_files_list) and ("contigs.fa" in worskpace_files_list))):
        sys.exit("ERROR : Workspace Architecture is not respected! The program execution will stop. Make sure you setup the Workspace directory properly. See the program README Annex for more details. After making the appropriate adjustments, run the program again! ")

    # create the variable samples_repertory_list that will contain the list of the samples data repertory in the workspace
    os.system("ls -1 {0} | grep 'G' | grep -v 'gene' | grep -v 'HGT' > {0}Samples_workspace.txt".format(workspace_path))
    samples_repertory_list=""
    with open("{0}Samples_workspace.txt".format(workspace_path)) as f:
        samples_repertory_list = f.readlines()

    #Number of samples being analysed simultaneously
    if (len(samples_repertory_list) % 10 == 0):
        nb_sim_process = 10
    elif (len(samples_repertory_list) % 9 == 0):
        nb_sim_process = 9
    elif (len(samples_repertory_list) % 8 == 0):
        if (len(samples_repertory_list) > 23):
            nb_sim_process = 24
        else:
            nb_sim_process = 8
    elif (len(samples_repertory_list) % 7 == 0):
        nb_sim_process = 7
    elif (len(samples_repertory_list) % 6 == 0):
        nb_sim_process = 6
    elif (len(samples_repertory_list) % 5 == 0):
        nb_sim_process = 5
    elif (len(samples_repertory_list) % 4 == 0):
        nb_sim_process = 4
    elif (len(samples_repertory_list) % 3 == 0):
        nb_sim_process = 3
    elif (len(samples_repertory_list) % 2 == 0):
        nb_sim_process = 2
    else:
        nb_sim_process = 1

    #Number of processors allowed
    nb_proc_total = 240

    #Number of threads to use for SNV profiling of 1 sample
    nb_threads_profiling = (nb_proc_total // nb_sim_process)

    #for each sample name, check if repertory contains the required files and that the sample name is valid
    index_line = 0
    for line in samples_repertory_list:
        samples_repertory_list[index_line] = line.replace("\n", "")
        current_sample_repertory= "{0}{1}/".format(workspace_path,samples_repertory_list[index_line])

        # create the variable workspace_files_list that will contain the list of the files in the current sample directory.
        os.system("ls -1 {0} > {0}Workspace_filesname.txt".format(current_sample_repertory))
        with open("{0}Workspace_filesname.txt".format(current_sample_repertory)) as f:
            current_sample_directory_files_list = f.readlines()
        indexx_line = 0
        for line in current_sample_directory_files_list:
            current_sample_directory_files_list[indexx_line] = line.replace("\n", "")
            # Make sure that the required repertory architecture is respected. Stop the program if it isn't respected.
            if (not("{0}.hgt.mapped.merged.99.50.bam".format(samples_repertory_list[index_line]) in current_sample_directory_files_list)):
                sys.exit("ERROR : Non-sorted or sorted_and_indexed bam files are missing in the Workspace {0}. The program execution will stop. Make sure you setup the Workspace directory properly. See the program README Annex for more details. After making the appropriate adjustments, run the program again! ".format(current_sample_repertory))

            indexx_line = indexx_line + 1

        index_line = index_line + 1
    
    #Create a file named "Appropriate_samples_list.txt" that will contain the list of acceptable sample to analyse, .i.e sample with enough genes annotated in HGTgenes_lookup.txt. Also create a file named "Excluded_Samples.txt" that will contain the excluded samples. Example G30796 will be excluded because it only has 14 genes that can be converted into genes fasta ids...
    if (not "Appropriate_samples_list.txt" in worskpace_files_list):
        os.system("touch {0}Appropriate_samples_list.txt".format(workspace_path))
    if (not "Excluded_Samples.txt" in worskpace_files_list):
        os.system("touch {0}Excluded_Samples.txt".format(workspace_path))

    #execute nb_sim_process samples at a time
    liste_index_sample = range(0,len(samples_repertory_list),nb_sim_process)
    for i in liste_index_sample:
        # create a list that will record the samples that have already been analysed
        appropriate_samples_list=""
        with open("{0}Appropriate_samples_list.txt".format(workspace_path)) as f:
            appropriate_samples_list = f.readlines()
        index_line = 0
        for line in appropriate_samples_list:
            appropriate_samples_list[index_line] = line.replace("\n", "")
            index_line = index_line + 1
        
        processes = []
        #define all 8 processes that will be executed parallelly (run only those that have not been executed yet)
    
        if (((i)<len(samples_repertory_list)) and (not (samples_repertory_list[i] in appropriate_samples_list))):
        
            p1 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i], nb_threads_profiling, anvio_bin_path,))
            processes.append(p1)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p1.start()
        if (((i+1)<len(samples_repertory_list)) and (not (samples_repertory_list[i+1] in appropriate_samples_list))):
            p2 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+1], nb_threads_profiling, anvio_bin_path,))
            processes.append(p2)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p2.start()
        if (((i+2)<len(samples_repertory_list)) and (not (samples_repertory_list[i+2] in appropriate_samples_list))):
            p3 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+2], nb_threads_profiling, anvio_bin_path,))
            processes.append(p3)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p3.start()
        if (((i+3)<len(samples_repertory_list)) and (not (samples_repertory_list[i+3] in appropriate_samples_list))):
            p4 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+3], nb_threads_profiling, anvio_bin_path,))
            processes.append(p4)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p4.start()
        if (((i+4)<len(samples_repertory_list)) and (not (samples_repertory_list[i+4] in appropriate_samples_list))):
            p5 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+4], nb_threads_profiling, anvio_bin_path,))
            processes.append(p5)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p5.start()
        if (((i+5)<len(samples_repertory_list)) and (not (samples_repertory_list[i+5] in appropriate_samples_list))):
            p6 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+5], nb_threads_profiling, anvio_bin_path,))
            processes.append(p6)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p6.start()
        if (((i+6)<len(samples_repertory_list)) and (not (samples_repertory_list[i+6] in appropriate_samples_list))):
            p7 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+6], nb_threads_profiling, anvio_bin_path,))
            processes.append(p7)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p7.start()
        if (((i+7)<len(samples_repertory_list)) and (not (samples_repertory_list[i+7] in appropriate_samples_list))):
            p8 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+7], nb_threads_profiling, anvio_bin_path,))
            processes.append(p8)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p8.start()

        if (((i+8)<len(samples_repertory_list)) and (not (samples_repertory_list[i+8] in appropriate_samples_list))):
            p9 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+8], nb_threads_profiling, anvio_bin_path,))
            processes.append(p9)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p9.start()
        if (((i+9)<len(samples_repertory_list)) and (not (samples_repertory_list[i+9] in appropriate_samples_list))):
            p10 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+9], nb_threads_profiling, anvio_bin_path,))
            processes.append(p10)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p10.start()
        if (((i+10)<len(samples_repertory_list)) and (not (samples_repertory_list[i+10] in appropriate_samples_list))):
            p11 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+10], nb_threads_profiling, anvio_bin_path,))
            processes.append(p11)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p11.start()
        if (((i+11)<len(samples_repertory_list)) and (not (samples_repertory_list[i+11] in appropriate_samples_list))):
            p12 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+11], nb_threads_profiling, anvio_bin_path,))
            processes.append(p12)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p12.start()
        if (((i+12)<len(samples_repertory_list)) and (not (samples_repertory_list[i+12] in appropriate_samples_list))):
            p13 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+12], nb_threads_profiling, anvio_bin_path,))
            processes.append(p13)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p13.start()
        if (((i+13)<len(samples_repertory_list)) and (not (samples_repertory_list[i+13] in appropriate_samples_list))):
            p14 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+13], nb_threads_profiling, anvio_bin_path,))
            processes.append(p14)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p14.start()
        if (((i+14)<len(samples_repertory_list)) and (not (samples_repertory_list[i+14] in appropriate_samples_list))):
            p15 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+14], nb_threads_profiling, anvio_bin_path,))
            processes.append(p15)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p15.start()
        if (((i+15)<len(samples_repertory_list)) and (not (samples_repertory_list[i+15] in appropriate_samples_list))):
            p16 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+15], nb_threads_profiling, anvio_bin_path,))
            processes.append(p16)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p16.start()
        if (((i+16)<len(samples_repertory_list)) and (not (samples_repertory_list[i+16] in appropriate_samples_list))):
            p17 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+16], nb_threads_profiling, anvio_bin_path,))
            processes.append(p17)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p17.start()
        if (((i+17)<len(samples_repertory_list)) and (not (samples_repertory_list[i+17] in appropriate_samples_list))):
            p18 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+17], nb_threads_profiling, anvio_bin_path,))
            processes.append(p18)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p18.start()
        if (((i+18)<len(samples_repertory_list)) and (not (samples_repertory_list[i+18] in appropriate_samples_list))):
            p19 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+18], nb_threads_profiling, anvio_bin_path,))
            processes.append(p19)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p19.start()
        if (((i+19)<len(samples_repertory_list)) and (not (samples_repertory_list[i+19] in appropriate_samples_list))):
            p20 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+19], nb_threads_profiling, anvio_bin_path,))
            processes.append(p20)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p20.start()
        if (((i+20)<len(samples_repertory_list)) and (not (samples_repertory_list[i+20] in appropriate_samples_list))):
            p21 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+20], nb_threads_profiling, anvio_bin_path,))
            processes.append(p21)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p21.start()
        if (((i+21)<len(samples_repertory_list)) and (not (samples_repertory_list[i+21] in appropriate_samples_list))):
            p22 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+21], nb_threads_profiling, anvio_bin_path,))
            processes.append(p22)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p22.start()
        if (((i+22)<len(samples_repertory_list)) and (not (samples_repertory_list[i+22] in appropriate_samples_list))):
            p23 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+22], nb_threads_profiling, anvio_bin_path,))
            processes.append(p23)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p23.start()
        if (((i+23)<len(samples_repertory_list)) and (not (samples_repertory_list[i+23] in appropriate_samples_list))):
            p24 = mp.Process(target=lancer_Rscript_Curve, args=(workspace_path, samples_repertory_list[i+23], nb_threads_profiling, anvio_bin_path,))
            processes.append(p24)
            # subsample_bam_v_df.r lance curve_Cov_dnds_Ne.r aussi!!!
            p24.start()

        #wait for all processes to finish
        if (len(processes) != 0):
            for p in processes:
                p.join()
            for p in processes:
                while (p.is_alive()):
                    time.sleep(5)

    os.system("Rscript {0}Samples_comparisons.r {0} {0}Appropriate_samples_list.txt".format(workspace_path))
    os.system("grep -vFxf {0}Appropriate_samples_list.txt {0}Samples_workspace.txt >> {0}Excluded_Samples.txt".format(workspace_path)) 
    print("--- Pipeline Execution started at %s ---" % (start_time))
    print("--- Pipeline Execution finished at %s ---" % (time.time()))
    print("--- Number of samples analysed : %s ---" % (len(samples_repertory_list)))
    print("--- Runtime for the pipeline with this machine is {0} seconds using {1} processors only for the anvio profiling step!!! ---".format(time.time() - start_time, nb_threads_profiling))

