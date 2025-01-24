
import os
import numpy as np
import pandas as pd
import sys
import time

#import utility scripts

import utility_functions as uf
import multiprocessing as mp
import spike_point_mass as sp_pointmass
import spike_normal as sp_normal

def main():

	DIR = os.path.realpath(os.path.dirname(__file__))

	args = uf.parse_arguments()

	if args.task == "count":

		## STEP 1: Read multi-fasta file and store them in a dict

		sequences = uf.read_fasta_file(args.sequence)

		## STEP 2: identify unique kmers from the sequences

		ALL_K_mers = uf.count_kmers_from_seq(sequences,args.k,args.gap)
		print("Finished counting unique k-mers. Identified %d kmers in total" %(len(ALL_K_mers)))

		## STEP 3: perform sequence based clustering on unique kmers (default = True)

		if args.cluster:
			kmer_clusters = uf.kmer_clustering(ALL_K_mers)

		## STEP 3: generate the kmer design matrix for each sequence, the number indicates the dosage / presence of the kmer

		sequence_names,dosage,presence = uf.generate_DM(sequences,ALL_K_mers,args.k,args.gap)

		dosage_pd = pd.DataFrame(dosage)
		presence_pd = pd.DataFrame(presence)
		dosage_pd.index = sequence_names
		presence_pd.index = sequence_names
		dosage_pd.to_csv(args.output+"_DosageMatrix.txt",header=ALL_K_mers)
		presence_pd.to_csv(args.output+"_IndicatorMatrix.txt",header=ALL_K_mers)

	elif args.task == "mapping":

		## The following script will perform the mapping algorithm to identify the causal 

		y, X, kmer_names,C = uf.read_input_files(args.geno,args.pheno,args.covar)

		trace_container = mp.Manager().dict()
		gamma_container = mp.Manager().dict()
		beta_container = mp.Manager().dict()
		alpha_container = mp.Manager().dict()

		processes = []

		if args.model == 1:
			for num in range(args.num):
				p = mp.Process(target = sp_normal.sampling,args=(args.verbose,y,C,X,args.s0,12000,args.output,num,trace_container,gamma_container,beta_container,alpha_container))
				processes.append(p)
				p.start()

		else:
			for num in range(args.num):
				p = mp.Process(target = sp_pointmass.sampling, args=(args.verbose,y,C,X,12000,args.output,num,trace_container,gamma_container,beta_container,alpha_container))
				processes.append(p)
				p.start()

		for process in processes:
			process.join()

		alpha_posterior = []
		alpha_posterior_sd = []
		beta_posterior = []
		beta_posterior_sd = []
		kmer_pip = []
		trace_posterior = []
		trace_posterior_sd = []

		for num in range(args.num):
			alpha_posterior.append(alpha_container[num]["avg"])
			alpha_posterior_sd.append(alpha_container[num]["sd"])
			beta_posterior.append(beta_container[num]["avg"])
			beta_posterior_sd.append(beta_container[num]["sd"])
			trace_posterior.append(trace_container[num]["avg"])
			trace_posterior_sd.append(trace_container[num]["sd"])
			kmer_pip.append(gamma_container[num]["kmer"])


		alpha_posterior_median = np.median(alpha_posterior,axis=0)
		alpha_posterior_sd_median = np.median(alpha_posterior_sd,axis=0)
		beta_posterior_median = np.median(beta_posterior,axis=0)
		beta_posterior_sd_median = np.median(beta_posterior_sd,axis=0)
		trace_posterior_median = np.median(trace_posterior,axis=0)
		trace_posterior_sd_median = np.median(trace_posterior_sd,axis=0)
		kmer_pip_median = np.median(kmer_pip,axis=0)


		## calculate FDR for different kmers
		index,kmer_fdr = uf.fdr_calculation(kmer_pip_median)
		
		## sort pip, kmer names, beta and beta_sd based on pip
		sorted_kmer_names = kmer_names[index]
		sorted_kmer_pip = kmer_pip_median[index]
		sorted_beta_median = beta_posterior_median[index]
		sorted_beta_sd_median = beta_posterior_sd_median[index]


		OUTPUT_KMER = open(args.output+"_kmer_pip.txt","w")
		print("%s\t%s\t%s" %("kmer_name","pip","fdr"),file = OUTPUT_KMER)
		for i in range(len(sorted_kmer_names)):
			print("%s\t%s\t%s" %(sorted_kmer_names[i],sorted_kmer_pip[i],kmer_fdr[i]),file = OUTPUT_KMER)

		OUTPUT_ALPHA = open(args.output+"_alpha.txt","w")
		print("%s\t%s" %("covariate_effect","covariate_effect_sd"),file = OUTPUT_ALPHA)
		for i in range(len(alpha_posterior_median)):
			print("%f\t%f" %(alpha_posterior_median[i],alpha_posterior_sd_median[i]),file = OUTPUT_ALPHA)

		OUTPUT_BETA = open(args.output+"_beta.txt","w")
		print("%s\t%s\t%s" %("kmer_name","kmer_effect","kmer_effect_sd"),file = OUTPUT_BETA)
		for i in range(len(sorted_kmer_names)):
			print("%s\t%f\t%f" %(sorted_kmer_names[i],sorted_beta_median[i],sorted_beta_sd_median[i]),file = OUTPUT_BETA)


		trace_variables = ["Sigma_1","Sigma_e","Large_Effect_Kmer_Proportion","Variance_Explained","Num_Causal_Kmer"]
		OUTPUT_TRACE = open(args.output+"_trace.txt","w")
		for i in range(len(trace_posterior_median)):
			print("%s\t%f\t%f" %(trace_variables[i],trace_posterior_median[i],trace_posterior_sd_median[i]),file = OUTPUT_TRACE)

	else:
		sys.exit("ERROR: Please provide the name of the task: count or mapping. Details see the manual (-h).")



## run program
if __name__ == "__main__":
    main()
