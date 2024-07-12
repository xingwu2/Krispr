
import os
import numpy as np
import pandas as pd
import sys

#import utility scripts

import utility_functions as uf
import multiprocessing as mp
import additive_gibbs as add

def main():

	DIR = os.path.realpath(os.path.dirname(__file__))

	args = uf.parse_arguments()

	if args.task == "count":

		## STEP 1: Read multi-fasta file and store them in a dict

		sequences = uf.read_fasta_file(args.sequence)

		## STEP 2: identify unique kmers from the sequences

		ALL_K_mers = uf.count_kmers_from_seq(sequences,args.k,args.gap)

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

		for num in range(args.num):
			p = mp.Process(target = add.sampling, args=(args.verbose,y,C,X,12000,args.output,num,trace_container,gamma_container,beta_container,alpha_container))
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

		OUTPUT_KMER = open(args.output+"_kmer_pip.txt","w")
		for i in range(len(kmer_pip_median)):
			print("%s\t%s" %(kmer_names[i],kmer_pip_median[i]),file = OUTPUT_KMER)

		OUTPUT_ALPHA = open(args.output+"_alpha.txt","w")
		for i in range(len(alpha_posterior_median)):
			print("%f\t%f" %(alpha_posterior_median[i],alpha_posterior_sd_median[i]),file = OUTPUT_ALPHA)

		OUTPUT_BETA = open(args.output+"_beta.txt","w")
		for i in range(len(beta_posterior_median)):
			print("%s\t%f\t%f" %(kmer_names[i],beta_posterior_median[i],beta_posterior_sd_median[i]),file = OUTPUT_BETA)

		OUTPUT_TRACE = open(args.output+"_trace.txt","w")
		for i in range(len(trace_posterior_median)):
			print("%f\t%f" %(trace_posterior_median[i],trace_posterior_sd_median[i]),file = OUTPUT_TRACE)

	else:
		sys.exit("ERROR: Please provide the name of the task: count or mapping. Details see the manual (-h).")



## run program
if __name__ == "__main__":
    main()
