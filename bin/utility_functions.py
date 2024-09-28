import re
import argparse
import numpy as np
import sys
import pandas as pd



def parse_arguments():
	"""
    Parse command line arguments.
    """
	parser = argparse.ArgumentParser()
	parser.add_argument('-f',type = str, action= 'store',dest='sequence',help='the multi-fasta file')
	parser.add_argument('-k',type = int, action= 'store',dest='k',default=7,help = "size of the kmer")
	parser.add_argument('-g',type = int, action = 'store', dest = 'gap',default=0,help = "the number of nucleotide gap between 2 kmers")
	parser.add_argument('-t',type = str, action = 'store', dest = 'task',help = "count | mapping")
	parser.add_argument('-x',type = str, action = 'store', dest = 'geno',help = "the input matrix (X) for the mapping step")
	parser.add_argument('-c',type = str, action = 'store', dest = 'covar',help = "the covariates (C) for the mapping step")
	parser.add_argument('-y',type = str, action = 'store', dest = 'pheno',help = "the response variable for the mapping step")
	parser.add_argument('-m',type = int, action = 'store', dest = 'model',default = 1, help = "the statistical model for kmer effect estimation. Krispr offers two spike priors 1 (default): small effect around 0; 2: point mass at 0")
	parser.add_argument('-s0',type = float, action = 'store', dest = 's0',default = 0.01, help = "value for sigma 0")
	parser.add_argument('-n',type = int, action = 'store', default = 8, dest = "num",help = "the number of threads for the MCMC chains. Recommend at least 5")
	parser.add_argument('-v',action = 'store_true', dest = 'verbose',default = False, help = "print out each MCMC iteration")
	parser.add_argument('-o',type = str, action = 'store', dest = 'output',help = "the prefix of the output files")
	args = parser.parse_args()

	return(args)

def read_fasta_file(file):

	sequences = {}

	with open(file,"r") as FILE:
		for line in FILE:

			line = line.strip("\n")

			## search for > for the header

			if line.startswith(">"):
				name = line[1:]
				if name not in sequences:
					sequences[name] = ""
				else:
					print("There are multiple %s sequences." %(name))
					sys.exit("ERROR: There are duplicated names in your fasta file. Please double check! ")

			else:
				if name is None:
					sys.exit("ERROR: The fasta file format is incorrect. No header line found before sequence.")
				sequences[name] += line
	print("Finished loading sequences.")
	return(sequences)



def count_kmers_from_seq(sequences,k,n):

	ALL_K_mers = []
	it = 1

	for key in sequences:

		sequence = sequences[key]
		
		l = len(sequence)

		start = 0

		end = start + k

		while( start < l - k + 1):

			kmer = sequence[start:end]

			### identify unique kmers 
			if kmer not in ALL_K_mers:
				ALL_K_mers.append(kmer)

			start = start + 1 + n
			end = start + k
		it += 1

		if it %100 ==0 :
			print("Processed %i sequences, found %i unique kmers so far" %(it,len(ALL_K_mers)))
	return(ALL_K_mers)

def generate_DM(sequences,ALL_K_mers,k,n):

	r = len(sequences)
	c = len(ALL_K_mers)

	DM_matrix = np.zeros((r,c),dtype=int)
	sequence_names = list(sequences.keys())

	for i in range(r):
		sequence = sequences[sequence_names[i]]

		for j in range(c):
			DM_matrix[i,j] = sequence.count(ALL_K_mers[j])

	presence_matrix= np.where(DM_matrix != 0, 1, 0)

	print("Finished counting unique kmer dosage for all sequences.")

	return(sequence_names,DM_matrix,presence_matrix)

def read_input_files(geno,pheno,covar):

	X = pd.read_csv(str(geno),sep="\t")
	n,p = X.shape
	kmer_names = np.array(X.columns.values.tolist())

	y = []
	with open(str(pheno),"r") as f:
		for line in f:
			line = line.strip("\n")
			y.append(float(line))

	y = np.asarray(y)

	if covar is None:
		C = np.ones(n)
		C = C.reshape(n, 1)
	else:
		C =  np.array(pd.read_csv(str(covar),sep="\t",header=None)) 

	return(y,X,kmer_names,C)

def fdr_calculation(kmer_pip_median):

	ordered_index = np.argsort(kmer_pip_median)[::-1]

	sorted_kmer_pip = kmer_pip_median[ordered_index]

	## FDR calculation from sorted pip values

	fdr = []

	for i in range(len(sorted_kmer_pip)):
		if i == 0:
			fdr.append( 1 - sorted_kmer_pip[i])
		else:
			fdr.append( 1 - np.mean(sorted_kmer_pip[:i+1]))

	return(ordered_index,fdr)



		
















