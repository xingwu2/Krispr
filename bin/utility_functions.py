import re
import argparse
import numpy as np
import sys


def parse_arguments():
	"""
    Parse command line arguments.
    """
	parser = argparse.ArgumentParser()
	parser.add_argument('-f',type = str, action= 'store',dest='sequence',help='the multi-fasta file')
	parser.add_argument('-k',type = int, action= 'store',dest='k',default=7)
	parser.add_argument('-g',type = int, action = 'store', dest = 'gap',default=0)
	parser.add_argument('-n',type = int, action = 'store', dest = 'threads',default=1)
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

		
















