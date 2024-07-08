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

			m = re.search(">(.*)",line)

			if m:
				name = m.group(1)
				
				if name not in sequences:
					sequences[name] = ""

			else:
				sequences[name] = sequences[name] + line

	return(sequences)



def count_kmers_from_seq(sequences,k,n):

	ALL_K_mers = {}

	for key in sequences:

		sequence = sequences[key]
		
		l = len(sequence)

		start = 0

		end = start + k

		while( start < l - k + 1):

			kmer = sequence[start:end]

			### count kmers 
			if kmer in ALL_K_mers:
				ALL_K_mers[kmer] += 1
			else:
				ALL_K_mers[kmer] = 1

			start = start + 1 + n
			end = start + k
	return(ALL_K_mers)

def generate_DM(sequences,ALL_K_mers,k,n):

	r = len(sequences)
	c = len(ALL_K_mers)

	DM_matrix = np.zeros((r,c),dtype=int)
	kmers = list(ALL_K_mers.keys())
	sequence_names = list(sequences.keys())

	for i in range(len(sequence_names)):
		sequence = sequences[sequence_names[i]]
		
		start = 0
		end = start + k

		l = len(sequence)

		while( start < l - k + 1):

			kmer = sequence[start:end]
			index = kmers.index(kmer)
			
			if type(index) is not list:
				j = index
				DM_matrix[i,j] += 1
			else:
				for j in index:
					DM_matrix[i,j] += 1

			start = start + 1 + n
			end = start + k
	
	## Check if DM matrix is constructed correctly
	for i in range(c):
		if ALL_K_mers[kmers[i]] != np.sum(DM_matrix[:,i]):
			sys.exit('SOMETHING WRONG.')

	presence_matrix= np.where(DM_matrix != 0, 1, 0)
	
	return(DM_matrix,presence_matrix)

		
















