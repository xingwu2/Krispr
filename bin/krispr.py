
import os
import numpy as np

#import utility scripts

import utility_functions as uf

def main():

	DIR = os.path.realpath(os.path.dirname(__file__))

	args = uf.parse_arguments()


	## STEP 1: Read multi-fasta file and store them in a dict

	sequences = uf.read_fasta_file(args.sequence)

	## STEP 2: count kmer from the sequences

	ALL_K_mers = uf.count_kmers_from_seq(sequences,args.k,args.gap)

	## STEP 3: generate the kmer design matrix for each sequence, the number indicates the dosage / presence of the kmer

	dosage,presence = uf.generate_DM(sequences,ALL_K_mers,args.k,args.gap)

	np.savetxt(args.output+"_DosageMatrix.txt",dosage,delimiter="\t",fmt='%i')
	np.savetxt(args.output+"_PresenceMatrix.txt",presence,delimiter="\t",fmt='%i')



## run program
if __name__ == "__main__":
    main()
