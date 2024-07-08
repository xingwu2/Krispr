
import os
import numpy as np
import pandas as pd

#import utility scripts

import utility_functions as uf

def main():

	DIR = os.path.realpath(os.path.dirname(__file__))

	args = uf.parse_arguments()


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


## run program
if __name__ == "__main__":
    main()
