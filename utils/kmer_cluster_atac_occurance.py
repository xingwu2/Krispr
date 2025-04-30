

### This python script will count how many times a kmer cluster is located under a atac peak


import argparse
import re
import pandas as pd
import numpy as np
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-i',action= 'store',dest='atac',help = 'atac peak fasta sequence')
parser.add_argument('-k',action= 'store',dest='uniquekmer',help = 'unique kmer fasta file')
parser.add_argument('-x',action= 'store',dest='dm',help = 'kmer cluster dosage matrix')
parser.add_argument('-c',action= 'store',dest='cluster',help = 'the kmer clustering output (.clstr)')
parser.add_argument('-o',action= 'store',dest='output',help='the output prefix')

args = parser.parse_args()

## load the unique kmer sequences

kmer_sequences = {}

with open(args.uniquekmer,"r") as FILE:
	for line in FILE:

		line = line.strip("\n")

			## search for > for the header

		if line.startswith(">"):
			name = line[1:]

		else:
			if line in kmer_sequences:
				print("There are multiple %s kmer sequences." %(line))
				sys.exit("ERROR: There are duplicated kmer in your kmer fasta file. Please double check! ")
			else:
				kmer_sequences[line] = [name,0]

print("Finished loading sequences.")

k = len(list(kmer_sequences.keys())[0])
print(k)


## load the atac sequence


atac_sequences = {}

with open(args.atac,"r") as FILE:
	for line in FILE:

		line = line.strip("\n")

			## search for > for the header

		if line.startswith(">"):
			name = line[1:]
			if name not in atac_sequences:
				atac_sequences[name] = ""
			else:
				print("There are multiple %s sequences." %(name))
				sys.exit("ERROR: There are duplicated names in your fasta file. Please double check! ")

		else:
			atac_sequences[name] += line
print("Finished loading atac sequences.")

## count kmers from the atac sequences

for key in atac_sequences:

	sequence = atac_sequences[key]
		
	l = len(sequence)

	start = 0

	end = start + k

	while( start < l - k + 1):

		kmer = sequence[start:end]
		if kmer in kmer_sequences:
			kmer_sequences[kmer][1] += 1
		else:
			sys.exit("ERROR: found a kmer in atac sequence but is not present in the unique kmer file! ")

		start = start + 1
		end = start + k

## load the kmer cluster file (.clstr)

kmer_cluster = {}
cluster_count = {}

with open(args.cluster,"r") as CLUSTER:
	for line in CLUSTER:
		line = line.strip("\n")
		if line.startswith(">"):
			m = re.match(r">Cluster\s(\d+)",line)
			index = m.group(1)
			key_name = "Cluster_" + str(index)
			cluster_count[key_name] = 0
		else:
			m = re.match(r".*>(kmer_\d+)\.\.\..*",line)
			kmer_index = m.group(1)
			kmer_cluster[kmer_index] = "Cluster_" + str(index)

## tally the unique kmer occurance under atac to kmer clusters

for key in kmer_sequences:
	kmer_id = kmer_sequences[key][0]
	count = kmer_sequences[key][1]
	cluster_id = kmer_cluster[kmer_id]
	cluster_count[cluster_id] += count

## load the kmer cluster dosage matrix


X = pd.read_csv(str(args.dm),sep=",")
n,p = X.shape

kmer_cluster_dm_names = np.array(X.columns.values.tolist())

print(len(kmer_cluster_dm_names))

total_cluster_dosage = np.sum(np.asarray(X),axis=0)

with open(args.output+".txt","w") as OUTPUT:
	for i in range(len(kmer_cluster_dm_names)):
		cluster_id = kmer_cluster_dm_names[i]
		print("%s\t%i\t%i" %(cluster_id,cluster_count[cluster_id],total_cluster_dosage[i]),file=OUTPUT)


