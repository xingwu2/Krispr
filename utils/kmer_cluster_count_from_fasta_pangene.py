### This python script will count the occurance of a kmer cluster in a fasta file


import argparse
import re
import pandas as pd
import numpy as np
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-i',action= 'store',dest='fasta',help = 'fasta sequence')
parser.add_argument('-k',action= 'store',dest='uniquekmer',help = 'unique kmer fasta file')
parser.add_argument('-c',action= 'store',dest='cluster',help = 'the kmer clustering output (.clstr)')
parser.add_argument('-p',action= 'store',dest='assignment',help = 'cluster pangene assignment')
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

## load the fasta sequence


fasta_sequences = {}

with open(args.fasta,"r") as FILE:
	for line in FILE:

		line = line.strip("\n")

			## search for > for the header

		if line.startswith(">"):
			name = line[1:]
			if name not in fasta_sequences:
				fasta_sequences[name] = ""
			else:
				sys.exit("There are multiple %s sequences." %(name))

		else:
			fasta_sequences[name] += line
print("Finished loading fasta sequences.")

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


## load kmer cluster pangene assignment 

assignment = {}

with open(args.assignment,"r") as ASSIGNMENT:
	for line in ASSIGNMENT:
		line = line.strip("\n")
		items = line.split("\t")
		if items[0] not in assignment:
			assignment[items[0]] = [items[1]]
		else:
			assignment[items[0]].append(items[1])


## count kmers from the fasta sequences

for key in fasta_sequences:

	sequence = fasta_sequences[key]
		
	l = len(sequence)

	start = 0

	end = start + k

	while( start < l - k + 1):

		kmer = sequence[start:end]

		if kmer in kmer_sequences:
			cluster_ = kmer_cluster[kmer_sequences[kmer][0]]
			if cluster_ in assignment and any( p in key for p in assignment[cluster_]):
				kmer_sequences[kmer][1] += 1			
		else:
			print(key,kmer,sequence,"WARNING: found a kmer in fasta sequence but is not present in the unique kmer file! ")
		start = start + 1
		end = start + k



## tally the unique kmer occurance under atac to kmer clusters

for key in kmer_sequences:
	kmer_id = kmer_sequences[key][0]
	count = kmer_sequences[key][1]
	cluster_id = kmer_cluster[kmer_id]
	cluster_count[cluster_id] += count

with open(args.output+".txt","w") as OUTPUT:
	for key in cluster_count:
		print("%s\t%i" %(key,cluster_count[key]),file=OUTPUT)





