import re
import argparse
import numpy as np
import sys
import pandas as pd
import os
from sklearn import metrics
from sklearn.cluster import DBSCAN
import scipy.cluster.hierarchy as sch
from sklearn.neighbors import kneighbors_graph
from scipy.spatial.distance import pdist, squareform,cdist
import time
import subprocess



def cdhit_cluster(unique_kmer_fa,output,similarity,wordsize):
	command = "cd-hit-est" + " -i " + str(unique_kmer_fa) + " -o " + str(output) + " -c " + str(similarity) + " -n " + str(wordsize) +" -r 1" +" -g 1 " + "-gap -20 " + "-gap-ext -10 " + "-l 5 " + "-M 0"

	#print(command)
	subprocess.check_call(command,shell=True)


def kmer_one_hot_encoding(kmer):
	kmer_numeric = np.zeros(len(kmer)*4)

	for i in range(len(kmer)):
		if kmer[i] == "A" or kmer[i] == "a":
			kmer_numeric[i*4] = 1
		elif kmer[i] == "T" or kmer[i] == "t":
			kmer_numeric[i*4 + 1] = 1
		elif kmer[i] == "G" or kmer[i] == "g":
			kmer_numeric[i*4 + 2] = 1
		elif kmer[i] == "C" or kmer[i] == "c":
			kmer_numeric[i*4 + 3] = 1
		else:
			sys.exit("ERROR: found non-conventional nucleotide")
	return(kmer_numeric)

def kmer_one_hot_decoding(kmer_numeric):
	kmer = []

	k = int(len(kmer_numeric) / 4)

	for i in range(k):
		if kmer_numeric[i*4] == 1:
			kmer.append("A")
		elif kmer_numeric[i*4 + 1] == 1:
			kmer.append("T")
		elif kmer_numeric[i*4 + 2] == 1:
			kmer.append("G")
		elif kmer_numeric[i*4 + 3] == 1:
			kmer.append("C")
	kmer_decode = "".join(kmer)

	return(kmer_decode)

def reorder_kmers(kmer_matrix):
	r,c = kmer_matrix.shape
	kmer_matrix = kmer_matrix.astype(np.float32)
	batch_size=100
	degree = []
	for i in range(0,r,batch_size):
		X = kmer_matrix[i:i+batch_size,:]
		distance_matrix = cdist(X, kmer_matrix, metric='euclidean')
		affinity_matrix = np.exp(- distance_matrix * distance_matrix)
		degree.extend(np.array(np.sum(affinity_matrix,axis=1) -1))

	order = np.argsort(degree)[::-1]
	return(order)

# def write_kmer_output(ordered_ALL_K_mers):
# 	with open("tmp_ordered_unique_kmer.fa","w") as OUTPUT:
# 		for i in range(len(ordered_ALL_K_mers)):
# 			print(">kmer_%i\n%s" %(i,ordered_ALL_K_mers[i]),file=OUTPUT)

# def kmer_clustering(ALL_K_mers):

# 	k = len(ALL_K_mers[0])
# 	nrow = len(ALL_K_mers)
# 	ncol = k*4
# 	kmer_matrix = np.empty((nrow,ncol))

# 	one_mismatch = []

# 	for I in range(len(ALL_K_mers)):
# 		kmer = ALL_K_mers[I]
# 		kmer_numeric = kmer_one_hot_encoding(kmer)
# 		kmer_matrix[I,:] = kmer_numeric

# 	order = reorder_kmers(kmer_matrix)

# 	ordered_ALL_K_mers = [ALL_K_mers[i] for i in order]

# 	write_kmer_output(ordered_ALL_K_mers)

# 	identity_threshold = np.round(1 - float(1/k),3)
# 	print(identity_threshold)
# 	cd_hit_command = "cd-hit-est -i tmp_ordered_unique_kmer.fa -o tmp_ordered_unique_kmer -l 5 -c " + str(identity_threshold) + " -M 0 -n 5 -g 1 -gap -20 -gap-ext -10"
# 	print(cd_hit_command)
# 	format_output_command = "cat tmp_ordered_unique_kmer.clstr | sed -E 's/.*>(kmer_.*)\.\.\..*/\\1/'"
# 	print(format_output_command)
	# ## Use manhattan distance to find reference kmers

	# manhattan_distances = np.array(squareform(pdist(kmer_matrix, metric='cityblock')),dtype=int)
	# one_mismatch = np.sum(manhattan_distances == 2,axis=1)

	# kmer_rank_index = np.array(np.argsort(one_mismatch)[::-1],dtype=int)

	# K_mer_clusters = {}

	# ALL_K_mers_tmp = [ALL_K_mers[i] for i in kmer_rank_index]

	# index_to_be_remove_1 = []

	# ## with 1 mismatch
	# for I in range(len(kmer_rank_index)):

	# 	index = kmer_rank_index[I]

	# 	kmer_ = ALL_K_mers[index]

	# 	if kmer_ in  ALL_K_mers_tmp:

	# 		#kmers with 1 mismatch to I only
	# 		index1 = np.where(manhattan_distances[index,:] ==2 )[0]
	# 		# initiate the cluster with a reference
	# 		K_mer_clusters[kmer_] = [kmer_]
	# 		for J in index1:
	# 			if sum(manhattan_distances[J,:] ==2) == 1:
	# 				K_mer_clusters[kmer_].append(ALL_K_mers[J])
	# 				ALL_K_mers_tmp.remove(ALL_K_mers[J])
	# 				index_to_be_remove_1.append(J)

	# kmer_rank_index = [item for item in kmer_rank_index if item not in index_to_be_remove_1]
	# manhattan_distances_removed = np.delete(manhattan_distances,index_to_be_remove_1,axis=1)
	# print(manhattan_distances_removed.shape)

	# ## with 2 mismatch

	# for I in range(len(kmer_rank_index)):

	# 	index = kmer_rank_index[I]

	# 	kmer_ = ALL_K_mers[index]
	# 	if kmer_ in  ALL_K_mers_tmp:

	# 		#kmers with 1 mismatch to I only
	# 		index1 = np.where(manhattan_distances[index,:] == 4 )[0]
	# 		index1 = [item for item in index1 if item not in index_to_be_remove_1]
	# 		# initiate the cluster with a reference
	# 		#key_name = "Kmer_cluster_"+str(I)
	# 		#K_mer_clusters[key_name] = [kmer_]
	# 		for J in index1:
	# 			#print(J,sum(manhattan_distances_removed[J,:] == 4))
	# 			if sum(manhattan_distances_removed[J,:] == 4) == 1:
	# 				print(J,kmer_,ALL_K_mers[J])
	# 				K_mer_clusters[kmer_].append(ALL_K_mers[J])
	# 				ALL_K_mers_tmp.remove(ALL_K_mers[J])
	
	# print(len(K_mer_clusters))

	# # Get the first 5 key-value pairs
	# first_five = list(K_mer_clusters.items())[:5]

	# # Print them
	# for key, value in first_five:
	# 	print(key, value)




	# for I in range(len(ALL_K_mers)):
		# one_mismatch.append(sum(manhattan_distances[I,:] == 2))

	#print(one_mismatch)


	
	# print(ALL_K_mers[0],ALL_K_mers[119])
	# linkage_matrix = sch.linkage(kmer_matrix, method='single')
	# print(linkage_matrix.shape)
	# clusters = sch.fcluster(linkage_matrix, t=1.5, criterion='distance')
	# print(len(clusters))

	# dists = squareform(pdist((kmer_matrix)))
	# affinity_matrix = np.exp( - dists * dists )
	# np.fill_diagonal(affinity_matrix, 0)

	#print(ALL_K_mers[:10],kmer_matrix[:10,:])

	##manhattan_distances = squareform(pdist(kmer_matrix, metric='cityblock'))
	#print(sum(manhattan_distances[0,:]==2))

	# connectivity = kneighbors_graph(kmer_matrix, n_neighbors=k*3, mode='connectivity', include_self=False)
	# connectivity_matrix = connectivity.toarray()
	# weighted_connectivity = np.multiply(affinity_matrix,connectivity_matrix)
	# #G = knn_graph.toarray()
	# d = np.sum(weighted_connectivity,axis=1)
	# print(np.unique(d))
	# print(np.sort(d)[::-1])
	# print(sum(d==max(d)))

