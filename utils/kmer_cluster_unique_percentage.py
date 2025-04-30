
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-c',action= 'store',dest='cluster',help = 'the kmer clustering output')
parser.add_argument('-m',action= 'store',dest='mapping',help='the kmer pip file')

args = parser.parse_args()


cluster = {}
total_kmer = 0

with open(args.cluster,"r") as CLUSTER:
	for line in CLUSTER:
		line = line.strip("\n")
		if line.startswith(">"):
			m = re.match(r">Cluster\s(\d+)",line)
			index = m.group(1)
			cluster[index] = []
		else:
			m = re.match(r".*>kmer_(\d+)\.\.\..*",line)
			kmer_index = m.group(1)
			cluster[index].append(kmer_index)
			total_kmer += 1

cluster_mapping = []
with open(args.mapping,"r") as MAPPING:
	header = MAPPING.readline()
	for line in MAPPING:
		items = line.split("\t")
		m = re.match(r"Cluster_(\d+)",items[0])
		cluster_mapping.append(m.group(1))

used_kmer = 0
for i in range(len(cluster_mapping)):
	count = len(cluster[cluster_mapping[i]])
	used_kmer += count

kmer_usage_percent = round(float(used_kmer) / total_kmer,5)
cluster_percent = round(float(len(cluster_mapping))/len(cluster),5)
print("total_kmer\tused_kmer\tkmer_usage_percentage\ttotal_cluster\tused_cluster\tcluster_usage_percentage")
print("%d\t%d\t%f\t%d\t%d\t%f" %(total_kmer,used_kmer,kmer_usage_percent,len(cluster),len(cluster_mapping),cluster_percent))
