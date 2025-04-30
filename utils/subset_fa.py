

### This python script will extract sequences in file A that also present in file B based on sequence header


import argparse
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-a',action= 'store',dest='filea',help = 'the file A .fa')
parser.add_argument('-b',action= 'store',dest='fileb',help = 'the file B .fa')
parser.add_argument('-o',action= 'store',dest='output',help = 'prefix of the output file')
args = parser.parse_args()


## read sequence B

sequences_B = []

with open(args.fileb,"r") as FILEB:
	for line in FILEB:

		line = line.strip("\n")

			## search for > for the header

		if line.startswith(">"):
			name = line[1:]
			if name in sequences_B:
				print("There are multiple %s sequences." %(name))
				sys.exit("ERROR: There are duplicated names in your fasta B file. Please double check! ")
			else:
				sequences_B.append(name)

print("Finished loading sequence file %s" %(args.fileb))

## read sequence A

sequences_A = {}

with open(args.filea,"r") as FILEA:
	for line in FILEA:

		line = line.strip("\n")

		## search for > for the header

		if line.startswith(">"):
			name = line[1:]
			if name not in sequences_A:
				sequences_A[name] = ""
			else:
				print("There are multiple %s sequences." %(name))
				sys.exit("ERROR: There are duplicated names in your fasta file. Please double check! ")

		else:
			sequences_A[name] += line
print("Finished loading sequence file %s" %(args.filea))


with open(args.output+".fa","w") as OUTPUT:
	for key in sequences_A:
		m = re.match(r"(.*)::.*",key)
		name = m.group(1)
		if name in sequences_B:
			print(">%s\n%s" %(key,sequences_A[key]),file=OUTPUT)

