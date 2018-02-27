# create all possible mismatches of a sequence
# input is the list of all experiment barcodes
# output is STDOUT
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("seq",help="Sequence to be mismatched", type=str)
args = parser.parse_args()
seq = args.seq
# dictionary #
NTPs=['A','C','T','G']
# mismatching #
k=0
# print the original sequence
print(seq)
for base in seq:
    k+=1
    for nt in NTPs:
        # print all sequences that do not match with the original one in only one position
        if base != nt:
            print(seq[0:k-1]+nt+seq[k:len(seq)])
