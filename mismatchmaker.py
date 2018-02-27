## create all possible mismatches of a sequence
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("seq",help="Sequence to be mismatched", type=str)
args = parser.parse_args()
seq = args.seq
# dictionary #
NTPs=['A','C','T','G']
# mismatching #
k=0
print(seq)
for base in seq:
    k+=1
    for nt in NTPs:
        if base != nt:
            print(seq[0:k-1]+nt+seq[k:len(seq)])
