from math import log
from collections import Counter
from Bio import SeqIO
fasta_read = open("aligned.fasta") # read a fasta  file
sequences= [i for i in SeqIO.parse(fasta_read,'fasta')] # read multiple sequences from the file

# store each sequence in a variable
sequence_1= sequences[0].seq
sequence_2=sequences[1].seq
sequence_3=sequences[2].seq
sequence_4=sequences[3].seq

seq1_str = str(sequence_1)
seq2_str = str(sequence_2)
seq3_str = str(sequence_3)
seq4_str = str(sequence_4)
seq_len = len(seq4_str) - 1

def stringToList(data):
   return list(data)

seq1 = stringToList(seq1_str)
seq2 = stringToList(seq2_str)
seq3 = stringToList(seq3_str)
seq4 = stringToList(seq4_str)
sequences_matrix = [seq1, seq2, seq3, seq4]



def MI(sequences,i,j):
    Pi = Counter(sequence[i] for sequence in sequences)
    Pj = Counter(sequence[j] for sequence in sequences)
    Pij = Counter((sequence[i],sequence[j]) for sequence in sequences)   

    return sum(Pij[(x,y)]*log(Pij[(x,y)]/(Pi[x]*Pj[y])) for x,y in Pij)

out = MI(sequences_matrix,seq_len,seq_len)
formatted_out = "{:.2f}".format(out)
print(formatted_out)


