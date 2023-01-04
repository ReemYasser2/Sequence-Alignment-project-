from Bio import SeqIO
fasta_read = open("aligned.fasta")
records= [i for i in SeqIO.parse(fasta_read,'fasta')] # loop over the sequences and place each in a variable
sequence_1= records[0].seq
sequence_2=records[1].seq
sequence_3=records[2].seq
sequence_4=records[3].seq
seq1_string = str(sequence_1)
seq2_string = str(sequence_2)
seq3_string = str(sequence_3)
seq4_string = str(sequence_4)
sequences = [seq1_string, seq2_string, seq3_string, seq4_string]
# sequences = [list(seq) for seq in sequences]

# Afunction to calculate the SOP for a pair of sequences
def calc_sop(seq1, seq2, pair_scores):
  sop = 0
  for i in range(len(seq1)):
    #if not a gap, get the score from the dictionary
    if seq1[i] != '-' and seq2[i] != '-':
        sop += pair_scores[(seq1[i]+seq2[i])]
    else:
        #sdet gap penalty with -2
        sop -=2
  return sop

# Calculate the SOP for all pairs of sequences
sop_total = 0
#setting match score with 3, and mismatch with -1
pair_scores = {'AA':3, 'CC':3, 'GG':3, 'TT':3, 'AT':-1, 'AG':-1, 'AC':-1, 'GA':-1, 'GC':-1, 'GT':-1, 'TC':-1, 'TA':-1, 'TG':-1, 'CA':-1, 'CT':-1, 'CG':-1}
#calculate sop for each pair of sequences
for i in range(len(sequences)):
  for j in range(i+1, len(sequences)):
    sop = calc_sop(sequences[i], sequences[j], pair_scores)
    sop_total += sop

print(f'Total SOP: {sop_total}')