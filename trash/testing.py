from Bio import SeqIO
fasta_read = open("Group_1.fasta")
sequences= [i for i in SeqIO.parse(fasta_read,'fasta')]
     
# print (len(sequences))
sequence_1= sequences[0].seq
print (sequence_1)
print("------------------------------------------------------------------------")
sequence_2=sequences[1].seq
print (sequence_2)
print("-------------------------------------------------------------------------")
sequence_3=sequences[2].seq
print (sequence_3)
print("-------------------------------------------------------------------------")
sequence_4=sequences[3].seq
print (sequence_4)
print("--------------------------------------------------------------------------")


                