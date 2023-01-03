from Bio import SeqIO
fasta_read = open("aligned.fasta")
# sequences= [i for i in SeqIO.parse(fasta_read,'fasta')] # loop over the sequences and place each in a variable

sequences = list(SeqIO.parse(fasta_read, "fasta"))
total_pairs=0
sop=0
count=0
countall=0 

class my_dictionary(dict): 
    # __init__ function 
    def __init__(self): 
        self = dict()   
    # Function to add key:value 
    def add(self, key, value): 
        self[key] = value 
residue_frequency = my_dictionary() 

def compare(a,b):
        identical_pairs_count = 0
        all_pairs_count=0
        mismatch_pairs_count =0
        gaps_count =0
        for x, y in zip(a, b):
            if x!='-' and y!='-':
                all_pairs_count+=1
                if x == y:
                    identical_pairs_count += 1
                else:
                    mismatch_pairs_count +=1
            else:
                gaps_count +=1
           
        return identical_pairs_count,all_pairs_count


for i in range(len(sequences)): 
 for j in range(i+1,len(sequences)): 
       seq1=sequences[i].seq 
       seq2=sequences[j].seq 

        #For percent identity analysis
       identical_pairs_count,all_pairs_count=compare(seq1,seq2)
       total_pairs+=all_pairs_count

       sop=sop+count
 percent_Identity=identical_pairs_count/total_pairs*100 

formatted_PI = "{:.3f}".format(percent_Identity)
print(formatted_PI)

