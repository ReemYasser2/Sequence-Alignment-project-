from PyQt5 import QtWidgets, uic, QtGui
import sys
import numpy as np
import subprocess
from Bio import SeqIO
from PyQt5.QtWidgets import *
from math import log
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import pandas as pd
def MI(sequences,i,j):
    # direct appplication for the rule 
    Pi = Counter(sequence[i] for sequence in sequences)
    Pj = Counter(sequence[j] for sequence in sequences)
    Pij = Counter((sequence[i],sequence[j]) for sequence in sequences)   

    return sum(Pij[(x,y)]*log(Pij[(x,y)]/(Pi[x]*Pj[y])) for x,y in Pij)

def stringToList(data):
   return list(data)
   
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

plt.rcParams["figure.autolayout"] = True
# plt.rcParams['axes.facecolor'] = 'black'
plt.rc('axes', edgecolor='w')
plt.rc('xtick', color='w')
plt.rc('ytick', color='w')

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('.\data\mainwindow.ui', self)
       
        self.setWindowIcon(QtGui.QIcon('./data/icons/icon.png')) # icon added
        self.setWindowTitle("Sequence Alignment Viewer")

        self.show()
        # alignment function calls using buttons
        self.global_align_button.clicked.connect(self.global_alignment)
        self.local_align_button.clicked.connect(self.local_alignment)
        self.align_msa_button.clicked.connect(self.multiple_sequence_alignment) 
        
        # clear button calls
        self.msa_clear_button.clicked.connect(self.clear_msa_fields)
        self.pairwise_clear_button.clicked.connect(self.clear_pairwise_fields)

        # data visualization function calls
        self.visualize_msa_data_button.clicked.connect(lambda: self.MSA_Visualization())
        self.visualize_msa_data_button.clicked.connect(lambda: self.tabWidget.setCurrentIndex(2))
        self.visualize_msa_data_button.clicked.connect(lambda: self.splitter_9.setSizes([0,100]))

        self.visualize_pairwise_data_button.clicked.connect(lambda: self.pairwise_plots(self.alignment_1_copy, self.alignment_2_copy, self.traceback_copy)) 
        self.visualize_pairwise_data_button.clicked.connect(lambda: self.tabWidget.setCurrentIndex(2)) 
        self.visualize_pairwise_data_button.clicked.connect(lambda: self.splitter_9.setSizes([100,0]))

        # menubar function calls
        self.actionAbout_Us_2.triggered.connect(self.about_us)
        self.actionOpen_Fasta.triggered.connect(self.browse_files)

        # MSA assessment function calls
        self.assess_msa_button.clicked.connect(lambda:self.msa_assessment(flag=True))
        
        self.sequence_list_1= None
        self.sequence_list_2= None
        self.sequence_list_3= None
        self.sequence_list_4= None
        self.traceback_copy= None
        self.alignment_1_copy =None   
        self.alignment_2_copy =None

        self.flag = 0
        self.msa_flag = 0

    def ShowPopUpMessage(self, popUpMessage):
        # Function to display a certain error message when an error occurs
        QMessageBox.warning(
            self, ' ERROR',popUpMessage)


    def browse_files(self): 
        try:
            # to choose between input in gui and from fasta file flag = 1 is for fasta file
            self.msa_flag = 1
            self.flag = 1 
            file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'Open File') 
            self.path = file_name[0] #variable to store opened file path
            fasta_read = open(self.path) # read a fasta  file
            sequences= [i for i in SeqIO.parse(fasta_read,'fasta')] # read multiple sequences from the file
            # store each sequence in a variable
            self.sequence_1= sequences[0].seq
            self.sequence_2=sequences[1].seq
            self.sequence_3=sequences[2].seq
            self.sequence_4=sequences[3].seq
        except:
            self.ShowPopUpMessage("An error has occured wile opening the file")

    def clear_msa_fields(self):
        # clear all fields that had anything in them 
        self.msa_seq1_in_line.clear()
        self.msa_seq2_in_line.clear()
        self.msa_seq3_in_line.clear()
        self.msa_seq4_in_line.clear()
        self.msa_output_line.clear()
        self.lcd_mutual_information.display(0)
        self.lcd_percent_identity.display(0)
        self.lcd_sum_of_pairs.display(0)
    
    def clear_pairwise_fields(self):
        # clear all fields that had anything in them 
        self.pairwise_seq1_in_line.clear()
        self.pairwise_seq2_in_line.clear()
        self.pairwise_output_line.clear()


    def get_input_pairwise(self):
        try:
            if self.flag == 0: # by default read from the gui
                # read user input sequences from the gui 
                seq1_in = self.pairwise_seq1_in_line.toPlainText()
                seq2_in = self.pairwise_seq2_in_line.toPlainText()
            elif self.flag == 1:
                self.flag = 0
                # read from the fasta file (take first two sequences for pairwise alignment)
                seq1_in = self.sequence_1
                seq2_in = self.sequence_2
                # display the read sequences for the user to see in the gui
                self.pairwise_seq1_in_line.clear()
                self.pairwise_seq2_in_line.clear()
                self.pairwise_seq1_in_line.setPlainText(str(self.sequence_1))
                self.pairwise_seq2_in_line.setPlainText(str(self.sequence_2))
            
            # convert the string to a list to deal with more easily in the alignment functions
            seq1 = stringToList(seq1_in)
            seq2 = stringToList(seq2_in)
            # read the scores from the user
            match = int(self.pair_match_input.value())
            mismatch = int(self.pair_mismatch_input.value())
            gap = int(self.pair_gap_input.value())
            return [seq1, seq2, match, mismatch, gap]
        except:
            self.ShowPopUpMessage("An error has occured wile obtaining the data")

    def get_input_msa(self):
        try:
            if self.msa_flag == 1:
                # read the input fasta file and display the sequences inside it in the input fields of the gui
                fasta_read = open(self.path)
                sequences= [i for i in SeqIO.parse(fasta_read,'fasta')] 
                # loop over the sequences and place each in a variable
                self.sequence_1= sequences[0].seq
                self.sequence_2=sequences[1].seq
                self.sequence_3=sequences[2].seq
                self.sequence_4=sequences[3].seq
                seq1_in = str(self.sequence_1)
                seq2_in = str(self.sequence_2)
                seq3_in = str(self.sequence_3)
                seq4_in = str(self.sequence_4)
                # display the sequences in the input fields
                self.msa_seq1_in_line.clear()
                self.msa_seq2_in_line.clear()
                self.msa_seq3_in_line.clear()
                self.msa_seq4_in_line.clear()
                self.msa_seq1_in_line.setPlainText(seq1_in)
                self.msa_seq2_in_line.setPlainText(seq2_in)
                self.msa_seq3_in_line.setPlainText(seq3_in)
                self.msa_seq4_in_line.setPlainText(seq4_in)
                
            elif self.msa_flag == 0:
                # read the sequences written by the user
                seq1 = self.msa_seq1_in_line.toPlainText()
                seq2 = self.msa_seq2_in_line.toPlainText()
                seq3 = self.msa_seq3_in_line.toPlainText()
                seq4 = self.msa_seq4_in_line.toPlainText()
                seq1= seq1.upper()
                seq2= seq2.upper()
                seq3= seq3.upper()
                seq4= seq4.upper()
                # check if a sequence is empty, user shouldn't enter less than 3 sequences
                #  an error message is displayed if less than 3 sequences are written
                if seq1 == "":
                    all_sequences = [seq2, seq3, seq4]
                if seq2 == "":
                    all_sequences = [seq1, seq3, seq4]
                if seq3 == "":
                    all_sequences = [seq1, seq2, seq4]
                if seq4 == "":
                    all_sequences = [seq1, seq2, seq3]
                else:
                    all_sequences = [seq1, seq2, seq3, seq4]

                file = open(r'.\user_seq.fasta', 'w+')
                out = '\n'.join(['>Sequence' + str(i+1) + "\n" + j for i,j in enumerate(all_sequences)])
                file.write(out)
                file.close()
        except:
            self.ShowPopUpMessage("An error has occured wile opening the file")

    def global_alignment(self):
        try:
            seq1, seq2, match, mismatch, gap = self.get_input_pairwise()
            # traceback step directions, assigning numerical values to each direction 
            done = 5
            up_dir = 1
            left_dir = 2
            diag_dir = 3
            # length of sequences
            len_seq_1 = len(seq1)
            len_seq_2 = len(seq2)
            # length for loops (zero-based)
            matrix_global_num_cols = len_seq_1 + 1 # seq1
            matrix_global_num_rows = len_seq_2 + 1  # seq2
            #-----------------------------------GLOBAL ALIGNMENT-----------------------------------

            # SCORE MATRIX INITIALIZATION
            # initialize first row and column with gapscore * index
            matrix_global = np.zeros((matrix_global_num_rows, matrix_global_num_cols))
            for i in range(matrix_global_num_cols):
                matrix_global[0][i]= gap * i
            for j in range(matrix_global_num_rows):
                matrix_global[j][0]= gap * j

            # TRACEBACK INITIALIZATION
            # initialize cell at (0,0) as end, columns as left,  rows as up, and the rest as zero
            traceback_matrix = np.zeros((matrix_global_num_rows, matrix_global_num_cols))
            traceback_matrix[0,0] = done
            for j in range(1, matrix_global_num_cols):
                traceback_matrix[0,j] = left_dir
            for i in range(1, matrix_global_num_rows):
                traceback_matrix[i,0] = up_dir

            # SCORE MATRIX AND TRACEBACK MATRIX FILL
            # starting from the second row and column fill according to max of match, mismatch, gap
            for i in range(1, matrix_global_num_rows):
                for j in range(1, matrix_global_num_cols):
                    if seq1[j-1] == seq2[i-1]:
                        score = match
                    else:
                        score = mismatch
                    res = max(matrix_global[i-1][j-1] + score, matrix_global[i-1][j] + gap, matrix_global[i][j-1] + gap)
                    # fill the matrix with the max score at each step
                    matrix_global[i][j] = res
                    # fill the traceback matrix with the direction or arrows 
                    if res == matrix_global[i-1][j-1] + score:
                        traceback_matrix[i][j] = diag_dir
                    if res == matrix_global[i-1][j] + gap:
                        traceback_matrix[i][j] = up_dir
                    if res == matrix_global[i][j-1] + gap:
                        traceback_matrix[i][j] = left_dir
            
            # save the score matrix and traceback matrix for later use (in plots)
            self.traceback_copy = traceback_matrix

            aligned_seq1 = []
            aligned_seq2 = []
            i = len_seq_2 # rows
            j = len_seq_1 # cols 
            # write into the aligned arrays, the nucleotide or added gap 
            # according to the direction in the traceback matrix
            while(i > 0 or j > 0):
                if traceback_matrix[i,j] == diag_dir:
                    aligned_seq1.append(seq1[j-1])
                    aligned_seq2.append(seq2[i-1])
                    i -= 1
                    j -= 1
                elif traceback_matrix[i,j] == left_dir:
                    aligned_seq1.append(seq1[j-1])
                    aligned_seq2.append("-" )
                    j -= 1
                elif traceback_matrix[i,j] == up_dir:
                    aligned_seq1.append("-")
                    aligned_seq2.append(seq2[i-1])
                    i -= 1
                elif traceback_matrix[i,j] == done:
                    break

            # format the output and display in the gui
            aligned_seq1.reverse()
            aligned_seq2.reverse()
            separator = ""
            output_seq1 = separator.join(aligned_seq1)
            output_seq2 = separator.join(aligned_seq2)
            alignment_score = str(matrix_global[-1,-1])
            output_seq1 += "\n"
            output_seq2 += "\n"
            alignment_score += "\n"

            # save the output for later use (in plots)
            self.alignment_1_copy = output_seq1
            self.alignment_2_copy = output_seq2

            output = output_seq1 + output_seq2 + alignment_score
            self.pairwise_output_line.clear()
            self.pairwise_output_line.setText(str(output))
        except:
            self.ShowPopUpMessage("An error has occured wile aligning the sequences")

    def local_alignment(self):
        try:
            seq1_in, seq2_in, match, mismatch, gap = self.get_input_pairwise()
            # input sequences
            seq1 = stringToList(seq1_in)
            seq2 = stringToList(seq2_in)
            # length of sequences
            len_seq_1 = len(seq1)
            len_seq_2 = len(seq2)
            # length for zero-based indexing
            matrix_local_num_cols = len_seq_1 + 1 # seq1
            matrix_local_num_rows = len_seq_2 + 1  # seq2
            # traceback step directions
            done = 5
            up_dir = 1
            left_dir = 2
            diag_dir = 3
            # # TRACEBACK INITIALIZATION
            traceback_matrix = np.zeros((matrix_local_num_rows, matrix_local_num_cols))
            traceback_matrix[0,0] = done
            for j in range(1, matrix_local_num_cols):
                traceback_matrix[0,j] = left_dir
            for i in range(1, matrix_local_num_rows):
                traceback_matrix[i,0] = up_dir
            # matrix initialization
            matrix_local = np.zeros((len_seq_2  + 1, len_seq_1  + 1))
            # Alignment Matrix and Traceback Matrix Fill
            for i in range(1, len_seq_2 + 1):
                for j in range(1, len_seq_1 + 1):
                    if seq1[j-1] == seq2[i-1]:
                        score = match
                    else:
                        score = mismatch
                    # fill the matrix with the max score at each step                   
                    res = max(0, matrix_local[i-1][j-1] + score, matrix_local[i-1][j] + gap, matrix_local[i][j-1] + gap)
                    matrix_local[i][j] = res
                    # fill the traceback matrix with the direction or arrows for traceback
                    if res == matrix_local[i-1][j-1] + score:
                        traceback_matrix[i][j] = diag_dir
                    if res == matrix_local[i-1][j] + gap:
                        traceback_matrix[i][j] = up_dir
                    if res == matrix_local[i][j-1] + gap:
                        traceback_matrix[i][j] = left_dir

            alignment_score = np.max(matrix_local)
            start_loc = np.where(matrix_local == alignment_score) 
            start_loc_col = start_loc[1]
            start_loc_row = start_loc[0]

            aligned_seq1 = []
            aligned_seq2 = []
            i = start_loc_row[0] # rows
            j = start_loc_col[0] # cols
            
            # save the score matrix and traceback matrix for later use (in plots)
            self.traceback_copy = traceback_matrix
            
            # write into the aligned arrays, the nucleotide or added gap 
            # according to the direction in the traceback matrix
            while(i > 0 or j > 0):
                if traceback_matrix[i,j] == diag_dir:
                    aligned_seq1.append(seq1[j-1])
                    aligned_seq2.append(seq2[i-1])
                    i -= 1
                    j -= 1
                elif traceback_matrix[i,j] == left_dir:
                    aligned_seq1.append(seq1[j-1])
                    aligned_seq2.append("-" )
                    j -= 1
                elif traceback_matrix[i,j] == up_dir:
                    aligned_seq1.append("-")
                    aligned_seq2.append(seq2[i-1])
                    i -= 1
                elif traceback_matrix[i,j] == done:
                    break

            # format the output and display in the gui
            aligned_seq1.reverse()
            aligned_seq2.reverse()
            separator = ""
            output_seq1 = separator.join(aligned_seq1)
            output_seq2 = separator.join(aligned_seq2)
            
            out_score = str(alignment_score)
            output_seq1 += "\n"
            output_seq2 += "\n"
            out_score += "\n"
            
            # save the output for later use (in plots)
            self.alignment_1_copy = output_seq1
            self.alignment_2_copy = output_seq2
            
            output = output_seq1 + output_seq2 + out_score
            self.pairwise_output_line.clear()
            self.pairwise_output_line.setText(str(output))
        except:
            self.ShowPopUpMessage("An error has occured wile aligning the sequences")

    def multiple_sequence_alignment(self):
        try:
            self.get_input_msa()
            # use muscle to align multiple sequences
            if self.msa_flag == 1:
                self.msa_flag = 0
                output = subprocess.check_output(
                ["F:\Installations\muscle5.1.win64.exe",
                "-align", self.path,
                "-output", r".\aligned.fasta"],
                text=True)
                # open the file to display the results of the multiple sequence alignment
                file =open("aligned.fasta")
                seq = ""
                for line in file:
                    if line.startswith(">"): 
                        seq+= "\n\n"
                        continue
                    seq += line.strip()
                self.msa_output_line.clear()
                self.msa_output_line.setText(seq)
            elif self.msa_flag == 0:
                output = subprocess.check_output(
                ["F:\Installations\muscle5.1.win64.exe",
                "-align", r".\user_seq.fasta",
                "-output", r".\user_aligned.fasta"],
                text=True)
                # open the file to display the results of the multiple sequence alignment
                file = open("user_aligned.fasta")
                seq = ""
                for line in file:
                    if line.startswith(">"): 
                        seq+= "\n\n"
                        continue
                    seq += line.strip()
                self.msa_output_line.clear()
                self.msa_output_line.setText(seq)

        except:
            self.ShowPopUpMessage("An error has occured wile aligning the sequences")

    def msa_assessment(self,flag=False):
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
        seq1_str = seq1_str.upper()
        seq2_str = seq2_str.upper()
        seq3_str = seq3_str.upper()
        seq4_str = seq4_str.upper()
        seq_len = len(seq4_str) - 1

        self.sequence_list_1 = stringToList(seq1_str)
        self.sequence_list_2 = stringToList(seq2_str)
        self.sequence_list_3 = stringToList(seq3_str)
        self.sequence_list_4 = stringToList(seq4_str)
        sequences_matrix = [self.sequence_list_1, self.sequence_list_2, self.sequence_list_3, self.sequence_list_4]
      
        
        # calculate the percent identity
        total_pairs=0
        for i in range(len(sequences)): 
            for j in range(i+1,len(sequences)): 
                seq1=sequences[i].seq 
                seq2=sequences[j].seq 

                    #For percent identity analysis
                identical_pairs_count,all_pairs_count=compare(seq1,seq2)
                total_pairs+=all_pairs_count

            percent_Identity=identical_pairs_count/total_pairs*100 


        # Calculate the SOP for all pairs of sequences
        sop_total = 0
        #setting match score with 3, and mismatch with -1
        pair_scores = {'AA':3, 'CC':3, 'GG':3, 'TT':3, 'AT':-1, 'AG':-1, 'AC':-1, 'GA':-1, 'GC':-1, 'GT':-1, 'TC':-1, 'TA':-1, 'TG':-1, 'CA':-1, 'CT':-1, 'CG':-1}
        #calculate sop for each pair of sequences
        for i in range(len(sequences_matrix)):
            for j in range(i+1, len(sequences_matrix)):
                sop = calc_sop(sequences_matrix[i], sequences_matrix[j], pair_scores)
                sop_total += sop

        formatted_SOP = "{:.2f}".format(sop_total)
        print(f'Total SOP: {sop_total}')

        formatted_PI = "{:.2f}".format(percent_Identity)
        print(formatted_PI)
        
        out = MI(sequences_matrix,seq_len,seq_len)



        formatted_mi = "{:.2f}".format(out)
        
        if flag == True:
            self.lcd_mutual_information.display(formatted_mi)
            self.lcd_percent_identity.display(formatted_PI)
            self.lcd_sum_of_pairs.display(formatted_SOP)

    def pairwise_plots(self, aligned_seq1,aligned_seq2, traceback_matrix):
        self.heatmap_parwise_plot(traceback_matrix)
        self.correlation_parwise_plot(aligned_seq1,aligned_seq2)

    def draw_MSA_canvas(self, image, layout):
        self.figure = plt.figure(figsize=(15,5))
        self.figure.patch.set_facecolor('#140826')
        self.Canvas = FigureCanvas(self.figure)
        layout.addWidget(self.Canvas,0, 0, 1, 1)
        plt.imshow(image, aspect= 5)
        self.Canvas.draw() 

    def heatmap_parwise_plot(self, traceback_matrix):
        plt.pcolormesh(traceback_matrix,cmap="ocean")
        plt.ylim(len(traceback_matrix),0)
        self.fig = plt.figure(figsize=(15,5))
        self.fig.patch.set_facecolor('#140826')
        self.Canvas = FigureCanvas(self.fig)
        self.pairwise_heatmap.addWidget(self.Canvas,0, 0, 1, 1)
        plt.imshow(traceback_matrix)
        self.Canvas.draw() 

    def correlation_parwise_plot(self, aligned_seq1, aligned_seq2):
        aligned_seq1_copy = []
        aligned_seq2_copy = []
        
        for i in range(len(aligned_seq1)):
            if aligned_seq1[i] == 'a' or aligned_seq1[i] == 'A':
                aligned_seq1_copy.append(1)
            elif aligned_seq1[i] == 'c' or aligned_seq1[i] == 'C':
                aligned_seq1_copy.append(4)
            elif aligned_seq1[i] == 'g' or aligned_seq1[i] == 'G':
                aligned_seq1_copy.append(2)
            elif aligned_seq1[i] == 't' or aligned_seq1[i] == 'T':
                aligned_seq1_copy.append(3)
            elif aligned_seq1[i] == '-':
                aligned_seq1_copy.append(0)
        
        for i in range(len(aligned_seq2)):
            if aligned_seq2[i] == 'a' or aligned_seq2[i] == 'A':
                aligned_seq2_copy.append(2)
            elif aligned_seq2[i] == 'c' or aligned_seq2[i] == 'C':
                aligned_seq2_copy.append(1)
            elif aligned_seq2[i] == 'g' or aligned_seq2[i] == 'G':
                aligned_seq2_copy.append(3)
            elif aligned_seq2[i] == 't' or aligned_seq2[i] == 'T':
                aligned_seq2_copy.append(4)
            elif aligned_seq2[i] == '-':
                aligned_seq2_copy.append(5)

        y = pd.Series(aligned_seq1_copy)
        x = pd.Series(aligned_seq2_copy)
        correlation = y.corr(x)
        correlation
        
        self.fig2 = plt.figure(figsize=(15,5))
        self.fig2.patch.set_facecolor('#140826')
        self.Canvas2 = FigureCanvas(self.fig2)
        self.pairwaise_correlation.addWidget(self.Canvas2,0, 0, 1, 1)
        plt.scatter(x, y)
        plt.plot(np.unique(x),
                np.poly1d(np.polyfit(x, y, 1))
                (np.unique(x)), color='red')
        
        self.Canvas2.draw()

    def MSA_Visualization(self):
        self.msa_assessment()
        mat = []
        for i in range(4):
            mat.append(self.sequence_list_1)
        for i in range(4):
            mat.append(self.sequence_list_2)
        for i in range(4):
            mat.append(self.sequence_list_3)
        for i in range(4):
            mat.append(self.sequence_list_4)

        seq_arr = np.asarray(mat)

        seq_arr[seq_arr=='A'] = '1' 
        seq_arr[seq_arr=='G'] = '2' 
        seq_arr[seq_arr=='T'] = '3'
        seq_arr[seq_arr=='C'] = '4' 
        seq_arr[seq_arr=='-'] = '0'
        # with gaps
        # a blue
        # - purple
        # g dark green
        # t light green
        # c yellow

        # no gaps
        # a purple
        # g blue
        # t green
        # c yellow

        int_seq = np.uint8(seq_arr)
        final = int_seq * 50
        self.draw_MSA_canvas(final, self.MSAVisualization)

    def about_us(self):
        QMessageBox.about(
            self,
            " About ",
            "This is an Sequence Alignment Viewer \nCreated by senior student from the faculty of Engineering, Cairo Uniersity, Systems and Biomedical Engineering department \n \nTeam members: \n- Abdullah Saeed \n- Ro'aa Ehan \n- Farah Gamal \n- Huda Saeed \n- Reem Yasser \nhttps://github.com/Abdullahsaeed10/Sequence-Alignment-project- ",
        )


    






app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()
