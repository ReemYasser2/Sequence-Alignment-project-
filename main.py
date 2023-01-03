from PyQt5 import QtWidgets, uic
import sys
import numpy as np
import subprocess # for MSA using muscle
from Bio import SeqIO
from PyQt5.QtWidgets import *


def stringToList(data):
   return list(data)

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('.\data\mainwindow.ui', self)
        self.show()
        self.global_align_button.clicked.connect(self.global_alignment)
        self.local_align_button.clicked.connect(self.local_alignment)
        self.align_msa_button.clicked.connect(self.multiple_sequence_alignment) 
        # self.align_msa_button.clicked.connect(self.get_input_msa)
        self.actionOpen_Fasta.triggered.connect(self.browse_files)
        self.flag = 0
        self.msa_flag = 0

    def ShowPopUpMessage(self, popUpMessage):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText(popUpMessage)
        msg.setWindowTitle("Error")
        msg.exec_()

    def browse_files(self): 
        try:
            self.msa_flag = 1
            self.flag = 1 # to choose between input in gui and from fasta file flag = 1 is for fasta file
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

    def get_input_pairwise(self):
        try:
            if self.flag == 0: # by default read from the gui
                # read user input sequences from the gui 
                seq1_in = self.pairwise_seq1_in_line.toPlainText()
                seq2_in = self.pairwise_seq2_in_line.toPlainText()
            elif self.flag == 1:
                self.flag = 0
                # read from the fasta file
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
                # self.msa_flag = 0
                fasta_read = open(self.path)
                sequences= [i for i in SeqIO.parse(fasta_read,'fasta')] # loop over the sequences and place each in a variable
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
                seq1 = self.msa_seq1_in_line.toPlainText()
                seq2 = self.msa_seq2_in_line.toPlainText()
                seq3 = self.msa_seq3_in_line.toPlainText()
                seq4 = self.msa_seq4_in_line.toPlainText()
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
            # traceback step directions
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
            matrix_global = np.zeros((matrix_global_num_rows, matrix_global_num_cols))
            for i in range(matrix_global_num_cols):
                matrix_global[0][i]= gap * i
            for j in range(matrix_global_num_rows):
                matrix_global[j][0]= gap * j

            # TRACEBACK INITIALIZATION
            traceback_matrix = np.zeros((matrix_global_num_rows, matrix_global_num_cols))
            traceback_matrix[0,0] = done
            for j in range(1, matrix_global_num_cols):
                traceback_matrix[0,j] = left_dir
            for i in range(1, matrix_global_num_rows):
                traceback_matrix[i,0] = up_dir

            # SCORE MATRIX AND TRACEBACK MATRIX FILL
            for i in range(1, matrix_global_num_rows):
                for j in range(1, matrix_global_num_cols):
                    if seq1[j-1] == seq2[i-1]:
                        score = match
                    else:
                        score = mismatch
                    res = max(matrix_global[i-1][j-1] + score, matrix_global[i-1][j] + gap, matrix_global[i][j-1] + gap)
                    # fill the matrix with the max score at each step
                    matrix_global[i][j] = res
                    # fill the traceback matrix with the direction or arrows for traceback
                    if res == matrix_global[i-1][j-1] + score:
                        traceback_matrix[i][j] = diag_dir
                    if res == matrix_global[i-1][j] + gap:
                        traceback_matrix[i][j] = up_dir
                    if res == matrix_global[i][j-1] + gap:
                        traceback_matrix[i][j] = left_dir

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
                "-align", r".\Group_1.fasta",
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
                        seq+= "\n\n\n"
                        continue
                    seq += line.strip()
                self.msa_output_line.clear()
                self.msa_output_line.setText(seq)

        except:
            self.ShowPopUpMessage("An error has occured wile aligning the sequences")

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()
