from PyQt5 import QtWidgets, uic
import sys
import numpy as np
import subprocess # for MSA using muscle

def stringToList(data):
   return list(data)

class Ui(QtWidgets.QMainWindow):
    def __init__(self):
        super(Ui, self).__init__()
        uic.loadUi('mainwindow.ui', self)
        self.show()
        self.global_align_button.clicked.connect(self.global_alignment)
        self.align_msa_button.clicked.connect(self.multiple_sequence_alignment)
        self.local_align_button.clicked.connect(self.local_alignment)
        # get_input_msa

    def get_input_pairwise(self):

        # read user input sequences from the gui 
        seq1_in = self.pairwise_seq1_in_line.text()
        seq2_in = self.pairwise_seq2_in_line.text()
        # convert the string to a list to deal with more easily in the alignment functions
        seq1 = stringToList(seq1_in)
        seq2 = stringToList(seq2_in)
        # read the scores from the user
        match = int(self.pair_match_input.value())
        mismatch = int(self.pair_mismatch_input.value())
        gap = int(self.pair_gap_input.value())
        return [seq1, seq2, match, mismatch, gap]

    def get_input_msa(self):
        # read user input sequences from the gui 
        seq1_in = self.msa_seq1_in_line.text()
        seq2_in = self.msa_seq2_in_line.text()
        # convert the string to a list to deal with more easily in the alignment functions
        seq1 = stringToList(seq1_in)
        seq2 = stringToList(seq2_in)
        # read the scores from the user
        match = int(self.msa_match_input.value())
        mismatch = int(self.msa_mismatch_input.value())
        gap = int(self.msa_gap_input.value())
        return [seq1, seq2, match, mismatch, gap]

    def global_alignment(self):
        seq1, seq2, match, mismatch, gap = self.get_input_pairwise()
        seq1 = "agc"
        seq2 = "aaac"
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
                matrix_global[i][j] = res
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
        self.pairwise_output_line.setText(str(output))

    def local_alignment(self):
        print(1)
    
    def multiple_sequence_alignment(self):
        seq1, seq2, match, mismatch, gap = self.get_input_msa()
        output = str(match) + " " +  str(mismatch) + " " + str(gap)
        self.msa_output_line.setText(str(output))     

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()




