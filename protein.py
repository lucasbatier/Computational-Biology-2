import os
from Bio import SeqIO
import numpy as np
from collections import Counter
import argparse

class Protein():
    def __init__(self, fasta):
        dirpath = os.getcwd()
        sequences = list(SeqIO.parse(dirpath +'/'+ fasta, "fasta"))
        MSA = []
        for seq in sequences:
            MSA.append(list(str(seq.seq)))
        (N, M) = np.shape(MSA)
        MSA = np.array(MSA)

    def MI(self, sequences, i, j):
        N = np.shape(sequences)[0]
        Pi = Counter(s[i] for s in sequences)
        Pj = Counter(s[j] for s in sequences)
        Pij = Counter((s[i], s[j]) for s in sequences)

        return sum(Pij[(x, y)] / N * np.log(Pij[(x, y)] * N / (Pi[x] * Pj[y])) for x, y in Pij)

    def m_prime(self, MI_matrix):
        c_mean = np.mean(MI_matrix, 0)  # columns mean
        r_mean = np.mean(MI_matrix, 1)  # rows mean
        MI_mean = [[np.mean([c, r]) for c in c_mean] for r in r_mean]  # matrix of mean
        MI_prime = np.subtract(MI_matrix, MI_mean)
        return MI_prime

    def export_cmatrix(self, name, MI_thresh):

        dirpath = os.getcwd()
        f = open(dirpath+ '/' + name+  '.txt', 'w')
        f.write(str(len(MI_thresh)) + '\n')
        for i in range(len(MI_thresh)):
            for j in range(len(MI_thresh)):
                f.write(str(int(MI_thresh[i][j])))
                f.write('\n')
        f.close()

if __name__ == '__main__':
    # define the program description
    text = "This tool takes a multiple sequence alignment in standard FASTA file format and outputs a protein structure"

    # initiate the parser with a description
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument("-i", "--fasta", help="Name of the FASTA file", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="Name of the file to store the contact matrix",
        required=False,
        default="cmatrix",
    )
    parser.add_argument("-t", "--threshold", help="Threshold", required=False, default=0.)
    args = parser.parse_args()