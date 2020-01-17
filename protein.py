import os
from Bio import SeqIO
import numpy as np
from collections import Counter
import argparse
from sys import platform
import subprocess


class Protein():
    def __init__(self, fasta, tau):
        dirpath = os.getcwd()
        self.sequences = list(SeqIO.parse(dirpath +'/'+ fasta, "fasta"))
        self.tau = tau

        MSA = []
        for seq in self.sequences:
            MSA.append(list(str(seq.seq)))
        self.MSA = np.array(MSA)
        (self.N, self.M) = np.shape(MSA)


    def MI(self, sequences, i, j):
        N = np.shape(sequences)[0]
        Pi = Counter(s[i] for s in sequences)
        Pj = Counter(s[j] for s in sequences)
        Pij = Counter((s[i], s[j]) for s in sequences)

        return sum(Pij[(x, y)] / N * np.log(Pij[(x, y)] * N / (Pi[x] * Pj[y])) for x, y in Pij)

    def calculate_m_prime(self, MI_matrix):
        c_mean = np.mean(MI_matrix, 0)  # columns mean
        r_mean = np.mean(MI_matrix, 1)  # rows mean
        MI_mean = [[np.mean([c, r]) for c in c_mean] for r in r_mean]  # matrix of mean
        MI_prime = np.subtract(MI_matrix, MI_mean)
        return MI_prime

    def tau(self, tau, MI_prime):
        MI_thresh = np.array(MI_prime)
        for i in range(len(MI_prime)):
            for j in range(len(MI_prime[i])):
                MI_thresh[i][j] = int(MI_prime[i][j] > tau)
        return MI_thresh

    def remove_zeros(self, MI_thresh):
        index = []
        for i in range(len(MI_thresh)):
            if np.sum(MI_thresh[i]) == 0:
                index.append(i)

        MI_thresh = np.delete(MI_thresh, index, axis=1)
        MI_thresh = np.delete(MI_thresh, index, axis=0)
        return MI_thresh

    def export_cmatrix(self, name, MI_thresh):

        dirpath = os.getcwd()
        f = open(dirpath+ '/' + name + '.txt', 'w')
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
    parser.add_argument("-t", "--threshold", help="Threshold", required=False, default=0.01)
    args = parser.parse_args()

    protein = Protein(args.fasta, args.threshold)
    print('Computing the MI Matrix, patiente this process would take time')
    MI_matrix = [[protein.MI(protein.sequences, i, j) for i in range(protein.M)] for j in range(protein.M)]
    print('Reducing the MI Matrix')
    MI_matrix = protein.calculate_m_prime(MI_matrix)
    MI_matrix = protein.tau(args.threshold, MI_matrix)
    MI_matrix = protein.remove_zeros(MI_matrix)
    print('Exporting the MI Matrix')
    protein.export_cmatrix(args.output, MI_matrix)
