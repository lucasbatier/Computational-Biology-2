import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter
import pickle


def save_obj(obj, name):
    '''
    @param obj: Dictionary of k-mers
    @param name: Name of the output file
    @return: Create a pickle file with the content the obj
    '''
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    '''
    @param name: Name of the pickle file to read
    @return: Load in a dictionary the content from the pickle file
    '''
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def k_mer_counting(reads, k=15):
  '''

  @param reads: List of reads
  @param k: Size of the k-mer by default is 15
  @return: A dictionary which contains the subsequences and the number of times
  that this sequence is repeated k-mer
  '''
  kmers = {}
  for N, read in enumerate(reads):
    read_str = str(read.seq)
    for start in range(len(read_str) - k):
      kmer = read_str[start:start+k]
      kmers[kmer] = kmers.get(kmer, 0) + 1
  return kmers

def max_coverage(C=49.82915, hist_kmers):
    '''

    @param C: Standard coverage
    @param hist_kmers: Histogram of k-mers
    @return: tuple which contains the index of the max coverage
    and the corresponding value.
    '''
    max = 0
    for i in range(round(C - C / 2), round(C + C / 2)):
        if hist_kmers[i] > max:
            max = hist_kmers[i]
            C_max = i
    return C_max

def mean_coverage(C, dict_kmers, start_value=5):
    '''

    @param C: Standard coverage
    @param dict_kmers: Dictionary of kmers
    @param start_value: Boundary for counting the mean values
    @return: the mean coverage
    '''
    C_mean = 0
    n = 0
    for v in list(dict_kmers.values()):
        if v > start_value:
            C_mean += v
            n += 1
    C_mean /= n


def pmf_poisson(l, k):
    '''

    @param l: lambda, parameters of Poisson law
    @param k: Value
    @return: Poisson probability mass function (PMF)
    '''
    return l ** k * np.exp(-l) / np.math.factorial(k)


def cdf_poisson(l, k):
    '''

    @param l: lambda, parameters of Poisson law
    @param k: Value
    @return: Poisson cumulative density function (CDF)
    '''
    return np.exp(-l) * np.sum([l ** i / np.math.factorial(i) for i in range(0, k + 1)])


def snp_with_cdf(kmers_0, kmers_1, cmean, error_rate=0.01):
    '''

    @param kmers_0: Dictionary of kmers to find SNPs
    @param kmers_1: Dictionary of kmers to find 
    @param cmean: Mean coverage for the corresponding kmer
    @param error_rate: Minimum error rate to add the SNP
    @return: array of SNPs 
    '''
    kmers_0_snp = []
    for kmer in kmers_0.keys():
    if kmer not in kmers_1.keys() and cdf_poisson(cmean, kmers_0[kmer]) > error_rate:
        kmers_0_snp.append(kmer)
    return kmers_0_snp
