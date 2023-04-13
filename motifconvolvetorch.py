#!/usr/bin/env python
import argparse
import numpy as np
import os
from pysam import FastaFile
import torch
from torch.nn import functional as F
from torch.utils.data import Dataset
import sys
import pandas as pd

torch.backends.cudnn.deterministic = True

def readbed(filename, up):
    with open(filename,"r") as file:
        data = file.readlines()
    data = [i.split('\t') for i in data]
    chrs = [i[0] for i in data]
    start = [int(i[1]) for i in data]
    end = [int(i[2]) for i in data]
    if(len(data[0])>4): #get the strand
        print("Strand detected")
        up = int(np.floor(up))
        strand = [i[4].strip() for i in data]
        #adjust the regions to acccount for strand and up
        start = [start[i]-up if strand[i]=="+" else start[i] for i in range(len(start))]
        end = [end[i]+up if strand[i]=="-" else end[i] for i in range(len(start))]
    return np.array(chrs), np.array(start, dtype = int), np.array(end, dtype = int)

def countlowercase(arr):
    return sum([1 for c in arr if c.islower()])

def stringstats(string):
    lowercaseratio = countlowercase(string)/len(string)
    string = string.upper()
    tmp = np.array(list(string))
    gccount = 0
    gcpattern = 0
    cgpattern = 0
    if len(tmp) >= 2:
        gccount = np.sum(np.logical_or(tmp == 'C', tmp == 'G'))/len(tmp)
        gcpattern = string.count("GC")/(len(tmp)-1)
        cgpattern = string.count("CG")/(len(tmp)-1)
    return np.array([gccount, gcpattern, cgpattern, lowercaseratio])

def returnonehot(string):
    string = string.upper()
    lookup = {'A':0, 'C':1, 'G':2, 'T':3}
    tmp = np.array(list(string))
    icol = np.where(tmp != 'N')[0]
    out = np.zeros((4,len(tmp)))
    irow = np.array([lookup[i] for i in tmp[icol]])
    if len(icol)>0:
        out[irow,icol] = 1
    return np.asarray(out)

def write_output(filename, mat, names):
    with open(filename, "w") as file:
        for i in range(len(names)):
            file.write(names[i])
            file.write("\t")
        file.write("GC_ratio\t")
        file.write("GC_pattern\t")
        file.write("CG_pattern\t")
        file.write("Masked_ratio\n")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                file.write("{:.4f}".format(mat[i,j]))
                if j != mat.shape[1]-1:
                    file.write("\t")
            file.write("\n")

def createmasks(len1, len2, mat):
    for i1, l1 in enumerate(len1):
        for i2, l2 in enumerate(len2):
            mat[i1, i2, int(np.abs(l1-l2+2)):, 0] = 0

class MEME():
    def __init__(self):
        self.version = 0
        self.alphabet = ""
        self.strands = ""
        self.headers = []
        self.background = []
        self.names = []
        self.nmotifs = 0

    def parse(self, text, transform, kernel_mode):
        with open(text,'r') as file:
            data = file.read()
        data = data.split("\n\n")
        data = data[:-1]

        offset_metadata = 4
        self.nmotifs = len(data) - offset_metadata
        self.version = int(data[0].split(' ')[-1])
        self.alphabet = data[1][10:].strip()
        self.strands = data[2][9:].strip()
        self.background = np.array(data[3].split('\n')[1].split(' ')[1::2],dtype=float)

        out_channels = self.nmotifs * 2

        lens = np.array([len(i.split('\n')[2:]) for i in data[offset_metadata:]])
        height = np.max(lens)
        maximumpadding = height - np.min(lens)
        width = len(self.alphabet)
        out = np.zeros((out_channels, width, height))

        data = data[offset_metadata:]
        for k, i in enumerate(data):
            tmp = i.split('\n')
            self.names.append(tmp[0].split()[-1])
            self.headers.append('\n'.join(tmp[:2]))
            kernel = np.array([j.split() for j in tmp[2:]],dtype=float).T
            if (transform == "constant"):
                bg=np.repeat(0.25,width).reshape(1,width)
            if (transform == "local"):
                bg=np.average(kernel,0).reshape(1,width)
            if (transform != "none"):
                offset=np.min(kernel[kernel>0])
                kernel=np.log((kernel+offset)/bg)
            if kernel_mode=='PPM':
                out[2*k  , :, :kernel.shape[1]] = kernel
                out[2*k+1, :, :kernel.shape[1]] = kernel[::-1, ::-1]
            elif kernel_mode=='IC':
                kernel = (kernel + 0.00001)/(1+4*0.00001) ### from this line onwards, convert PPM to PWM
                info_content = -kernel*np.log2(kernel)
                info_content_total = np.log2(4) ### DNA motifs have alphabet length of 4 (ACGT)
                uncertainty = np.sum(info_content, axis=0)
                ic_position = info_content_total - uncertainty
                kernel_out = 0*kernel
                for idx in range(kernel.shape[1]):
                    kernel_out[:,idx] = ic_position[idx]*kernel[:,idx]
                out[2*k  , :, :kernel.shape[1]] = kernel_out
                out[2*k+1, :, :kernel.shape[1]] = kernel_out[::-1, ::-1]
        return torch.from_numpy(out)

class SegmentData:
    def __init__(self, bed, batchsize, genome, windowsize, up):
        self.chrs, self.starts, self.ends = readbed(bed, up)
        self.midpoints = np.asarray(np.ceil((self.starts + self.ends)/2),dtype=int)
        self.starts = self.midpoints - windowsize
        self.ends = self.midpoints + windowsize
        self.batchsize = batchsize
        self.n = len(self.chrs)
        self.seqs = FastaFile(genome)
        self.padding = windowsize
        refs = self.seqs.references
        lengths = self.seqs.lengths
        self.limits = {refs[i]: lengths[i] for i in range(len(refs))}
        self.out = open("coordinatesUsed.bed", "w")

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def __getitem__(self, i):
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        height = np.max(self.ends[i1:i2] - self.starts[i1:i2]) + self.padding
        width = 4
        batch = np.zeros((batchsize, width, height)) 
        stats = np.empty((batchsize, 4))
        for i, c, s, e in zip(range(i2-i1), self.chrs[i1:i2], self.starts[i1:i2], self.ends[i1:i2]):
            self.out.write(c+"\t"+str(s)+"\t"+str(e)+"\n")
            if s>0 and e<self.limits[c]:
                seg = self.seqs.fetch(c, s, e)
            else:
                seg = "N"*(self.padding*2)
            stats[i] = stringstats(seg)
            batch[i, :, :(e-s)] = returnonehot(seg)
        return torch.from_numpy(batch), stats

    def __del__(self):
        self.out.close()


class OrthologousSequences:
    def __init__(self, orth_seqs, batchsize, windowsize):
        self.seqs = FastaFile(orth_seqs)
        self.padding = windowsize
        self.batchsize = batchsize
        self.labels = self.seqs.references
        self.n = len(self.labels)
        seq_lengths = np.zeros(self.n).astype(int)
        for j in range(len(seq_lengths)):
            seq = self.seqs.fetch(self.labels[j])
            seq_lengths[j] = len(seq)
        self.seq_lengths = seq_lengths
    
    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))
    
    def __getitem__(self, i):
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        refs_batch = self.labels[i1:i2]
        seq_lengths_batch = np.array(self.seq_lengths[i1:i2])
        height = np.max(np.array(seq_lengths_batch)) + self.padding
        width = 4
        batch = np.zeros((batchsize, width, height))
        stats = np.empty((batchsize, 4))
        for j in range(len(refs_batch)):
            stats[j] = stringstats(self.seqs.fetch(refs_batch[j]))
            batch[j, :, :int(seq_lengths_batch[j])] = returnonehot(self.seqs.fetch(refs_batch[j]))
        return torch.from_numpy(batch), stats


class OrthologousSequences_v1:
    def __init__(self, orth_seqs, batchsize, windowsize):
        seqs_info = pd.read_csv(orth_seqs, sep='\t')
        idx_na = seqs_info.loc[pd.isna(seqs_info.sequence),:].index
        seqs_info.drop(idx_na, axis=0, inplace=True)
        self.seqs_info = seqs_info
        self.padding = windowsize
        self.batchsize = batchsize
        self.n = len(self.seqs_info.index)
        seq_lengths = np.zeros(self.n).astype(int)
        for j in range(len(seq_lengths)):
            seq = self.seqs_info.iloc[j].sequence
            seq_lengths[j] = len(seq)
        self.seq_lengths = seq_lengths

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def __getitem__(self, i):
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        height = np.max(np.array(self.seq_lengths[i1:i2])) + self.padding
        width = 4
        seqs_batch = self.seqs_info.iloc[i1:i2,:]
        seq_lengths_batch = np.array(self.seq_lengths[i1:i2])
        batch = np.zeros((batchsize, width, height))
        stats = np.empty((batchsize, 4))
        for j in range(len(batch)):
            stats[j] = stringstats(seqs_batch.sequence.iloc[j])
            batch[j, :, :int(seq_lengths_batch[j])] = returnonehot(seqs_batch.sequence.iloc[j])
        return torch.from_numpy(batch), stats


def motif_convolve(orth_seqs, meme, batch=128, window=30, transform='none', kernel_mode='IC', mode='max'):
    motif = MEME()
    kernels = motif.parse(meme, transform, kernel_mode)
    sequences = OrthologousSequences(orth_seqs, batch, window)
    out = np.empty((sequences.n, motif.nmotifs+4))
    print("Calculating convolutions")
    for i in range(len(sequences)):
        i1, i2 = i*batch, (i+1)*batch
        if i2 >= sequences.n: i2 = sequences.n
        mat, out[i1:i2, motif.nmotifs:] = sequences[i]
        tmp = F.conv1d(mat, kernels)
        if mode == "average":
            tmp = F.avg_pool1d(tmp, tmp.shape[2]).numpy()
        if mode == "max":
            tmp = F.max_pool1d(tmp, tmp.shape[2]).numpy()
        out[i1:i2, :motif.nmotifs] = np.max(tmp.reshape(tmp.shape[0],-1,2), axis=2)
    out = pd.DataFrame(out)
    colnames = []
    for i in range(len(motif.names)):
        colnames.append(motif.names[i])
    colnames.append('GC_ratio')
    colnames.append('GC_pattern')
    colnames.append('CG_pattern')
    colnames.append('Masked_ratio')
    out.columns=colnames
    #out.index = sequences.seqs_info.sequence_fa.replace('.fa', '')
    out.index = sequences.labels
    return(out)


def motif_convolve_with_pval(orth_seqs, meme, null_motif_scores, batch=128, window=240, transform='none', kernel_mode='IC', mode='max', compute_pval=False):
    motif = MEME()
    kernels = motif.parse(meme, transform, kernel_mode)
    sequences = OrthologousSequences(orth_seqs, batch, window)
    out_scores = np.zeros((sequences.n, motif.nmotifs+4))
    out_pval = np.zeros((sequences.n, motif.nmotifs+4))
    print("Calculating convolutions")
    for i in range(len(sequences)):
        i1, i2 = i*batch, (i+1)*batch
        if i2 >= sequences.n:i2 = sequences.n
        mat, out_scores[i1:i2, motif.nmotifs:] = sequences[i]
        tmp = F.conv1d(mat, kernels)
        if mode == "average":
            tmp = F.avg_pool1d(tmp, tmp.shape[2]).numpy()
        if mode == "max":
            tmp = F.max_pool1d(tmp, tmp.shape[2]).numpy()
        out_scores[i1:i2, :motif.nmotifs] = np.max(tmp.reshape(tmp.shape[0],-1,2), axis=2)
        if compute_pval==True:
            idx_complement = np.argmax(tmp.reshape(tmp.shape[0],-1,2), axis=2)
            for r in range(idx_complement.shape[0]):
                for c in range(idx_complement.shape[1]):
                    null_scores_c = null_motif_scores[:, (c*2)+idx_complement[r,c]]
                    percentile = getPercentile(null_scores_c, out_scores[i1+r,c])
                    out_pval[i1+r,c] = 1-percentile
    out_scores = pd.DataFrame(out_scores)
    #out_pval = pd.DataFrame(out_pval)
    colnames = []
    for i in range(len(motif.names)):
        colnames.append(motif.names[i])
    colnames.append('GC_ratio')
    colnames.append('GC_pattern')
    colnames.append('CG_pattern')
    colnames.append('Masked_ratio')
    out_scores.columns=colnames
    out_scores.index = sequences.seqs_info.sequence_fa.replace('.fa', '')
    if compute_pval==False:
        return out_scores
    elif compute_pval==True:
        out_pval = pd.DataFrame(out_pval)
        out_pval.columns=colnames
        out_pval.index = sequences.seqs_info.sequence_fa.replace('.fa', '')
        #out_scores.index = sequences.seqs_info.sequence_fa.replace('.fa', '')
        out_pval.index = sequences.seqs_info.sequence_fa.replace('.fa', '')
        return out_scores, out_pval


def simulateKmers(sequence_len, nucleotides, background):
    random_numbers = np.random.uniform(0,1,sequence_len)
    background_bin_boundaries = np.array([background[0], np.sum(background[:2]), np.sum(background[:3]),np.sum(background[:4])])
    kmer = np.empty(sequence_len, dtype=str)
    for i in range(len(random_numbers)):
        if random_numbers[i] <= background_bin_boundaries[0]:
            kmer[i] = nucleotides[0]
        if random_numbers[i] <= background_bin_boundaries[1] and random_numbers[i] > background_bin_boundaries[0]:
            kmer[i] = nucleotides[1]
        if random_numbers[i] <= background_bin_boundaries[2] and random_numbers[i] > background_bin_boundaries[1]:
            kmer[i] = nucleotides[2]
        if random_numbers[i] >= background_bin_boundaries[2]:
            kmer[i] = nucleotides[3]
    return np.str.join('',kmer)

class NullSequences:
    def __init__(self, sequences, batchsize, max_len):
        self.sequences = sequences
        self.max_len = max_len
        self.batchsize = batchsize
        self.n = len(sequences)

    def __len__(self):
        return int(np.ceil(self.n / self.batchsize))

    def __getitem__(self, i):
        i1, i2 = i*self.batchsize, (i+1)*self.batchsize
        if i2 >= self.n: i2 = self.n
        batchsize = int(i2 - i1)
        height = self.max_len
        width = 4
        seqs_batch = self.sequences[i1:i2]
        batch = np.zeros((batchsize, width, height))
        for j in range(len(batch)):
            batch[j, :, :] = returnonehot(seqs_batch[j])
        return torch.from_numpy(batch)

def compute_null_motif_score_distributions(meme, numnulls, nucleotides, background, batchsize=128, transform='none', kernel_mode='IC'):
    motif = MEME()
    kernels = motif.parse(meme, transform=transform, kernel_mode = kernel_mode)
    # get unique lengths across all motifs
    kernels_colsums = torch.sum(kernels, axis=1)
    kernels_nonzero = kernels_colsums > 0
    kernel_lens = np.sum(kernels_nonzero.numpy().astype(int), axis=1)
    kernel_lens = np.unique(kernel_lens)
    # simulate Kmers
    max_kmer_len = np.max(kernel_lens)
    simulated_kmers_maxlen = np.array([], dtype=str)
    for j in range(numnulls):
        simulated_sequence = simulateKmers(max_kmer_len, nucleotides, background)
        simulated_kmers_maxlen = np.append(simulated_kmers_maxlen, simulated_sequence)
    kmers = NullSequences(simulated_kmers_maxlen, batchsize, max_kmer_len)
    out = torch.zeros((numnulls, len(kernels)))
    for i in range(len(kmers)):
        i1, i2 = i*batchsize, (i+1)*batchsize
        mat = kmers[i]
        tmp = F.conv1d(mat, kernels)
        out[i1:i2,:] = tmp[:,:,0]
    return out

def getPercentile(distribution, value):
    percentile = (distribution<=value).numpy().astype(int).mean()
    return percentile

