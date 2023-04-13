#!/net/dali/home/mscbio/els191/anaconda3/envs/afconverge/bin/python
import sys
from os import listdir
from os.path import isfile, join
import numpy as np
import argparse
from motifconvolvetorch import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--orthseqfolder', required=True, help='folder containing FASTA files of orthologous sequences')
parser.add_argument('-o', '--outputfolder', required=True, help='output folder')
parser.add_argument('-m', '--meme', required=True, help='meme file')
parser.add_argument('-s', '--suffix', required=True, help='suffix of fasta file (example: .fa)')
parser.add_argument('-b', '--batch', required=False, nargs='?', const=128, default=128, type=int, help='batch size')
parser.add_argument('-w', '--window', required=False, nargs='?', const=30, default=30, type=int, help='window size')
parser.add_argument('-x', '--mode', required=True, nargs='?', const='max', default='max', type=str, help='pool mode')

args = parser.parse_args()
orth_seq_folder = args.orthseqfolder
output_folder = args.outputfolder
meme = args.meme
suffix = args.suffix
batch_size = args.batch
window_size = args.window
mode = args.mode

all_files = [f for f in listdir(orth_seq_folder) if np.logical_and(isfile(join(orth_seq_folder,f)), '.fai' not in f)]

for i in range(len(all_files)):
    element_i = all_files[i].replace(suffix, '')
    print(i+1,'/',len(all_files), element_i)
    outputfile = join(output_folder, 'motif_scores_'+element_i+'.bed')
    orth_seqs = join(orth_seq_folder, all_files[i])
    scores = motif_convolve(orth_seqs, meme, mode=mode, batch=batch_size, window=window_size)
    scores.to_csv(outputfile, sep='\t', header=True, index=True)


