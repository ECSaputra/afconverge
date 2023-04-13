# AFconverge

AFconverge (alignment-free converge) is a comparative genomics method that predicts the regulatory adaptations underlying the evolution of convergent phenotypes. Given a set of "comparable" (not necessarily alignable) orthologous sequences of regulatory non-coding elements, AFconverge computes phenotypic convergence signals at individual transcription factor (TF) motif level, identifying motif gains and losses associated with the phenotype.


## Installation and dependencies

AFconverge was developed on Linux and requires the installation of both Python and R. The following lists the dependencies required for running AFconverge:

- Python:

    - Python 3.6.13

    - [PyTorch](https://pytorch.org/) 1.4.0

    - [Pysam](https://pysam.readthedocs.io/en/latest/index.html) 0.15.3

    - pandas 1.1.5
    
    - numpy 1.19.2
    
- R:
    - R 4.2.2
    
    - [RERconverge](https://github.com/nclark-lab/RERconverge)
    
After cloning this Github repository to your machine, add the path to the `afconverge` folder to both the PATH and PYTHONPATH variables in `.bashrc`:


```
export PATH="/path/to/afconverge:$PATH"
export PYTHONPATH="/path/to/afconverge"
```


## Running example analysis

AFconverge requires the following input data:

* A phylogenetic tree file in the Newick format (see example in `data/tree/`)

* FASTA file(s) containing the set(s) of orthologous sequences to be analyzed. 

    - The folder `data/element-fasta/` contains some example FASTA files with the required format. Each file contains the orthologous sequences of a particular element (e.g., the AP3B2 promoter) across species, and each sequence in the file is given a sequence identifier (SeqID) that is the species name (make sure that the species names are ***identical*** to the species names in the tree file).
    
    - As AFconverge is an alignment-free method, sequences in each FASTA files do not have to be aligned or transformed to have the appropriate polarity.
    
* Consensus sequences of transcription factor (TF) motifs in the MEME format (see example in `data/motif-meme/`)

* Phenotype values, represented as an R vector named with the species names (see example in `data/phenotype-values/`). Make sure that the species names are ***identical*** to the species names in the tree file. AFconverge works with both binary and continuous phenotypes.



### Step 1: Quantifying TF binding with motif convolution

As a first step, AFconverge performs one-dimensional convolution with TF motifs as convolutional filters to quantify the binding affinity of each TF on each sequence, estimated as sequence similarity to consensus motifs. **This step only has to be run ONCE for each new set of FASTA files. The motif convolution outputs from this run can then be used for any downstream analysis with any phenotype.**

The function `score_orthologs.py` performs this motif convolution step on the GPU. The following are the required input arguments and their assigned flags:

* `-i, --orthseqfolder`: Folder path containing FASTA files of orthologous sequences

* `-o, --outputfolder`: Folder path to store outputs

* `-m, --meme`: meme file

* `-s, --suffix`: file extension of FASTA file (e.g., .fa, .fasta)

* `-b, --batch`: batch size (default 128)

* `-w, --window`: window size for buffer; must be greater than the longest motif (default 30)

* `-x, --mode`: pooling mode (default max)

The function can be run on the command line as follows:

```
score_orthologs.py -i data/element-fasta -m data/motif-meme/consensus_pwms.meme -o exoutput/motif-scores -s .fa -x max
```
