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
    
After cloning this Github repository to your machine, add the path to the "afconverge" folder to both the PATH and PYTHONPATH variables:


```
export PATH="/path/to/afconverge:$PATH"
export PYTHONPATH="/path/to/afconverge$PATH"
```


## Running example analysis

AFconverge requires the following input data:

* Fasta file(s) containing the set(s) of orthologous sequences to be analyzed. The folder `data/element-fasta/` 



