# AFconverge

AFconverge (alignment-free converge) is a comparative genomics method that predicts the regulatory adaptations underlying the evolution of convergent phenotypes. Given a set of "comparable" (not necessarily alignable) orthologous sequences of regulatory non-coding elements, AFconverge computes phenotypic convergence signals at individual transcription factor (TF) motif level, identifying motif gains and losses associated with the phenotype.


## Installation and dependencies

Running AFconverge requires the installation of both Python and R. The following lists the environment that AFconverge was developed on and the necessary packages:

- Python:

    - Python 3.6.13

    - [PyTorch](https://pytorch.org/) 1.4.0

    - [Pysam](https://pysam.readthedocs.io/en/latest/index.html) 0.15.3

    - pandas 1.1.5
    
    - numpy 1.19.2
    
- R:
    - R 4.2.2
    
    - [RERconverge](https://github.com/nclark-lab/RERconverge)