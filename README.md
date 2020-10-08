# Peptide Alignment

This script written in CUDA C++ calculates the average aligment of petides for each protein in a query set, against a reference seq. The two files should both be in FASTA format and contain amino acids. The output is a tab-delimited file containing the average PAM30 distance for each k-mer in the query set, along with its protein of origin.

## Requirements

The tool requires the NVIDIA C compiler (https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html). The script was tested on CUDA 9.1 but it should work on most versions.

## Installation

To install, simply clone this repository (```git clone github.com/whalleyt/Peptide-Alignment```) and change into the directory. Then run ```make```.

## Usage

After compiling the software one can run it with the ```alignPep``` command. The usage is as follows:


```./alignPep <query_fasta> <reference_fasta> <peptide_length>```

where ```<query_fasta>``` and ```<reference_fasta>``` are both FASTA files containing amino acids in which query refers to the dataset for which to take peptides from; and reference is the dataset for which the query is to be scored against. ```<peplen>``` refers to the length of peptide.


## License and contact

This repository is licensed under the GNU GPL, for more information see the LICENSE file. For any queries either contact through this repository or email whalleyt@cardiff.ac.uk.