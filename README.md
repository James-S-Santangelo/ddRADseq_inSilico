# _In silico_ dual digest of DNA sequences
### Author: James S. Santangelo

This script performs a dual digest of DNA sequences contained in a multi-fasta file using two restriction enzymes chosen by the user. The user must provive:

1. The directory containing the fasta file with the sequences to be digested
2. The name of the first restriction enzyme (e.g. EcoRI)
3. The name of the second restriction enzyme (e.g. XbaI)

The script produces two outputs. First, a plot showing the distribution the fragment lengths for fragments that are bound on either side by cuts from a different enzyme. These are the fragments most relevant for ddRADseq studies. Second, the script will ask the user for their desired fragment size range by asking 1) the minimum desired fragment size; 2) the maximum desired fragment size and 3) the size of the bins (in bp) into which fragments should be lumped. It will then output a text file containing summary statistics (e.g. total number of fragments, number of fragments within each size class, etc.)
