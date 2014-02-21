2kplus2
=======

2kplus2 is a c++ source code used for the detection and the classification of single nucleotide polymorphisms in transformed De Bruijn graphs using Cortex assembler.

Introduction

2kplus2 is a new algorithm that search graphs produced from an in-house assembler called cortex [1]. The Complex Bubble Detection (CBD) algorithm use concepts taken from the graph theory in order to search for possible SNPs. Cortex uses de bruijn graphs to assemble genomes from sequencing reads, and in the process, an interesting phenomenon emerges from these graphs in forms of "bubbles" [1]. These bubbles maybe a formulation of potential SNPs and Indels which, if found in the graph, can be extracted for further analysis [2]. The proposed algorithm construct a simplied undirected graph from the cortex's de bruijn graph then search for specic bubbles that match 2k+2 cycles where K is the k-mer size of the assembled genome.  2K+2 cycles have high prediction rate to be potential SNPs.



* [1] De novo assembly and genotyping of variants using colored de Bruijn graphs. Iqbal Z, Caccamo M, Turner I, Flicek P, McVean G. Nature Genetics. 2012 Feb; 44(2):226-32
* [2] Identifying and classifying trait linked polymorphisms in non reference species by walking coloured de Bruijn graphs Leggett R, et al, plos one, 2013.
