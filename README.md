2kplus2
=======

2kplus2 is a c++ source code used for the detection and the classification of single nucleotide polymorphisms in transformed De Bruijn graphs.

Introduction

2kplus2 is a new algorithm that search graphs produced from an in-house assembler called cortex [ref]. The Complex Bubble Detection (CBD) algorithm use concepts taken from the graph theory in order to search for possible SNPs. Cortex uses de bruijn graphs to assemble genomes from sequencing reads, and in the process, an interesting phenomenon emerges from these graphs in forms of "bubbles" [ref]. These bubbles maybe a formulation of potential SNPs and Indels which, if found in the graph, can be extracted for further analysis [ref]. The proposed algorithm construct a simplied graph from the cortex's de bruijn graph then search for specic bubbles that t the description of a potential SNP.


