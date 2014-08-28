# Getting a list of SNPs from Cortex

This package runs the cortex con assembler on FASTQ input files from a test and a reference genome. The instructions here will allow you to replicate our analysis. 

## Setting up

You will need:

1. Cortex compiled at k = 95 in the `$PATH` on your system - available on [SourceForge](http://sourceforge.net/projects/cortexassembler/files/cortex_con/cortex_con_beta_0.04c.zip/download)
2. NCBI BLAST+ in the `$PATH` on your system - available from [NCBI](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
3. The Ruby `rake` gem, and Ruby > 2.0.0 - install `rake` with `gem install rake`



The `cortex_compare` is a self-contained directory (i.e you can move and copy it wherever you want) that contains a `Rakefile` (kinda like a makefile) that will automate the tasks to generate the list of SNPs found by Cortex, as done in our manuscript. 

Check that the package can find your binaries with `rake check_kit`. It will complain if it can't find anything.

We always use Paired End reads, so we assume that you have a left and right read file for both test and control sequencing experiments. Add those (or link them to the `cortex_compare`)

## Creating the SNP list
1. `cd cortex_compare`
2. Add left and right fastq files from the test, rename to `test_left.fq` and `test_right.fq`
3. Add left and right fastq files from the control, rename to `control_left.fq` and `control_right.fq` 
4. Kick off the whole process with `rake make_snp_list` which will run cortex, get the reference genome (Arabidopsis TAIR10) over the internet, run BLAST using the cortex generated contigs against the reference and from this generate the list of SNPs cortex found.

## Running specific parts of the pipeline
Run `rake tasks` to see a list of tasks in the `Rakefile`. The two most important are `rake bubbles.fasta` which generates cortex output and `rake cortex_bubbles_matched_to_snps.txt` which matches cortex output to the reference and creates SNP lists.

## Starting from the start
Running `rake clean` will delete all intermediate files and databases, leaving only your reads and the `bubbles.fasta` and `cortex_bubbles_matched_to_snps.txt` results. 

