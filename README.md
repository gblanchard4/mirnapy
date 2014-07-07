mirnapy
=======


Gene Blanchard
me@geneblanchard.com

mirna_pipe.py

This script is a python wrapper for the miRNAkey software.
Roy Ronen; Ido Gan; Shira Modai;Alona Sukacheov; Gideon Dror; Eran Halperin; Noam Shomron. miRNAkey: a software for microRNA Deep Sequencing analysis. Bioinformatics 2010; doi: 10.1093/bioinformatics/btq493

This script requires the fastx_toolkit and bwa to be installed on your machine. 

If using ubuntu you can install with 'sudo apt-get install bwa fastx-toolkit', you may also need to install libgtextutils. 

If using a mac I recomend recomend to intall homebrew from 'http://brew.sh/'. Next install the homebrew-science repo with 'brew tap homebrew/science'. 
Install the fastx-toolkit and bwa with 'brew install bwa fastx_toolkit'

An example usage would be:
mirna_pipe.py -i A1.fastq,B1.fastq,A2.fastq,B2.fastq -o A1_B1_A2_B2_output

This would create the A1_B1_A2_B2_output folder in your current working directory.
That folder contains the folowing:
commands.txt: 		A file that lists all commands run
clipped_files: 		The output of the fastx clipping 
alignment_files: 	The output of the BWA alignment
SEQEM: 			The results of the SEQEM command
counts:			Raw counts from the SEQEM command
RSEM:			Groomed tab-seperated files that are ready to input into RSEM
