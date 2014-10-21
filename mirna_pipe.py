#! /usr/bin/env python

#~~~~~~~~~~~~~~~~~~~~~~
# Imports
#~~~~~~~~~~~~~~~~~~~~~~
import argparse
import os
import commands
import subprocess
import sys
import distutils.spawn
import platform

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

'''
Flow to get to SEQEM fils...

Arguments needed are:
1) Input ` files 
2) Adapter sequence 
3) output directory
4) SEQEM iterations
5) Database to align to 


#fastx clipper
fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -n -v -l 16 -i A3.fastq -o clipped_files/A3.fastq_clipped_sequence.txt

#align
bwa aln -n 2 /home/gene/Code/mirnakey/miRNAkey/DB_mature/mouse/mature_dna_mouse.fa clipped_files/A3.fastq_clipped_sequence.txt > alignment_files/A3.fastq_aligned.sai

#sam_gen_cmd
bwa samse -n 100000 /home/gene/Code/mirnakey/miRNAkey/DB_mature/mouse/mature_dna_mouse.fa alignment_files/A3.fastq_aligned.sai clipped_files/A3.fastq_clipped_sequence.txt > alignment_files/A3.fastq_aligned.sam

#sam_stats
/home/gene/Code/mirnakey/miRNAkey/pipeline/sam_stats.pl new-format < alignment_files/A3.fastq_aligned.sam
A3.fastq - adaptor clipping and alignment to mature_dna_mouse.fa

#SEQEM
/home/gene/Code/mirnakey/miRNAkey/pipeline/make_multi_input.pl /home/gene/Code/mirnakey/miRNAkey/DB_mature/mouse/mature_dna_mouse.fa clipped_files/A3.fastq_clipped_sequence.txt alignment_files/A3.fastq_aligned.sam SEQEM/A3stq_clipped_sequence 1
/home/gene/Code/mirnakey/miRNAkey/pipeline/SEQEM_linux/SEQEM SEQEM/A3stq_clipped_sequence_regions.txt SEQEM/A3stq_clipped_sequence_reads.txt SEQEM/A3stq_clipped_sequence_mappings.txt 1000 SEQEM/A3stq_clipped_sequence_seqem_out.txt 
'''

def main():

	# Build a path from the executed script, Argv[0], to the bin directory
	bin_dir = os.path.realpath(__file__).rstrip('mirna_pipe.py')

	# Default mouse mature database
	#database =  bin_dir+'DB_mature/mouse/mature_dna_mouse.fa'

	# Find out what OS we are running to set the proper SEQEM binary
	if platform.system() == 'Darwin':
		seqem_binary = '%s/SEQEM_MAC' % bin_dir
	else:
		seqem_binary = '%s/SEQEM' % bin_dir

	# Usage string for the help text 
	usage_text = """

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

This script uses some default parameters that can be seen in the options help below. 

"""

	#~~~~~~~~~~~~~~~~~~~~~~
	# Parameters
	#~~~~~~~~~~~~~~~~~~~~~~

	#Create the argument parser
	parser = argparse.ArgumentParser(description=usage_text)

	# input file
	# -i --input
	parser.add_argument("-i", "--input", dest="input", required=True, help="The input fasta, for multiple files seperate with a comma")
	# adapter sequence
	# -a --adapter
	parser.add_argument("-a", "--adapter", dest="adapter", default="TGGAATTCTCGGGTGCCAAGG" , help="The adapter sequence, DEFAULT=TGGAATTCTCGGGTGCCAAGG")
	# output file
	# -o --output
	parser.add_argument("-o", "--output", dest="output", required=True, help="The root output directory")
	# SEQEM iterations 
	# -s --seqem
	parser.add_argument("-s", "--seqem", dest="seqem", default="1000" , help="Seqem iterations, DEFAULT=1000")
	# Alignment Fasta
	# -d --db
	parser.add_argument("-d", "--db", dest="db", required=True, help="The .fa to align to, DEFAULT=mature_dna_mouse.fa")

	# Grab command line args
	args = parser.parse_args()
	# Set argument values
	input_name = args.input
	adapter = str.upper(args.adapter)
	output = args.output
	seqem = args.seqem
	db = args.db


	#~~~~~~~~~~~~~~~~~~~~~~
	# Error checking
	#~~~~~~~~~~~~~~~~~~~~~~
	ERROR = False

	# Make sure the input name is valid
	if input_name == None:
		print 'ERROR: You need to enter an input file or files\n'
		ERROR = True
	else:	
		path_list = [os.path.abspath(fastq) for fastq in input_name.split(',')]

	if output == None:
		print 'ERROR: You need to enter an output directory\n'
		ERROR = True
	else:
		output = os.path.abspath(output)+'/'

	bases = "ACGT"
	if not any(base in adapter for base in bases):
		print "Your adapter sequence contains invalid bases"
	 	ERROR = True

	# Select the mouse mature database if no other option has been supplied
	if db == None:
		db = bin_dir+'DB_mature/mouse/mature_dna_mouse.fa'
	else:
		db = options.db

	# Check if the output directory exists
	if not output == None:
		if not os.path.exists(output):
			os.makedirs(output)

	# Quit on error
	if ERROR == True:
		print "\n Errors found: Use the -h option for more information"
		sys.exit()


	# Run on each file
	master_command_list = []
	# Make output directories
	dirs = ["alignment_files/","clipped_files/","SEQEM/", "counts/", "RSEM/"]
	for directory in dirs:
		master_command_list.append("mkdir -p %s%s" % (output, directory))

	for fastq in path_list:
		 master_command_list.extend(mirna_command_builder(bin_dir, output, adapter, fastq, seqem, db, seqem_binary) )

	shell_file =  output+'commands.txt'
	with open(shell_file, 'wb') as shell:
		# Write the options we used to the file
		input_string =  'Input Files: %s\n' % (input_name)
		adapter_string = 'Adapter: %s\n' % (adapter)
		database_string =  'Database: %s\n' % (db)
		output_string = 'Output %s\n' % (output)
		shell.write(input_string+adapter_string+database_string+output_string)
		for command in master_command_list:
			try:	
				shell.write(command+'\n')
				proc = subprocess.Popen(command, shell=True)
				proc.wait()
			except OSError:
				print 'OH NOES SOMETHING BROKE'

	parse_seqem_counts(output)





def mirna_command_builder(bin_dir, output, adapter, fastq, seqem, db, seqem_binary):
	command_list = []
	fastq_name = fastq.split('/')[-1] 

	# fastx clipper
	fastx_cmd = "fastx_clipper -a %s -n -v -l 16 -i %s -o %sclipped_files/%s_clipped_sequence.txt" % (adapter, fastq, output, fastq_name)
	command_list.append(fastx_cmd)

	# align
	align_cmd = "bwa aln -n 2 %s %sclipped_files/%s_clipped_sequence.txt > %salignment_files/%s_aligned.sai" % (db, output, fastq_name, output, fastq_name)
	command_list.append(align_cmd)

	# sam_gen_cmd
	sam_gen_cmd = "bwa samse -n 100000 %s %salignment_files/%s_aligned.sai %sclipped_files/%s_clipped_sequence.txt > %salignment_files/%s_aligned.sam" % (db, output, fastq_name, output, fastq_name, output, fastq_name)
	command_list.append(sam_gen_cmd)

	# sam_stats
	#sam_stats_cmd = "/home/gene/Code/mirnakey/miRNAkey/pipeline/sam_stats.pl new-format < alignment_files/A3.fastq_aligned.sam"
	#command_list.append(sam_stats_cmd)
	
	# SEQEM
	make_multi_cmd = "%smake_multi_input.py -d %s -r %sclipped_files/%s_clipped_sequence.txt -b %salignment_files/%s_aligned.sam -o %sSEQEM/%s_clipped_sequence " % (bin_dir, db, output, fastq_name, output, fastq_name, output, fastq_name )
	command_list.append(make_multi_cmd)

	seqem_cmd = "%s %sSEQEM/%s_clipped_sequence_regions.txt %sSEQEM/%s_clipped_sequence_reads.txt %sSEQEM/%s_clipped_sequence_mappings.txt %s %scounts/%s_clipped_sequence_seqem_out.txt " % (seqem_binary, output, fastq_name, output, fastq_name, output, fastq_name, seqem, output, fastq_name)
	command_list.append(seqem_cmd)

	return command_list

def parse_seqem_counts(output):
	counts_path = output+"counts/"
	rsem_path =  output+"RSEM/"
	counts_files = os.listdir(counts_path)
	for counts_file in counts_files:
		if counts_file.endswith("_clipped_sequence_seqem_out.txt"):
			name = rsem_path+counts_file.strip("_clipped_sequence_seqem_out.txt")+".tsv"
			with open(counts_path+counts_file, 'r') as counts:
				with open(name, 'w') as outfile:
					outfile.write("miRNA_name\tread_counts\tmiRNA_name\tread_counts\tread_counts\n")
					for line in counts:
						if not line.startswith('#'):
							split_line = line.strip('\n').split(' ')
							mir = split_line[1]
							count = split_line[5]
							outfile.write("%s\t%s\t%s\t%s\t%s\n" %(mir, count, mir, count, count) )

if __name__ == '__main__':
	main()




