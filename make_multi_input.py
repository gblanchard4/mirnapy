#! /usr/bin/env python
from optparse import OptionParser
from itertools import islice

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

'''
Python implementation of the broken make_multi_input.pl

output files are:
	_regions.txt
	_reads.txt
	_mappings.txt

'''

def main():

	#Create the argument parser
	parser = OptionParser(usage="Usage: ")

	# db file
	# -d --db
	parser.add_option("-d", "--db", action="store", type="string", dest="db", help="The database file")

	# reads file
	# -r --reads
	parser.add_option("-r", "--reads", action="store", type="string", dest="reads", help="The reads file")

	# bwa output file
	# -b --bwa
	parser.add_option("-b", "--bwa", action="store", type="string", dest="bwa", help="The bwa file")

	# output prefix
	# -o --out 
	parser.add_option("-o", "--output", action="store", type="string", dest="output", help="The output prefix")	

	# Grab command line args
	(options, args) = parser.parse_args()

	# Set argument values
	db = options.db
	reads = options.reads
	bwa = options.bwa
	output = options.output


	# Error checking
	ERROR = False

	# Make sure all arguments have values
	if db == None:
		ERROR = True
	if reads ==  None:
		ERROR = True
	if bwa == None:
		ERROR = True
	if output == None:
		ERROR = True

	# Quit on error
	if ERROR == True:
		sys.exit()

	regions = make_regions_outfile(db, output)
	reads = make_reads_outfile(reads, output)

	make_mappings_out(bwa, output, reads, regions)

	print "Done!"


def make_regions_outfile(db, output):
	regions_dict = dict()
	outfile_name = "%s_regions.txt" % (output)
	regions_list = []
	with open(db, 'r') as fasta:
		for line in fasta:
			# Grab header lines
			if line.startswith('>'):
				region =  line.split(' ')[0].lstrip('>')
				regions_list.append(region)
	with open(outfile_name, 'w') as outfile:
	 	for count, region in enumerate(regions_list):
	 		outstring = "%s %s 0 0\n" % (count, region)
	 		# Write the file
	 		outfile.write(outstring)
	 		# Build the dictionary
	 		regions_dict[region] = count
	return regions_dict

def make_reads_outfile(reads, output):
	reads_dict = dict()
	outfile_name = "%s_reads.txt" % (output)
	reads_linecount = line_count(reads)
	reads_name_linenumber = range(1,reads_linecount, 4)
	reads_name_list = []
	with open(reads, 'r' ) as fastq:
		fastq_lines = fastq.readlines()
	for linenumber in reads_name_linenumber:
		read_name =  fastq_lines[linenumber - 1].lstrip('@').rstrip('\n').split(' ')[0]
		reads_name_list.append(read_name)
	with open(outfile_name, 'w') as outfile:
		for count, name in enumerate(reads_name_list):
			outstring = "%s %s\n" % (count, name)
			# Write the file
			outfile.write(outstring)
			# build the dictionary
			reads_dict[name] = count
	return reads_dict

def make_mappings_out(bwa, output, reads, regions):
	outfile_name = "%s_mappings.txt" % (output)
	with open(bwa, 'r') as sam_file:
		with open(outfile_name, 'w') as outfile:
			for line in sam_file:
				if not line.startswith('@'): 
					split = line.rstrip('\n').split('\t')
					if not (split[1] == '4' or split[1] == '20'):
						number_of_optimal = int(cigar_finder('X0:i:', line))
						#number_of_suboptimal = int(cigar_finder('X1:i:', line))
						if cigar_finder('X1:i:', line) == None:
							number_of_suboptimal = 0
						else:
							number_of_suboptimal = int(cigar_finder('X1:i:', line))
						read_name = split[0]
						edit_distance = cigar_finder('NM:i:', line)
						mirna = split[2]
						
						outstring = "%s %s %s\n" % (reads[read_name], regions[mirna], edit_distance)
						outfile.write(outstring)

						if number_of_optimal > 1 or not number_of_suboptimal == 0:
							try:
								more_mappings_list_raw = cigar_finder('XA:Z:', line).rstrip('\n').split(';')
							except AttributeError:
								print "No proper cigar string found for line: %s" % (line)
								pass

							more_mappings_list = filter(None, more_mappings_list_raw)
							# Sanity check he included
							if not (number_of_optimal + number_of_suboptimal -1) == len(more_mappings_list):
								inconsistant_bwa_line = (number_of_optimal + number_of_suboptimal -1) - len(more_mappings_list)

							# And I quote... get the FRIGGIN-MAPPINGS and write them to file
							for mapped_mirna in more_mappings_list:
								mapped_mirna_values = mapped_mirna.split(',')
								mirna = mapped_mirna_values[0]
								edit_distance = mapped_mirna_values[3]
								

								outstring = "%s %s %s\n" % (reads[read_name], regions[mirna], edit_distance)
								outfile.write(outstring)

def line_count(filename):
	with open(filename) as openfile:
		for index, line in enumerate(openfile):
			pass
	return index + 1

def cigar_finder(cigar, line):
	for element in line.split('\t'):
		if element.startswith(cigar):
			mapping_integer = element.split(':')[-1]
			return mapping_integer




if __name__ == '__main__':
	main()

