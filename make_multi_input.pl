#!/usr/bin/perl
use warnings;
use strict;

# This script receives as input: 
#
# 1) list of all genes/ regions that a mapped file can map to.
#    If the mapping was to mirBase - then a list of all mir names in mirBase.
#	 If the mapping was to a genome, this requirs some work - get a list of all known genes/regions in the genome. 
#	 Anyhow - basically a FASTA file - all headers will be taken and numnered. 
#
# 2) Reads file. This is the s_2_sequence.txt FASTQ file that Illumina GA outputs (befor/after clipping, no matter)
#
# 3) Output of BWA mapping of 2 to 1 (generation done using the samse -n option!). 


# It will then output all 3 files needed to run multiread! yey. 

###########################################################################################################
################################## Command line args & Global vars ########################################
###########################################################################################################

if (scalar @ARGV != 4 and scalar @ARGV != 5){
	#print join " ", @ARGV;
	print "\nUsage: make_multi_input.pl <genes_or_regions_file> <reads_file> <bwa_output_file> <out-prefix> [new_bwa_format? 0/1]\n";
	print "\n\tOutput files will be: out-prefix_regions.txt, out-prefix_reads.txt and out-prefix_mappings.txt\n";
	die "\tNote: BWA aligner must be run with the 'samse -n' option in order for this script to run properly.\n\n";
	exit 1;
}

my ($regions_filename, $reads_filename, $bwa_out_filename, $out_prefix, $bwa_new_format) = @ARGV;

# global vars
my $regions_out_filenmae = $out_prefix."_regions.txt";
my $reads_out_filename = $out_prefix."_reads.txt";
my $mappings_out_filename = $out_prefix."_mappings.txt";
my %regions_2_regions_nums;
my %reads_2_read_nums;
my $show_inconsistant_mappings_details = 0;
#my $format;

###########################################################################################################
########################################## Main Program ###################################################
###########################################################################################################

# check all input files open and exits
open (my $REGIONS_INPUT,'<', $regions_filename) or die "can't open $regions_filename. $!\n";
#open (my $READS_INPUT  ,'<', $reads_filename) or die "can't open $reads_filename. $!\n";
#$format = check_file_format($READS_INPUT);
#close $READS_INPUT;

open (my $READS_INPUT  ,'<', $reads_filename) or die "can't open $reads_filename. $!\n";
open (my $BWA_INPUT    ,'<', $bwa_out_filename) or die "can't open $bwa_out_filename. $!\n";

# open all output files, regions/genes, reads and mappings
open (my $REGIONS_OUT ,'>', $regions_out_filenmae) or die "can't open out file $regions_out_filenmae. $!\n";
open(my $READS_OUT    ,'>', $reads_out_filename) or die "can't open out file $reads_out_filename. $!\n";
open(my $MAPPINGS_OUT ,'>', $mappings_out_filename) or die "can't open out file $mappings_out_filename. $!\n";

########################## make regions out file ############################
my $region_num = 0;
while (my $region_line = <$REGIONS_INPUT>){
	if ($region_line =~ /^>/){
		# this is a fasta header, make it a region name
		chomp $region_line;
		my @region_line_arr = split (/\s+/, $region_line);
		my $region_name = substr($region_line_arr[0], 1, (length $region_line_arr[0]) -1);
		print $REGIONS_OUT "$region_num $region_name 0 0\n";
		$regions_2_regions_nums{$region_name} = $region_num;
		$region_num++;
	}
}
close ($REGIONS_OUT);

#########################  make reads out file #############################

my $read_num = 0;
while (my $read_line = <$READS_INPUT>){
	if ($read_line =~ /^[@|>]/){
		# this is a fastq header, make it a read name
		my $read_name = substr($read_line, 1, length($read_line) - 15);
		print $READS_OUT "$read_num ";
		print $READS_OUT $read_name;
		print $READS_OUT "\n";
		$reads_2_read_nums{$read_name} = $read_num;
		$read_num++;
	}
}
close $READS_OUT;

######################### make mappings out file ###########################

my $inconsist_count = 0;
my $inconsist_deg_sum = 0;
if ($bwa_new_format){
	
	# read new format SAM file
	while(my $line = <$BWA_INPUT>){
		 
		 next if ($line =~/^@/); # header line
		 my @line_arr = split(/\s+/,$line);
		 next if ($line_arr[1] eq "4" or $line_arr[1] eq "20"); # non mapped line
		 my $inconsist_deg = read_mapping_line_and_generate_new_format($line);
		 if($inconsist_deg != 0){
		 	$inconsist_count ++;
		 	$inconsist_deg_sum += $inconsist_deg;
		 }
		 
	}
}
else{
	
	# read old format SAM file
	while (my $line = <$BWA_INPUT>){		
		# supposed to only read lines that start with '>' in this while, problem if not
		if ($line !~ /^>/){
			die "[make_multi_input.pl] Mapping-line (expecting only mapping headers). Quitting...\n";
		}
		my @line_arr = split(/\s+/, $line);
		my $curr_read_name = substr($line_arr[0], 1, (length $line_arr[0]) -1);
		read_mappings_and_generate($BWA_INPUT, $line_arr[1], $curr_read_name);
	}
}

if ($show_inconsistant_mappings_details){
	if ($inconsist_count > 0){
		print "[make_multi_input.pl] Encountered $inconsist_count inconsistent BWA mappings".
					" (#optimal + #subopt -1 != #mappings in XA-tag). Average diff: ". $inconsist_deg_sum / $inconsist_count ."\n";	
	}	
}

close $MAPPINGS_OUT;

######################### close all input files ############################
close $READS_INPUT;
close $REGIONS_INPUT;
close $BWA_INPUT;

#####################################################################################################
##################################### Auxillery Sub Routines ########################################
#####################################################################################################
 
sub read_mappings_and_generate {
	
	my ($BWA_IN_FILE, $number_of_mappings, $curr_read_name) = @_;
	foreach (1..$number_of_mappings) {
		my $map_line = <$BWA_IN_FILE>;
		my @map_line_arr = split(/\s+/, $map_line);
		
		# DEBUG stuff
		#print "curr_read_name is: $curr_read_name\n";
		#print "number for that read: $reads_2_read_nums{$curr_read_name}\n";
		#print "region is $map_line_arr[0]\n";
		#print "number for that region is: $regions_2_regions_nums{$map_line_arr[0]}\n";
		#print "edit dist is: $map_line_arr[1]\n";
		
		print $MAPPINGS_OUT "$reads_2_read_nums{$curr_read_name} $regions_2_regions_nums{$map_line_arr[0]} $map_line_arr[2]\n";
	}
}

sub read_mapping_line_and_generate_new_format{
	
	 my ($line) = @_;
	 my @line_arr = split(/\s+/, $line);
	 my $how_inconsistant_bwa_line = 0;
	 
	 # get mapping details
	 $line =~ /X0:i:([0-9]+)/;
	 my $num_of_opt = $1;
	 $line =~ /X1:i:([0-9]+)/;
	 my $num_of_subopt = $1;
	 
	 # get the mappings
	 my $curr_read_name = $line_arr[0];
	 
	 # first (1st optimal) mapping
	 my $mapped_to = $line_arr[2];
	 $line =~ /NM:i:([0-9]+)/;
	 my $edit_dist = $1;
	 
	 print $MAPPINGS_OUT "$reads_2_read_nums{$curr_read_name} $regions_2_regions_nums{$mapped_to} $edit_dist\n";
	 
	 # get more mappings if needed
	 if ($num_of_opt > 1 or $num_of_subopt != 0){
	 	$line =~ /XA:Z:(.+)/;
	 	my $more_mappings = $1;
	 	chomp $more_mappings;
	 	my @more_mappings_list = split(/;/,$more_mappings);
	 	
	 	# sanity check
	 	if ( ($num_of_opt + $num_of_subopt -1 ) != scalar @more_mappings_list){
	 		$how_inconsistant_bwa_line = ( ($num_of_opt + $num_of_subopt -1 ) - (scalar @more_mappings_list) );
	 	}
	 	
	 	foreach my $mapping_inf (@more_mappings_list){
	 		# get the FRIGGIN-MAPPINGS and write them to file
	 		my @mapping_inf_list = split(/,/,$mapping_inf);
	 		$mapped_to = $mapping_inf_list[0];
	 		$edit_dist = $mapping_inf_list[3];
	 		print $MAPPINGS_OUT "$reads_2_read_nums{$curr_read_name} $regions_2_regions_nums{$mapped_to} $edit_dist\n";
	 	}
	 }
	 return $how_inconsistant_bwa_line;
}

sub check_file_format{
	
	my ($FH) = @_;
	my $ans = "";
	my $line1 = <$FH>;
	
	if ($line1 =~ /^@/){
		# FASTQ
		return "Q";
	}
	elsif($line1 =~ /^>/){
		# FASTA
		return "A";
	}
	
	# not FASTA or FASTQ
	return "";
}
