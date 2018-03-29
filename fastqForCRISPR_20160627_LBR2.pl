#!usr/bin/perl;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use Text::LevenshteinXS qw(distance);
use Getopt::Long;


######### Initializations #####################

my $ind_len = 0;  # length of index
my $bc_len = 6;
my $pat1 = ""; # put the first constant sequence here from the start of the read untill the barcode
my $pat2 = "TAAAGCTG";    # the first 8 nucleotides of the constant
my $kind = "exp";


my $crispr_target = "CTTACCACTTCACCATCGGC";  #gRNA sequence
my $ds_target_1 = "GAATAGTT";    # a  small piece of sequence downstream of the cut site
my $gap_1 = 118;    # gap between $pat2 and $ds_target_1
my $gap_2 = 110;   # gap between $pat2 and $ds_target_2
my $ds_target_2 = "ATTAACCT"; 




my $regex = make_regex ($ind_len, $bc_len, $pat1, $pat2, "exp");
my $regex_1 = $pat2 . "(.+)" . $crispr_target;
my $regex_2 = $pat2 . "(.+)" . $ds_target_1;
my $regex_3 = $pat2 . "(.+)" . $ds_target_2;
my $BC_this;
my $status = "unclear";
my $score = "NA";

#print "$regex\n";


my $fastq_file = $ARGV[0];

die "Cannot open the file $fastq_file\n" unless open (FASTQ , "<" , $fastq_file);



while (my $line = <FASTQ>) {
	
    die "$fastq_file is not a proper fastq file\n" if (substr $line, 0, 1) ne "@";
    $line = <FASTQ>;
    chomp ($line);
    
    if ($line =~ $regex) {
        $BC_this = $1;

        if ($line =~ $regex_1) {
            $status = "wt";
            $score = 0;
        } elsif ($line =~ $regex_2) {
            $score = (length $1) - $gap_1;
            if ($score > 0) {
                $status = "ins";
            } elsif ($score < 0) {
                $status = "del";
            } else {
                $status = "wt_point_mut";
            }
            
        } elsif ($line =~ $regex_3) {
            $score = (length $1) - $gap_2;
            if ($score > 0) {
                $status = "ins";
            } elsif ($score < 0) {
                $status = "del";
            } else {
                $status = "wt_point_mut";
            }

        } else { $status = "not_clear" }
     print $BC_this . "\t" . $status . "\t" . $score . "\t" . $line . "\n";
     # print $BC_this . "\t" . $status . "\t" . $score . "\n";
        
    }
    <FASTQ>;
    <FASTQ>;
    
}


















##########################################################################################
#################***************   SUBROUTINE make_regex  ******************##############
##########################################################################################
# DESCRIPTION:
# A subroutine to make a regular expression string from given arguments for proper matching
# in cases of different formats of read structure
#
# INPUT:
# It takes five arguments
# 	a) $ind_len  -> The length of the index sequence (0 if no index is used)
# 	b) $bc_len   -> The length of the barcode
# 	c) $pat1	 -> The sequence of first constant part ('NA' if no constant part 1 is used)
# 	d) $pat1	 -> The sequence of second constant part (must be DNA string)
# 	e) $kind	 -> Either "exp" (for normalization/expression reads where only barcode is to be extracted)
# 					or "map" (for mapping reads where both the barcode and the neighboring genomic DNA
# 					is to be extracted)
#
# DEPENDS:
# no dependencies
#
# WORKING:
# From each of $pat1 and $pat2, it finds if they are longer than 10 bases. If so, it takes
# only the last 10 bases of $pat2 and the first 10 bases of $pat2 are saved into two additional
# variables $short_pat1 and $short_pat2.
# Next, it determines min and max values in order to put into the regex, as shown below
# ^{ $min1, $max1} $short_pat1 (.{(barcode length -1) , (barcode length + 1)} $short_pat2 .+
# Then it tries to make sure that the $min1 is not a negative number (as it can happen
# if no index is used). Also if it is above 2 (in cases of index absence where $pat1 is >1 longer
# than $short_pat1) it subtracts 2 from min1 to allow for some positional flexibility in
# patten match.
# In the end it makes the regex from all the different components generated.
# If it is for normalization/expression reads then the last .+ does not get parenthesis. e.g.,
# ^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA.+
# If it is for mapping reads then the last .+ gets parthesis around it. e.g.,
# ^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA(.+)
#
# OUTPUT:
# The regular expression string
##########################################################################################

sub make_regex {
	my ($ind_len, $bc_len, $pat1, $pat2, $kind) = @_;
	
	my $len1 = length $pat1;	#get the length of pat1
	my $short_pat1 = $pat1;
	if  ( $len1  > 10 ) {
		$short_pat1 = substr($pat1, ($len1 -10)) ;
	}
	
	my $len2 = length $pat2;	#get the length of pat2
	my $short_pat2 = $pat2;
	if ( $len2  > 10 ) {
		$short_pat2 = substr($pat2, 0, 10) ;
	}
	
	# initialization of some important numbers for generating the final regexp
	my $regex;
	my $min1 = $ind_len + $len1 - (length $short_pat1) - 2;
	my $max1 = $ind_len + $len1 - (length $short_pat1) + 2;
	my $min2 = $bc_len - 1;
	my $max2 = $bc_len + 1;
	
	# dealing with special conditions where no index is used ± no first constant part
	if ($ind_len == 0 ) {
		if ($short_pat1) {
			$min1 = $len1 - length ($short_pat1);
			if ($min1 >= 2 ) { $min1 -= 2 }
			$max1 = $len1 - (length $short_pat1) + 2;
		} else { $min1 = $max1 = 0 }
	}
	
	# now depending on the $kind it makes two different types of regex's
	if ($kind eq "exp" ) {
		$regex = "^.{" . $min1 . "," . $max1 . "}" . $short_pat1 .
        "(.{" . $min2 . "," . $max2 . "})" . $short_pat2 . ".+";
	} elsif ($kind eq "map" ) {
        $regex = "^.{" . $min1 . "," . $max1 . "}" . $short_pat1 .
        "(.{" . $min2 . "," . $max2 . "})" . $short_pat2 . "(.+)";
    }
	
    return qr/$regex/;		# this qr makes a regex object which can be directly used
}

