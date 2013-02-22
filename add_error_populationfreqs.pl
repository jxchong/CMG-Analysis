#!/usr/bin/env perl
#
# Description: Add systematic error frequency and highest MAF observed in ESP/1000 Genomes to a SSAnnotation.txt file.
#
#
#
# Created by Jessica Chong on 2012-11-20.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $capturearray);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	# 'capture=s' => \$capturearray,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} 
# old code, very slow; reads from multiple files
# elsif (!defined $capturearray) {
	# optionUsage("option --capture not defined\n");
# }

my $commonvarpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
my $errorMAFfile = "$commonvarpath/errorMaxAltAlleleFreq.tsv.gz";

if (!-e $errorMAFfile) {
	die "Cannot read from reference data: $errorMAFfile :$?\n";
}


open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";

if ($inputfile !~ m/SSAnnotation/) {
	die "Input file ($inputfile) isn't an SSAnnotation file\n";
}

my %chr_contents = ();
my $workingchr = "NA";
open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
print OUT join("\t", @header)."\tPrctAltFreqinCMG\tPrctAltFreqinOutsidePop\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($chr, $pos, $vartype, $ref, $alt) = @line[0..4];

	# if ($pos < 207015957) {
	# 	next;
	# }
	# if ($pos > 207304900) {
	# 	exit;
	# }
	
	if ($chr ne $workingchr) {
		my $maxmafs_ref = readData($chr);
		%chr_contents = %{$maxmafs_ref};
		$workingchr = $chr;
		# my $fakelookup = "207304900.T.C";
		# print "static 207304900.T.C: $chr_contents{'207304900.T.C'}{'pop'}\n";
		# print "fakelookup = $fakelookup : $chr_contents{$fakelookup}{'pop'}\n";
		# print "\n";
	}
	
	# print "\npos=$pos, ref=$ref, alt=$alt\n";
	# while (my ($lookup, $val) = each %chr_contents) {
	# 	print "$lookup -> $val ($chr_contents{$lookup}{'pop'})\n";
	# }
	# my $fakelookup = "207304900.T.C";
	# print "1 static 207304900.T.C: $chr_contents{'207304900.T.C'}{'pop'}\n";
	# print "1 fakelookup = $fakelookup : $chr_contents{$fakelookup}{'pop'}\n";
	# print "\n";
	
	my ($errorfreq, $maxpopmaf) = (0, 0);
	my $lookup = "$pos.$ref.$alt";
	if (exists $chr_contents{$lookup}{'pop'}) {
		$errorfreq = $chr_contents{$lookup}{'error'};
		$maxpopmaf = $chr_contents{$lookup}{'pop'};
	}	
	# print "2 lookup = $lookup : $chr_contents{$lookup}{'pop'}\n";
	# print "2 fakelookup = $fakelookup : $chr_contents{$fakelookup}{'pop'}\n";
	# if (defined $chr_contents{$lookup}{'pop'}) {
	# 	print "lookup exists\n";
	# }
	# if (defined $chr_contents{$fakelookup}{'pop'}) {
	# 	print "fake lookup exists\n";
	# }
	
	print OUT join("\t", @line)."\t".sprintf("%.4f", $errorfreq)."\t".sprintf("%.4f", $maxpopmaf)."\n";
}
close FILE;


close OUT;







sub readData {
	my $currchr = $_[0];
	
	my %maxmafs = ();
	my $exit_value;

	print "Reading in error/MAF data for chr $currchr\n";
	my @popdata = `tabix $errorMAFfile $currchr:1-300000000`;
	$exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@popdata\n";
	}
	foreach (@popdata) {
		$_ =~ s/\s+$//;					# Remove line endings
		# print "storing $_\n";
		my ($chr, $varstart, $varend, $ref, $alt, $freqinCMG, $freqinOutside) = split("\t", $_);
		my $lookup = "$varstart.$ref.$alt";
		$maxmafs{$lookup}{'error'} = $freqinCMG;
		$maxmafs{$lookup}{'pop'} = $freqinOutside;
		# print "stored $maxmafs{$lookup}{'pop'} in $lookup\n";
	}
	return \%maxmafs;
}







sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput SSAnnotation file\n";
	print "\t--out\toutput file\n";
	# print "\t--capture\tcapture array used in sequencing (bigexome or v2)\n";
	die;
}