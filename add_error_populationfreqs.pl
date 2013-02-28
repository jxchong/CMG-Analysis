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
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} 


my $commonvarpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
my $errorMAFfile = "$commonvarpath/errorMaxAltAlleleFreq.tsv.gz";

if (!-e $errorMAFfile) {
	die "Cannot read from reference data: $errorMAFfile :$?\n";
}

my $firstline = `head -1 $inputfile`;
$firstline =~ s/\s+$//;

my $inputfiletype = determineInputType($firstline, $inputfile);

open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";

my %chr_contents = ();
my $headerline;
my $workingchr = "NA";

open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
if ($inputfiletype eq 'vcf') {						# skip all the metadata lines at the top of the file
	while (<FILE>) {
		$headerline = $_;
		$headerline =~ s/\s+$//;					# Remove line endings
		if ($headerline !~ /^#CHROM/) {
			print OUT "$headerline\n";
			next;
		} else {
			last;
		}
	}
} else {
	$headerline = <FILE>;
	$headerline =~ s/\s+$//;					# Remove line endings	
}
close FILE;

my @header = split("\t", $headerline);
print OUT join("\t", @header)."\tPrctAltFreqinCMG\tPrctAltFreqinOutsidePop\n";

open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	next if ($_ =~ /^#/);
	my @line = split ("\t", $_);
	my ($chr, $pos, $vartype, $ref, $alt) = @line[0..4];
	if ($inputfiletype eq 'vcf') {
		$vartype = determinevartype($ref, $alt);
	}
	
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



sub determinevartype {
	my ($ref, $alt) = @_;
	my $vartype = 'SNP';
	if (length($ref)>1 || length($alt)>1) {
		$vartype = 'indel';
	}
	return $vartype;
}

sub determineInputType {
	my ($firstline, $filename) = @_;
	my $inputtype = 'SSAnnotation';
	if ($filename =~ '.vcf' || $firstline =~ 'VCFv4') {
		$inputtype = 'vcf';
	} elsif ($filename =~ 'SeattleSeqAnnotation' || $firstline =~ /^# inDBSNPOrNot/) {
		$inputtype = 'SeattleSeqAnnotation';
	}
	return $inputtype;
}

sub readData {
	my $currchr = $_[0];
	
	my %maxmafs = ();
	my $exit_value;

	print "Reading in error/MAF data for chr $currchr\n";
	my @popdata = `tabix $errorMAFfile $currchr:1-300000000`;
	# my @popdata = `tabix $errorMAFfile $currchr:1-300000`;					# DEBUG
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