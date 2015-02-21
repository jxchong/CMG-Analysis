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
my $errorMAFfile = "$commonvarpath/errorMaxAltAlleleFreq.2015-02-16.tsv.gz";
# my $errorMAFfile = "MYH3.tsv.gz";

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
my $usetabix = 0;

my ($countinputlines, $blah) = split(" ", `wc -l $inputfile`);
if ($countinputlines < 50) {
	$usetabix = 1;
	print "Only $countinputlines variants requested; using tabix for each variant\n";
}

open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
if ($inputfiletype eq 'vcf') {						# skip all the metadata lines at the top of the file
	while (<FILE>) {
		$headerline = $_;
		$headerline =~ s/\s+$//;					# Remove line endings
		if ($headerline !~ /^#CHROM/) {
			if ($headerline =~ "##fileformat=VCFv4") {
				print OUT "$headerline\n";
				print OUT '##INFO=<ID=AFPOP,Number=1,Type=Float,Description="Percent alt allele freq in outbred populations in ExAC, ESP6500, or 1000 Genomes phase 3 (errorMaxAltAlleleFreq.2015-02-16)">'."\n";
				print OUT '##INFO=<ID=AFPOP.DB,Number=1,Type=String,Description="Source database for AFPOP value (highest observed alt allele frequency in outbred population)">'."\n";
				print OUT '##INFO=<ID=AFCMG,Number=1,Type=Float,Description="Percent alt allele freq in exomes sequenced by CMG">'."\n";
				##INFO=<ID=HA,Number=1,Type=Float,Description="AfricanHapMapFreq">
			} else {
				print OUT "$headerline\n";
			}
		} else {
			last;
		}
	}
} else {
	$headerline = <FILE>;
	$headerline =~ s/\s+$//;					# Remove line endings	
}
close FILE;

my $vcfinfocol;
my @header = split("\t", $headerline);
if ($inputfiletype eq 'vcf') {
	for (my $i=0; $i<=$#header; $i++) {
		if ($header[$i] eq 'INFO') {
			$vcfinfocol = $i;
		}
	}
	print OUT join("\t", @header)."\n";
} else {
	print OUT join("\t", @header)."\tPrctAltFreqinCMG\tPrctAltFreqinOutsidePop\n";
}

open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	next if ($_ =~ /^#/);
	my @line = split ("\t", $_);
	my ($chr, $pos, $vartype, $ref, $alt) = @line[0..4];
	# if ($inputfiletype eq 'vcf') {
	# 	$vartype = determinevartype($ref, $alt);
	# }
	
	if ($usetabix == 1) {
		my $maxmafs_ref = readDataTabix($chr, $pos, $pos);
		%chr_contents = %{$maxmafs_ref};
	} elsif ($chr ne $workingchr && $usetabix == 0) {
		my $maxmafs_ref = readData($chr);
		%chr_contents = %{$maxmafs_ref};
		$workingchr = $chr;
	} 
	
	my (@errorfreq, @maxpopmaf, @popmafDB);

	my @altalleles = split(",", $alt);									# in case there are multi-allelic SNPs
	for (my $i=0; $i<=$#altalleles; $i++) {
		my $lookup = "$pos.$ref.$altalleles[$i]";
		if (exists $chr_contents{$lookup}{'pop'}) {
			push(@errorfreq, sprintf("%.4f", $chr_contents{$lookup}{'error'}));
			push(@popmafDB, $chr_contents{$lookup}{'DB'});
			push(@maxpopmaf, sprintf("%.4f", $chr_contents{$lookup}{'pop'}));
		} else {
			push(@errorfreq, '0');
			push(@maxpopmaf, '0');
			push(@popmafDB, '');
		}
	}
	
	my $afpop = join(",", @maxpopmaf);
	my $afpopDB = join(",", @popmafDB);
	my $afcmg = join(",", @errorfreq);
	
	if ($inputfiletype eq 'vcf') {
		print OUT join("\t", @line[0..($vcfinfocol-1)]);
		print OUT "\t$line[$vcfinfocol];AFCMG=$afcmg;AFPOP=$afpop;AFPOP.DB=$afpopDB\t";
		print OUT join("\t", @line[($vcfinfocol+1)..$#line])."\n";
	} else {
		print OUT join("\t", @line)."$afcmg\t$afpop\n";
	}
	
}
close FILE;


close OUT;



# sub determinevartype {
# 	my ($ref, $alt) = @_;
# 	my $vartype = 'SNP';
# 	if (length($ref)>1 || length($alt)>1) {
# 		$vartype = 'indel';
# 	}
# 	return $vartype;
# }

sub readDataTabix {
	my $currchr = $_[0];
	my $startbp = $_[1];
	my $endbp = $_[2];
	
	my %maxmafs = ();
	my $exit_value;

	print "Reading in error/MAF data for chr $currchr\n";
	my @popdata = `tabix $errorMAFfile $currchr:$startbp-$endbp`;
	$exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@popdata\n";
	}
	foreach (@popdata) {
		$_ =~ s/\s+$//;					# Remove line endings
		# print "storing $_\n";
		my ($chr, $varstart, $varend, $ref, $alt, $freqinCMG, $freqinOutside, $DBsourceOutside) = split("\t", $_);
		my $lookup = "$varstart.$ref.$alt";
		$maxmafs{$lookup}{'error'} = $freqinCMG;
		$maxmafs{$lookup}{'pop'} = $freqinOutside;
		$maxmafs{$lookup}{'DB'} = $DBsourceOutside;
		# print "stored $maxmafs{$lookup}{'pop'} in $lookup\n";
	}
	return \%maxmafs;
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
		my ($chr, $varstart, $varend, $ref, $alt, $freqinCMG, $freqinOutside, $DBsourceOutside) = split("\t", $_);
		my $lookup = "$varstart.$ref.$alt";
		$maxmafs{$lookup}{'error'} = $freqinCMG;
		$maxmafs{$lookup}{'pop'} = $freqinOutside;
		$maxmafs{$lookup}{'DB'} = $DBsourceOutside;
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