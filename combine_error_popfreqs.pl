#!/usr/bin/env perl
#
# Description: Make a merged systematic error and highest MAF in outside populations file.
#	Note that Y chromosome calls not available for 1KG or ESP
#
#
# Created by Jessica Chong on 2012-11-20.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $capturearray);

GetOptions(
	'out=s' => \$outputfile,
);


if (!defined $outputfile) {
	optionUsage("option --out not defined\n");
}

my $commonvarpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
my $thousandgenomesfile = "$commonvarpath/phase1_release_v3.20101123.snps_indels_svs.sites.vcf.gz";
my $espSNPsfile = "$commonvarpath/ESP6500.snps.vcf.gz";
my $errorpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/systematic_error/2012_oct';
my @errorinputs = ("$errorpath/snv.bigexome.vcf.gz", "$errorpath/indels.bigexome.vcf.gz", "$errorpath/snv.v2.vcf.gz", "$errorpath/indels.v2.vcf.gz");

# temporary files
# my $commonvarpath = '/nfs/home/jxchong/jessica_annovar';
# my $thousandgenomesfile = "$commonvarpath/1kg.chrXY.vcf.gz";
# my $espSNPsfile = "$commonvarpath/ESP.chrXYsnps.vcf.gz";


open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "#chr\tstart\tend\tref\talt\tPrctFreqinCMG\tPrctFreqinOutsidePop\n";
foreach my $currchr (((1..22), "X", "Y")) {
# foreach my $currchr ((("X", "Y"))) {
	my $maxmafs_ref = readData($currchr);
	my %chr_contents = %{$maxmafs_ref};
	foreach my $pos ( sort {$a<=>$b} keys %chr_contents ) { 
		my $thispos_ref = $chr_contents{$pos};
 		while (my($lookup, $freq_ref) = each %{$thispos_ref}) {
			my ($ref, $alt) = split(/\./, $lookup);
			my $errorfreq = 0;
			if (defined ${$freq_ref}{'error'}) {
				$errorfreq = ${$freq_ref}{'error'};
			}
			my $maxmaf = 0;
			if (defined ${$freq_ref}{'pop'}) {
				$maxmaf = ${$freq_ref}{'pop'};
			}
			print OUT "$currchr\t$pos\t$pos\t$ref\t$alt\t".sprintf("%.4f", $errorfreq*100)."\t".sprintf("%.4f", $maxmaf*100)."\n";	
		}
	}
}
close OUT;








sub readData {
	my $currchr = $_[0];
	
	my %maxmafs;
	my $exit_value;

	print "Reading in ESP data for chr $currchr\n";
	my @espdata = `tabix $espSNPsfile $currchr:1-300000000`;
	$exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@espdata\n";
	}
	foreach (@espdata) {
		$_ =~ s/\s+$//;	
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		$vardata[2] =~ /MAF=((\.|\d|,)+);/;
		my @popMAFs = split(",", $1);
		my $maxpopmaf = 0;
		foreach my $varmaf (@popMAFs) {
			my $actualmaf = $varmaf/100;						# in ESP, allele freqs reported as percentages
			if ($actualmaf > $maxpopmaf) {
				$maxpopmaf = $actualmaf;
			}
		}

		my $lookup = "$ref.$alt";
		if (!defined $maxmafs{$varpos}{$lookup}{'pop'}) {
			$maxmafs{$varpos}{$lookup}{'pop'} = $maxpopmaf;
		} elsif ($maxpopmaf > $maxmafs{$varpos}{$lookup}{'pop'}) {
			$maxmafs{$varpos}{$lookup}{'pop'} = $maxpopmaf;
		}
		# if ($varpos > 1000000) { last; }
	}

	print "Reading in 1KG data for chr $currchr\n";
	my @thousandgenomedata = `tabix $thousandgenomesfile $currchr:1-1000000000`;
	$exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@espdata\n";
	}
	foreach (@thousandgenomedata) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		my @populations = qw(ASN AMR AFR EUR);
		my @popMAFs;
		foreach my $population (@populations) {
			my $maffield = $population."_AF";
			$vardata[2] =~ m/;$maffield=((\.|\d)+)/;
			if (defined $1) {
				push(@popMAFs, $1);
			}
		}
		my $maxpopmaf = 0;
		foreach my $varmaf (@popMAFs) {
			if ($varmaf > $maxpopmaf) {
				$maxpopmaf = $varmaf;
			}
		}

		my $lookup = "$ref.$alt";
		if (!defined $maxmafs{$varpos}{$lookup}{'pop'}) {
			$maxmafs{$varpos}{$lookup}{'pop'} = $maxpopmaf;
		} elsif ($maxpopmaf > $maxmafs{$varpos}{$lookup}{'pop'}) {
			$maxmafs{$varpos}{$lookup}{'pop'} = $maxpopmaf;
		}
		# if ($varpos > 1000000) { last; }
	}

	foreach my $file (@errorinputs) {
		print "Reading in error data for chr $currchr\n";
		my @errordata = `tabix $file $currchr:1-1000000000`;
		$exit_value = $? >> 8;
		if ($exit_value != 0) {
			die "Error: exit value $exit_value\n@espdata\n";
		}
		foreach (@errordata) {
			if ($_ =~ '#') { next; }
			$_ =~ s/\s+$//;					# Remove line endings
			my ($chr, $varpos, $type, $ref, $alt, @errordata) = split("\t", $_);
			my @errorfreqinfo = split("/", $errordata[2]);
			$errorfreqinfo[0] =~ s/OBSERVED=//;
			my $errorfreq = 0;
			if ($errorfreqinfo[1] != 0) {
				$errorfreq = $errorfreqinfo[0]/$errorfreqinfo[1];
			}

			my $lookup = "$ref.$alt";
			if (defined $maxmafs{$varpos}{$lookup}{'error'}) {
				if ($errorfreq > $maxmafs{$varpos}{$lookup}{'error'}) {
					$maxmafs{$varpos}{$lookup}{'error'} = $errorfreq;
				}
			} else {
				$maxmafs{$varpos}{$lookup}{'error'} = $errorfreq;
				$maxmafs{$varpos}{$lookup}{'pop'} = 0;
			}
			# if ($varpos > 1000000) { last; }
		}
	}
	
	return \%maxmafs;
}




sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	print "\t--capture\tcapture array used in sequencing (bigexome or v2)\n";
	die;
}