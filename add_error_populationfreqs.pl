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
	'capture=s' => \$capturearray,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $capturearray) {
	optionUsage("option --capture not defined\n");
}





open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";

if ($inputfile !~ m/SSAnnotation/) {
	die "Input file ($inputfile) isn't an SSAnnotation file\n";
}

open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
print OUT join("\t", @header)."\tFreqinCMG\tFreqinOutsidePop\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my ($chr, $pos, $vartype, $ref, $alt) = ($line[0], $line[1], $line[2], $line[3], $line[4]);

	my $errorfreq = isSystematicError($chr,$pos,$vartype,$ref,$alt,$capturearray);
	my $maxmaf = isCommonVar($chr,$pos,$vartype,$ref,$alt);

	print OUT join("\t", @line)."\t".sprintf("%.2f", $errorfreq)."\t".sprintf("%.2f", $maxmaf)."\n";
}
close FILE;


close OUT;



sub isCommonVar {
	# remove common variation (ESP and 1kG)
	my ($targetchr, $targetpos, $targettype, $targetref, $targetalt, $mafcutoff) = @_;
	
	my $commonvarpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
	my $thousandgenomesfile = "$commonvarpath/phase1_release_v3.20101123.snps_indels_svs.sites.vcf.gz";
	my $espSNPsfile = "$commonvarpath/ESP6500.snps.vcf.gz";
	my $espindelsfile = "$commonvarpath/esp6500_indels.frq";
		
	my $maxmaf = 0;
	if ($targettype =~ m/snp/i) {
		## to do: implement indel check against ESP when Josh releases the data
		my @varatsamepos = `tabix $espSNPsfile $targetchr:$targetpos-$targetpos`;
		foreach my $variant (@varatsamepos) {
			my ($varchr, $varpos, $rsid, $varref, $varalt, @vardata) = split("\t", $variant);
			$vardata[2] =~ /MAF=((\.|\d|,)+);/;
			my @popMAFs = split(",", $1);
			if ($varpos == $targetpos && $varref eq $targetref && $varalt eq $targetalt) {
				foreach my $varmaf (@popMAFs) {
					my $actualmaf = $varmaf/100;						# in ESP, allele freqs reported as percentages
					if ($actualmaf > $maxmaf) {
						$maxmaf = $actualmaf;
					}
				}
			}
		}
	}
	
	my @varatsamepos = `tabix $thousandgenomesfile $targetchr:$targetpos-$targetpos`;
	foreach my $variant (@varatsamepos) {
		my ($varchr, $varpos, $rsid, $varref, $varalt, @vardata) = split("\t", $variant);
		my @populations = qw(ASN AMR AFR EUR);
		my @popMAFs;
		foreach my $population (@populations) {
			my $maffield = $population."_AF";
			$vardata[2] =~ m/;$maffield=((\.|\d)+)/;
			if (defined $1) {
				push(@popMAFs, $1);
			}
		}
		if ($varpos == $targetpos && $varref eq $targetref && $varalt eq $targetalt) {
			foreach my $varmaf (@popMAFs) {
				if ($varmaf > $maxmaf) {
					$maxmaf = $varmaf;
				}
			}
		}
	}
	
	return $maxmaf;
}

sub isSystematicError {
	# remove systematic errors
	my ($targetchr, $targetpos, $targettype, $targetref, $targetalt, $capturearray) = @_;
		
	my $errorpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/systematic_error/2012_oct';
	my ($snperrorsfile, $indelerrorsfile);
	if ($capturearray eq 'bigexome' || $capturearray eq 'v3') {
		$snperrorsfile = "$errorpath/snv.bigexome.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.bigexome.vcf.gz";	
	} elsif ($capturearray eq 'v2') {
		$snperrorsfile = "$errorpath/snv.v2.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.v2.vcf.gz";
	} else {
		$snperrorsfile = "$errorpath/snv.$capturearray.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.$capturearray.vcf.gz";
	}

	my @errors;
	if ($targettype =~ m/snp/i) {
		@errors = `tabix $snperrorsfile $targetchr:$targetpos-$targetpos`;
	} elsif ($targettype =~ m/indel/i) {
		@errors = `tabix $indelerrorsfile $targetchr:$targetpos-$targetpos`;
	}
	
	my $maxerrorfreq = 0;
	foreach my $error (@errors) {
		my ($errorchr, $errorpos, $errortype, $errorref, $erroralt, @errordata) = split("\t", $error);
		my @errorfreqinfo = split("/", $errordata[2]);
		$errorfreqinfo[0] =~ s/OBSERVED=//;
		my $errorfreq = $errorfreqinfo[0]/$errorfreqinfo[1];
		if ($errorpos == $targetpos && $errorref eq $targetref && $erroralt eq $targetalt) {
			if ($errorfreq >= $maxerrorfreq) {
				$maxerrorfreq = $errorfreq;
			}
		}
	}
	
	return $maxerrorfreq;
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