#!/usr/bin/env perl
#
# Description: Make a merged systematic error and highest MAF in outside populations file.
#
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

my $databasedir = '/net/grc/vol1/mendelian_projects/mendelian_analysis/software/annovar/humandb';
my @populationinputs = qw(hg19_ASN.sites.2012_04.txt hg19_AFR.sites.2012_04.txt hg19_EUR.sites.2012_04.txt hg19_esp6500si_aa.txt hg19_esp6500si_ea.txt hg19_esp6500si_aa.txt);

my %maxmafs;
foreach my $file (@populationinputs) {
	open (FILE, "$databasedir/$file") or die "Cannot read $databasedir/$file: $!.\n";
	while ( <FILE> ) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $start, $stop, $ref, $alt, $newfreq) = split ("\t", $_);
		if (defined $maxmafs{"$chr.$start.$stop.$ref.$alt"}) {
			if ($newfreq > $maxmafs{"$chr.$start.$stop.$ref.$alt"}) {
				$maxmafs{"$chr.$start.$stop.$ref.$alt"} = $newfreq;
			}
		} else {
			$maxmafs{"$chr.$start.$stop.$ref.$alt"} = $newfreq;
		}
	}
	close FILE;
}



my $errorpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/systematic_error/2012_oct';
my @errorinputs = qw(snv.bigexome.vcf.gz indels.bigexome.vcf.gz snv.v2.vcf.gz indels.v2.vcf.gz);
# indel in systematic error and annovar 1kg file 1       1289367 CT      2       0.82    .

my %errorfreqs;
foreach my $file (@errorinputs) {
	open (FILE, "$errorpath/$file") or die "Cannot read $errorpath/$file: $!.\n";
	while ( <FILE> ) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $pos, $type, $ref, $alt, @errordata) = split("\t", $_);
		my @errorfreqinfo = split("/", $errordata[2]);
		$errorfreqinfo[0] =~ s/OBSERVED=//;
		my $errorfreq = 0;
		if ($errorfreqinfo[1] != 0) {
			$errorfreq = $errorfreqinfo[0]/$errorfreqinfo[1];
		}
		
		if (defined $errorfreqs{"$chr.$pos.$pos.$ref.$alt"}) {
			if ($errorfreqs > $errorfreqs{"$chr.$pos.$pos.$ref.$alt"}) {
				$errorfreqs{"$chr.$pos.$pos.$ref.$alt"} = $newfreq;
			}
		} else {
			$errorfreqs{"$chr.$pos.$pos.$ref.$alt"} = $newfreq;
		}
	}
	close FILE;
}



open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "chr\tstart\tend\tref\talt\tPrctFreqinV2\tPrctFreqinV3\tPrctFreqinCMG\tPrctFreqinOutsidePop\n";




print OUT join("\t", @line)."\t".sprintf("%.4f", $errorfreq*100)."\t".sprintf("%.4f", $maxmaf*100)."\n";
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	print "\t--capture\tcapture array used in sequencing (bigexome or v2)\n";
	die;
}