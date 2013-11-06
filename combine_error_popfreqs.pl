#!/usr/bin/env perl
#
# Description: Make a merged systematic error and highest altAF in outside populations file.
#	Note that Y chromosome calls not available for 1KG or ESP
#
#
# Created by Jessica Chong on 2012-11-20.


use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum);

my ($inputfile, $outputfile, $capturearray, $tempstartbp, $tempendbp);

GetOptions(
	'out=s' => \$outputfile,
	# 's:i' => \$tempstartbp,
	# 'e:i' => \$tempendbp,
);


if (!defined $outputfile) {
	optionUsage("option --out not defined\n");
}

my $commonvarpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
my $thousandgenomesfile = "$commonvarpath/phase1_release_v3.20101123.snps_indels_svs.sites.vcf.gz";
my $errorpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/systematic_error/2013_april';
my @errorinputs = ("$errorpath/bigexome.vcf.gz", "$errorpath/v2.vcf.gz");

# temporary files
# my $commonvarpath = '/nfs/home/jxchong/jessica_annovar';
# my $thousandgenomesfile = "$commonvarpath/1kg.chrXY.vcf.gz";
# my $espSNPsfile = "$commonvarpath/ESP.chrXYsnps.vcf.gz";


open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "#chr\tstart\tend\tref\talt\tPrctAltFreqinCMG\tPrctAltFreqinOutsidePop\n";
# foreach my $currchr (((1..22), "X", "Y")) {
foreach my $currchr (("X", "Y")) {
	print "Reading in chr $currchr data\n";
	my $maxaltAFs_ref = readData($currchr);
	my %chr_contents = %{$maxaltAFs_ref};
	foreach my $pos ( sort {$a<=>$b} keys %chr_contents ) { 
		my $thispos_ref = $chr_contents{$pos};
 		while (my($lookup, $freq_ref) = each %{$thispos_ref}) {
			my ($ref, $alt) = split(/\./, $lookup);
			my $errorfreq = 0;
			if (defined ${$freq_ref}{'error'}) {
				$errorfreq = ${$freq_ref}{'error'};
			}
			my $maxaltAF = 0;
			if (defined ${$freq_ref}{'pop'}) {
				$maxaltAF = ${$freq_ref}{'pop'};
			}
			print OUT "$currchr\t$pos\t$pos\t$ref\t$alt\t".sprintf("%.4f", $errorfreq*100)."\t".sprintf("%.4f", $maxaltAF*100)."\n";	
		}
	}
}
close OUT;








sub readData {
	my $currchr = $_[0];
	
	my %maxaltAFs;
	my $exit_value;
	my $startbp = 1;
	my $endbp = 270000000;
	my ($start1kg, $end1kg) = (1, 270000000);
	my $incrementbp = 30000000;
	
	# note that 1KG Allele Frequency is for the ALTERNATE allele
	for (my $i=0; $i<=270000000; $i+=$incrementbp) {
		$start1kg = $i+1;
		$end1kg = $i+$incrementbp;
		save_parse_1KG($currchr, $start1kg, $end1kg, \%maxaltAFs);
	}
	
	# note that ESP altAF is for MINOR allele
	##INFO=<ID=altAF,Number=.,Type=String,Description="Minor Allele Frequency in percent in the order of EA,AA,All">
	##INFO=<ID=EA_AC,Number=.,Type=String,Description="European American Allele Count in the order of AltAlleles,RefAllele">
	##INFO=<ID=AA_AC,Number=.,Type=String,Description="African American Allele Count in the order of AltAlleles,RefAllele">
	print STDERR "Reading in ESP data for chr $currchr\n";
	# my $espSNPsfile = "$commonvarpath/ESP6500.snps.vcf.gz";	
	# my @espdata = `tabix $espSNPsfile $currchr:$startbp-$endbp`;
	my $espfile = "/nfs/home/jxchong/references/ESP6500SI-V2.chr$currchr.snps_indels.vcf.gz";	
	my @espdata = `tabix $espfile $currchr:$startbp-$endbp`;
	# 1       802311  .       ATCCCTGACG      A       .       PASS    DBSNP=.;EA_AC=191,5105;AA_AC=177,2753;TAC=368,7858;MAF=3.6065,6.041,4.4736;GTS=A1A1,A1R,RR;EA_GTC=88,15,2545;AA_GTC=79,19,1367;GTC=167,34,3912;DP=216;GL=.;CP=0.0;CG=-
	# 2.5;AA=.;CA=.;EXOME_CHIP=no;GWAS_PUBMED=.;GM=.;FG=intergenic;AAC=.;PP=.;CDP=.;GS=.;PH=.;EA_AGE=.;AA_AGE=.
	# 1       802314  .       C       T       .       PASS    DBSNP=.;EA_AC=0,3180;AA_AC=9,1375;TAC=9,4555;MAF=0.0,0.6503,0.1972;GTS=TT,TC,CC;EA_GTC=0,0,1590;AA_GTC=0,9,683;GTC=0,9,2273;DP=216;GL=.;CP=0.0;CG=0.2;AA=C;CA=.;EXOME_CHIP=no;
	# GWAS_PUBMED=.;GM=.;FG=intergenic;AAC=.;PP=.;CDP=.;GS=.;PH=.;EA_AGE=.;AA_AGE=.
	
	$exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@espdata\n";
	}
	foreach (@espdata) {
		$_ =~ s/\s+$//;	
		my @popaltAFs;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		my @altalleles = split(",", $alt);									# in case there are multi-allelic SNPs
		for (my $i=0; $i<=$#altalleles; $i++) {
			{ 
				$vardata[2] =~ /EA_AC=([\d,]+);/;
				my @allelecounts = split(",", $1);
				my $totalleles = sum @allelecounts;
				push(@popaltAFs, $allelecounts[$i]/$totalleles);
			}
			{
				$vardata[2] =~ /AA_AC=([\d,]+);/;
				my @allelecounts = split(",", $1);
				my $totalleles = sum @allelecounts;
				push(@popaltAFs, $allelecounts[$i]/$totalleles);
			}

			my $maxpopaltAF = 0;
			foreach my $varaltAF (@popaltAFs) {
				if ($varaltAF > $maxpopaltAF) {
					$maxpopaltAF = $varaltAF;			}
			}

			my $lookup = "$ref.$altalleles[$i]";
			if (!defined $maxaltAFs{$varpos}{$lookup}{'pop'}) {
				$maxaltAFs{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
				$maxaltAFs{$varpos}{$lookup}{'error'} = 0;
			} elsif ($maxpopaltAF > $maxaltAFs{$varpos}{$lookup}{'pop'}) {
				$maxaltAFs{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
				$maxaltAFs{$varpos}{$lookup}{'error'} = 0;
			}
		}
		# if ($varpos > 1000000) { last; }
	}

	foreach my $file (@errorinputs) {
		print STDERR "Reading in error data for chr $currchr\n";
		my @errordata = `tabix $file $currchr:$startbp-$endbp`;
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
			if (defined $maxaltAFs{$varpos}{$lookup}{'error'}) {
				if ($errorfreq > $maxaltAFs{$varpos}{$lookup}{'error'}) {
					$maxaltAFs{$varpos}{$lookup}{'error'} = $errorfreq;
				}
			} else {
				$maxaltAFs{$varpos}{$lookup}{'error'} = $errorfreq;
				$maxaltAFs{$varpos}{$lookup}{'pop'} = 0;
			}
			# if ($varpos > 1000000) { last; }
		}
	}
	
	return \%maxaltAFs;
}







sub save_parse_1KG {
	my ($currchr, $start1kg, $end1kg, $maxaltAF_ref) = @_;
	my %maxaltAFs = %{$maxaltAF_ref};
	
	my @thousandgenomedata;		
	print STDERR "Reading in 1KG data for chr $currchr with command tabix $thousandgenomesfile $currchr:$start1kg-$end1kg\n";
	@thousandgenomedata = `tabix $thousandgenomesfile $currchr:$start1kg-$end1kg`;
	print STDERR "Finished reading 1KG data for chr $currchr\n";
	my $exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@thousandgenomedata\n";
	}
	print STDERR "Parsing 1KG data for chr $currchr:$start1kg-$end1kg\n";

	foreach (@thousandgenomedata) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		my @populations = qw(ASN AMR AFR EUR);
		my @popaltAFs;
		foreach my $population (@populations) {
			my $altAFfield = $population."_AF";
			$vardata[2] =~ m/;$altAFfield=((\.|\d)+)/;
			if (defined $1) {
				push(@popaltAFs, $1);
			}
		}
		my $maxpopaltAF = 0;
		foreach my $varaltAF (@popaltAFs) {
			if ($varaltAF > $maxpopaltAF) {
				$maxpopaltAF = $varaltAF;
			}
		}

		my $lookup = "$ref.$alt";
		if (!defined $maxaltAFs{$varpos}{$lookup}{'pop'}) {
			$maxaltAFs{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
			$maxaltAFs{$varpos}{$lookup}{'error'} = 0;
		} elsif ($maxpopaltAF > $maxaltAFs{$varpos}{$lookup}{'pop'}) {
			$maxaltAFs{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
			$maxaltAFs{$varpos}{$lookup}{'error'} = 0;
		}
	}
}




sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	# print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	# print "\t--capture\tcapture array used in sequencing (bigexome or v2)\n";
	die;
}