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
my $thousandgenomesfile = "$commonvarpath/1KG_phase3_download-2015-02-16/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz";
my $exacpath = "$commonvarpath/broad_exome_63k/ExAC.r0.1.sites.vep.vcf.gz";
my $errorpath = "$commonvarpath/systematic_error/2013_nov/systematic_errors.vcf.gz";

# temporary files
# my $commonvarpath = '/nfs/home/jxchong/jessica_annovar';
# my $thousandgenomesfile = "$commonvarpath/1kg.chrXY.vcf.gz";
# my $espSNPsfile = "$commonvarpath/ESP.chrXYsnps.vcf.gz";


open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
print OUT "#chr\tstart\tend\tref\talt\tPrctAltFreqinCMG\tPrctAltFreqinOutsidePop\n";
foreach my $currchr (((1..22), "X", "Y", "M")) {
# foreach my $currchr ((22)) {
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
			my $altAF_sourceDB = "";
			if (defined ${$freq_ref}{'pop'}) {
				$maxaltAF = ${$freq_ref}{'pop'};
				$altAF_sourceDB = ${$freq_ref}{'DB'};
			}
			print OUT "$currchr\t$pos\t$pos\t$ref\t$alt\t".sprintf("%.4f", $errorfreq*100)."\t".sprintf("%.4f", $maxaltAF*100)."\t$altAF_sourceDB\n";	
		}
	}
}
close OUT;








sub readData {
	my $currchr = $_[0];
	
	my %maxaltAFs;
	my $maxaltAFs_ref = \%maxaltAFs;
	my $startbp = 1;
	my $endbp = 270000000;
	my ($startsubset, $endsubset) = (1, 270000000);
	my $incrementbp = 30000000;
	
	# ($startsubset, $endsubset) = (16287372, 16287372);
	# ($startsubset, $endsubset) = (1, 100000);
	
	for (my $startsubset=$startbp; $startsubset<=$endbp; $startsubset+=$incrementbp) {
		$endsubset = $startsubset+$incrementbp;
		print STDOUT "Reading 1KG data: tabix $errorpath $currchr:$startsubset-$endsubset\n";
		save_parse_1KG($currchr, $startsubset, $endsubset, $maxaltAFs_ref);
		print STDOUT "Reading ExAC data: tabix $errorpath $currchr:$startsubset-$endsubset\n";
		save_parse_exac($currchr, $startsubset, $endsubset, $maxaltAFs_ref);
		print STDOUT "Reading ESP data: tabix $errorpath $currchr:$startsubset-$endsubset\n";
		save_parse_esp($currchr, $startsubset, $endsubset, $maxaltAFs_ref);
		print STDOUT "Reading systematic_errors data: tabix $errorpath $currchr:$startsubset-$endsubset\n";
		save_parse_errors($currchr, $startsubset, $endsubset, $maxaltAFs_ref);
		print STDOUT "\n";
	}
	
	# note that 1KG Allele Frequency is for the ALTERNATE allele
	# note that ESP altAF is for MINOR allele
	
	return $maxaltAFs_ref;
}


sub save_parse_esp {
	my ($currchr, $startbp, $endbp, $maxaltAFs_ref) = @_;
	my $DBname = 'ESP';
	
	##INFO=<ID=EA_AC,Number=.,Type=String,Description="European American Allele Count in the order of AltAlleles,RefAllele">
	##INFO=<ID=AA_AC,Number=.,Type=String,Description="African American Allele Count in the order of AltAlleles,RefAllele">
	my $espfile = "/nfs/home/jxchong/references/ESP6500SI-V2.chr$currchr.snps_indels.vcf.gz";	
	my @espdata = `tabix $espfile $currchr:$startbp-$endbp`;
	my $exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@espdata\n";
	}
	print STDOUT "Parsing $DBname data for chr $currchr:$startbp-$endbp\n";
	foreach (@espdata) {
		$_ =~ s/\s+$//;	
		my @popaltAFs;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		my @altalleles = split(",", $alt);								
		my %popaltAFs;
		@popaltAFs{@altalleles} = (0) x scalar(@altalleles);
		my @populations = qw(EA AA);
		foreach my $population (@populations) {
			my $ACfield = $population."_AC";
			my ($AC) = $vardata[2] =~ m/;$ACfield=([\d,\.]+)/;
			my @allelecounts = split(",", $AC);
			my $totchrcount = sum @allelecounts;
			my @altACvals = @allelecounts[0..($#allelecounts-1)];
			determine_maxpopaltAF($chr, $varpos, \@altalleles, \@altACvals, $totchrcount, \%popaltAFs, $population, $DBname);
		}
		store_maxpopaltAF($ref, $varpos, \@altalleles, \%popaltAFs, $maxaltAFs_ref, $DBname);
	}
}


sub save_parse_exac {
	my ($currchr, $startbp, $endbp, $maxaltAFs_ref) = @_;
	my $DBname = 'ExAC';
	
	my @exacdata = `tabix $exacpath $currchr:$startbp-$endbp`;
	my $exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@exacdata\n";
	}
	print STDOUT "Parsing $DBname data for chr $currchr:$startbp-$endbp\n";

	foreach (@exacdata) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		my @altalleles = split(",", $alt);
		my %popaltAFs;
		@popaltAFs{@altalleles} = (0) x scalar(@altalleles);
		my @populations = qw(AFR AMR EAS NFE SAS);
		foreach my $population (@populations) {
			my $altACfield = "AC_".$population;
			my $chrcountfield = "AN_".$population;
			my ($altAC) = $vardata[2] =~ m/;$altACfield=([\d,.]+)/;
			my ($totchrcount) = $vardata[2] =~ m/;$chrcountfield=(\d+)/;
			my @altACvals = split(",", $altAC);
			determine_maxpopaltAF($chr, $varpos, \@altalleles, \@altACvals, $totchrcount, \%popaltAFs, $population, $DBname);
		}
		store_maxpopaltAF($ref, $varpos, \@altalleles, \%popaltAFs, $maxaltAFs_ref, $DBname);
	}
}

sub save_parse_1KG {
	my ($currchr, $start1kg, $end1kg, $maxaltAFs_ref) = @_;
	my $DBname = '1KG';
	
	my @thousandgenomedata = `tabix $thousandgenomesfile $currchr:$start1kg-$end1kg`;
	my $exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@thousandgenomedata\n";
	}
	print STDOUT "Parsing $DBname data for chr $currchr:$start1kg-$end1kg\n";

	foreach (@thousandgenomedata) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;
		my ($chr, $varpos, $rsid, $ref, $alt, @vardata) = split("\t", $_);
		# print STDERR "Checking $chr:$varpos-$varpos alt alleles = $alt\n";
		my @altalleles = split(",", $alt);
		my %popaltAFs;
		@popaltAFs{@altalleles} = (0) x scalar(@altalleles);
		my @populations = qw(EAS SAS AFR EUR AMR);
			
		foreach my $population (@populations) {
			my $altAFfield = $population."_AF";
			my ($altAFfieldval) = $vardata[2] =~ m/;$altAFfield=([\d,.]+)/;
			if (!defined $altAFfieldval) {
				next;
			}
			my @altAFvals = split(",", $altAFfieldval);
			
			# 1000 genomes already has AF calculated
			for (my $altallele_idx=0; $altallele_idx<=$#altalleles; $altallele_idx++) {
				my $altAF = $altAFvals[$altallele_idx];
				if ($altAF > $popaltAFs{$altalleles[$altallele_idx]}) {
					$popaltAFs{$altalleles[$altallele_idx]} = $altAF;
					# print STDERR "$DBname: AF for $chr:$varpos allele $altalleles[$altallele_idx] (#$altallele_idx) in $population = $altAF\n";
				}
			}	
		}
		store_maxpopaltAF($ref, $varpos, \@altalleles, \%popaltAFs, $maxaltAFs_ref, $DBname);
	}
}

sub save_parse_errors {
	my ($currchr, $startbp, $endbp, $maxaltAFs_ref) = @_;
	my $DBname = 'systematic_errors';
	
	my @errordata = `tabix $errorpath $currchr:$startbp-$endbp`;
	my $exit_value = $? >> 8;
	if ($exit_value != 0) {
		die "Error: exit value $exit_value\n@errordata\n";
	}
	print STDOUT "Parsing $DBname data for chr $currchr:$startbp-$endbp\n";
	
	foreach (@errordata) {
		if ($_ =~ '#') { next; }
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $varpos, $type, $ref, $alt, @vardata) = split("\t", $_);
		# alternate alleles are not a problem because Martin just makes them into a new entry
		my @errorfreqinfo = split("/", $vardata[2]);
		$errorfreqinfo[0] =~ s/OBSERVED=//;
		my $errorfreq = 0;
		if ($errorfreqinfo[1] != 0) {
			$errorfreq = $errorfreqinfo[0]/$errorfreqinfo[1];
		}

		my $lookup = "$ref.$alt";
		if (!defined $maxaltAFs_ref->{$varpos}{$lookup}{'error'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'error'} = $errorfreq;
		} elsif ($errorfreq > $maxaltAFs_ref->{$varpos}{$lookup}{'error'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'error'} = $errorfreq;
		}
		if (!defined $maxaltAFs_ref->{$varpos}{$lookup}{'pop'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'pop'} = 0;
		}
	}
}



sub determine_maxpopaltAF {
	my ($chr, $varpos, $altalleles_ref, $altACvals_ref, $totchrcount, $popaltAFs_ref, $population, $DBname) = @_;
	
	my @altalleles = @{$altalleles_ref};
	# print STDERR "$DBname: $chr:$varpos the alt AC for $population = @{$altACvals_ref}\n";
	for (my $altallele_idx=0; $altallele_idx<=$#altalleles; $altallele_idx++) {
		# print STDERR "$DBname: allele $altalleles[$altallele_idx] (#$altallele_idx), altAF value will = ${$altACvals_ref}[$altallele_idx]/$totchrcount\n";
		my $altAF;
		if ($totchrcount == 0) {
			$altAF = 0
		} else {
			$altAF = ${$altACvals_ref}[$altallele_idx]/$totchrcount;
		}
		if ($altAF > $popaltAFs_ref->{$altalleles[$altallele_idx]}) {
			$popaltAFs_ref->{$altalleles[$altallele_idx]} = $altAF;
			# print STDERR "$DBname: AF for $chr:$varpos allele $altalleles[$altallele_idx] (#$altallele_idx) in $population = $altAF\n";
		}
	}	
}

sub store_maxpopaltAF {
	my ($ref, $varpos, $altalleles_ref, $popaltAFs_ref, $maxaltAFs_ref, $DBname) = @_;
	
	my @altalleles = @{$altalleles_ref};
	for (my $altallele_idx=0; $altallele_idx<=$#altalleles; $altallele_idx++) {
		my $lookup = "$ref.$altalleles[$altallele_idx]";
		my $maxpopaltAF = $popaltAFs_ref->{$altalleles[$altallele_idx]};
		if (!defined $maxaltAFs_ref->{$varpos}{$lookup}{'pop'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
			$maxaltAFs_ref->{$varpos}{$lookup}{'DB'} = $DBname;
		} elsif ($maxpopaltAF > $maxaltAFs_ref->{$varpos}{$lookup}{'pop'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'pop'} = $maxpopaltAF;
			$maxaltAFs_ref->{$varpos}{$lookup}{'DB'} = $DBname;
		}
		if (!defined $maxaltAFs_ref->{$varpos}{$lookup}{'error'}) {
			$maxaltAFs_ref->{$varpos}{$lookup}{'error'} = 0;
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