#!/usr/bin/env perl
#
# Description: Add OMIM number and phenotype name, if any
#
#
#
# Created by Jessica Chong on 2014-12-09.


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


my $omimfile = "/nfs/home/jxchong/references/morbidmap.parsed.txt";
my %omimgenes;
if (!-e $omimfile) {
	die "Cannot read from reference data: $omimfile :$?\n";
} else {
	open (MIMDATA, "$omimfile") or die "Cannot read $omimfile: $?\n";
	<MIMDATA>;
	while (<MIMDATA>) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @omimdata = split("\t", $_);
		my @genealiases = split(", ", $omimdata[3]);
		foreach my $gene (@genealiases) {
			if (defined $omimgenes{$gene}{'pheno'}) {
				$omimgenes{$gene}{'pheno'} .= "|$omimdata[0]";
			} else {
				$omimgenes{$gene}{'pheno'} .= "$omimdata[0]";
			}
			$omimgenes{$gene}{'geneMIM'} = $omimdata[4];
		}
	}
	close MIMDATA;
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

open (FILE, "zcat $inputfile |") or die "Cannot read $inputfile file: $!.\n";
if ($inputfiletype eq 'vcf') {						# skip all the metadata lines at the top of the file
	while (<FILE>) {
		$headerline = $_;
		$headerline =~ s/\s+$//;					# Remove line endings
		if ($headerline !~ /^#CHROM/) {
			if ($headerline =~ "##fileformat=VCFv4") {
				print OUT "$headerline\n";
				print OUT '##INFO=<ID=geneMIM,Number=1,Type=String,Description="MIM number for gene entry">'."\n";
				print OUT '##INFO=<ID=phenoOMIM,Number=1,Type=String,Description="Phenotypes associated with gene in OMIM">'."\n";
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
	print OUT join("\t", @header)."\tgeneMIM\tphenoOMIM\n";
}

open (FILE, "zcat $inputfile |") or die "Cannot read $inputfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	next if ($_ =~ /^#/);
	my @line = split ("\t", $_);
	my $infofields_ref = retrieveVCFfields($line[7], 'info');
	my $gene = ${$infofields_ref}{'GL'};
	
	my ($mimpheno, $mimgene) = ('.', '.');
	if (defined $omimgenes{$gene}{'pheno'}) {
		$mimpheno = $omimgenes{$gene}{'pheno'};
	}
	if (defined $omimgenes{$gene}{'geneMIM'}) {
		$mimgene =  $omimgenes{$gene}{'geneMIM'};
	} 
	
	if ($inputfiletype eq 'vcf') {
		print OUT join("\t", @line[0..($vcfinfocol-1)]);
		print OUT "\t$line[$vcfinfocol];phenoOMIM=$mimpheno;geneMIM=$mimgene\t";
		print OUT join("\t", @line[($vcfinfocol+1)..$#line])."\n";
	} else {
		print OUT join("\t", @line)."$mimpheno\t$mimgene\n";
	}
	
}
close FILE;


close OUT;



sub retrieveVCFfields {
	my ($vcfcolumn, $vcfcoltype) = @_;
	my @columncontents;
	my %vcffields;
	foreach my $name (qw(GL FG CP CG PP AAC PP CDP AFCMG AFPOP CA KP PH)) {
		$vcffields{$name} = '.';
	}
	
	if ($vcfcoltype eq 'genotype') {
		@columncontents = split(":", $vcfcolumn);
	} elsif ($vcfcoltype eq 'info') {
		@columncontents = split(";", $vcfcolumn);
	} else {
		@columncontents = split(/[;:]/, $vcfcolumn);
	}
	foreach my $column (@columncontents) {
		my ($fieldname, $value) = split("=", $column);
		if (!defined $value) {
			$value = '.';
		}
		$vcffields{$fieldname} = $value;
	}
	return \%vcffields;
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








sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput SSAnnotation file\n";
	print "\t--out\toutput file\n";
	# print "\t--capture\tcapture array used in sequencing (bigexome or v2)\n";
	die;
}