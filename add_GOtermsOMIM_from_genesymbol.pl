#!/usr/bin/env perl
#
# Description: Process ENSEMBL Biomart output (GO Term Name, MIM Morbid Description, Description, Associated Gene Name) to make one line per gene
#
# Usage: perl untitled
#
#
# Created by Jessica Chong on 2012-11-16.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $martfile);

GetOptions(
	'in=s' => \$inputfile, 
	'mart=s' => \$martfile,
	'out=s' => \$outputfile,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $martfile) {
	optionUsage("option --mart not defined\n");
}


my %martannotations;
my $currentgene = 'NA';
my $currentdescription;
my %goannotations;
my %mimdescription;
open (MART, $martfile) or die "Cannot read $martfile: $!.\n";
while (<MART>) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($goterm, $omim, $description, $genesymbol) = split ("\t", $_);
	$description =~ s/ \[Source:HGNC Symbol;Acc:\d+\]//;
	if ($genesymbol ne $currentgene) {
		# output old gene's info
		if ($currentgene ne 'NA') {
			$martannotations{$currentgene} = "$currentgene\t$currentdescription\t".join(";", keys %goannotations)."\t".join(";", keys %mimdescription);			
		}
		# starting a new gene
		%goannotations = ();
		%mimdescription = ();
		$goannotations{$goterm} = 1;
		$mimdescription{$omim} = 1;
		$currentdescription = $description;
		$currentgene = $genesymbol;
	} else {
		if ($goterm ne '') {
			$goannotations{$goterm} = 1;
		}
		if ($omim ne '') {
			$mimdescription{$omim} = 1;
		}
	}
}
close MART;



open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
my $header = <FILE>;
$header =~ s/\s+$//;
print OUT "$header\tSymbol\tDescription\tGO Terms\tOMIM\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split("\t", $_);
	my $thisgene = $line[7];
	my $annotation;
	if (!exists $martannotations{$thisgene}) {
		$annotation = "\t\t\t";
	} else {
		$annotation = $martannotations{$thisgene};
	}
	print OUT join("\t", @line)."\t$annotation\n";
}
close FILE;
close OUT;





sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	die;
}