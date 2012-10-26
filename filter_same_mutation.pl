#!/usr/bin/env perl
#
# Description: Check for heterozygosity for same mutation
#
# Usage: perl untitled
#
#
# Created by Jessica Chong on 2012-10-15.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $subjectdeffile);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'subjectreq=s' => \$subjectdeffile,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $subjectdeffile) {
	optionUsage("option --subjectreq not defined\n");
}

my @orderedsubjects;
my %subjects;
open (SUBJECTS, "$subjectdeffile");
while (<SUBJECTS>) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($subjectid, $desiredgeno) = split("\t", $_);
	$subjects{$subjectid} = $desiredgeno;
	push(@orderedsubjects, $subjectid);
}
close SUBJECTS;



my $countoutputvariants = 0;
my $countinputvariants = 0;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
my @subjectcolumns;
for (my $i=0; $i<=$#header; $i++) {
	if (defined $subjects{$header[$i]}) {
		push(@subjectcolumns, $i);
	}
}
print OUT "$headerline\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	if ($_ =~ '#') {
		next;
	}
	$countinputvariants++;
	my @line = split ("\t", $_);
	my $countmatches = 0;
	my ($vartype, $ref, $alt);
	
	if ($inputfile =~ m/SeattleSeqAnnotation134/) {
		if ($inputfile =~ m/snps/) {
			my $sampleallelestring = $line[5];
			my @samplealleles = split("/", $sampleallelestring);
			if (!defined $samplealleles[1]) {
				$samplealleles[1] = 'N';
				next;
			}
			$vartype = 'SNP';
			($ref, $alt) = selectRefAlt($line[3], $samplealleles[0], $samplealleles[1]);
			my @subjectgenotypes = split(",", $line[4]);
			for (my $i=0; $i<=$#subjectgenotypes; $i++) {
				my $genotype = $subjectgenotypes[$i];
				my $subjectid = $orderedsubjects[$i];
				my $desiredgeno = $subjects{$subjectid};
				$countmatches += checkGenoMatch($vartype, $ref, $alt, $genotype, $desiredgeno);
				# print "$line[2]\t$vartype\tref=$ref\talt=$alt\tgeno=$genotype\t$subjectid\tdesired=$desiredgeno\t$countmatches\n";
			}
			# print "\n";
		}
		if ($inputfile =~ m/indels/) {
			# ($vartype, $ref, $alt) = ('indel', $line[3], $line[4]);		# complicated parsing; not implemented
		}
	} elsif ($inputfile =~ m/SSAnnotation/) {
		($vartype, $ref, $alt) = ($line[2], $line[3], $line[4]);
		for (my $i=0; $i<=$#subjectcolumns; $i++) {
			my $columnidx = $subjectcolumns[$i];
			my $genotype = $line[$columnidx];
			my $subjectid = $header[$columnidx];
			my $desiredgeno = $subjects{$subjectid};
			$countmatches += checkGenoMatch($vartype, $ref, $alt, $genotype, $desiredgeno);
		}
	}
	
	# if ($line[2] == 17538) {
	# 	exit;																		# DEBUG
	# }
	
	if ($countmatches == scalar(keys %subjects)) {
		print OUT join("\t", @line)."\n";
		$countoutputvariants++;
	}
}
close FILE;
close OUT;

print "Matched $countoutputvariants out of $countinputvariants\n";






sub checkGenoMatch {
	my ($vartype, $ref, $alt, $genotype, $desiredgeno) = @_;
		
	my @alleles;
	my $ismatch = 0;
	
	if ($genotype eq 'N') {
		$ismatch = 1;
	} else {
		if ($vartype eq 'SNP') {
			if ($genotype eq $ref) {
				@alleles = ($ref, $ref);
			} elsif ($genotype eq $alt) {
				@alleles = ($alt, $alt);
			} else {
				@alleles = ($ref, $alt);
			}
		} elsif ($vartype eq 'indel') {
			@alleles = split("/", $genotype);
		}
	
		# print "$vartype\t$genotype\t$desiredgeno\t$ref\t$alt\t@alleles\n";				# DEBUG
	
		if ($desiredgeno eq 'ref') {
			if ($alleles[0] eq $ref && $alleles[1] eq $ref) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'alt') {
			if ($alleles[0] eq $alt && $alleles[1] eq $alt) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'het') {
			if ($alleles[0] eq $ref && $alleles[1] eq $alt) {
				$ismatch = 1;
			}
		}
	}
	return $ismatch;
}


sub selectRefAlt {
	my ($origref, $allele1, $allele2) = @_;
	my ($ref, $alt);
	if ($origref eq $allele1) {
		$ref = $allele1;
		$alt = $allele2;
	} else {
		$ref = $allele2;
		$alt = $allele1;
	}
	return $ref, $alt;
}




sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	die;
}