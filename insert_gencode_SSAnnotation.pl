#!/usr/bin/env perl
#
# Description:
#
# Usage: perl untitled
#
#
# Created by Jessica Chong on 2013-01-24.


use strict;
use warnings;
use Getopt::Long;


my ($ssannotationfile, $variantfunctionfile, $exonicfunctionfile, $outputfile);

GetOptions(
	'SSAnnotation=s' => \$ssannotationfile, 
	'variant_function=s' => \$variantfunctionfile, 
	'exonic_function=s' => \$exonicfunctionfile, 
	'out=s' => \$outputfile,
);

if (!defined $ssannotationfile) {
	optionUsage("option --SSAnnotation not defined\n");
} elsif (!defined $variantfunctionfile) {
	optionUsage("option --variant_function not defined\n");
} 

my %annovar;
my $linecounter = 1;
open (FILE, "$variantfunctionfile") or die "Cannot read $variantfunctionfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	$annovar{$linecounter} = [@line[0..1]];		# exonic or not, gene
	$linecounter++;
}
close FILE;


open (FILE, "$exonicfunctionfile") or die "Cannot read $exonicfunctionfile file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my $linenum = $line[0];
	$linenum =~ s/line//;
	push(@{$annovar{$linenum}}, @line[1..2]);	# variant function, gene-related info
}
close FILE;



open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
$linecounter = 1;
open (FILE, "$ssannotationfile") or die "Cannot read $ssannotationfile file: $!.\n";
my $header = <FILE>;
print OUT "$header";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my @annovarannot = @{$annovar{$linecounter}};
	my ($newgene, $newfunc, $newaa, $newproteinpos, $newcdnapos, $ensembl, $exon);
	$newfunc = $annovarannot[0];
	$newgene = $annovarannot[1];
	
	if ($annovarannot[0] =~ 'exonic' && $annovarannot[0] !~ 'ncRNA') {
		# if ($annovarannot[1] !~ 'synonymous') {
		# 	print "@annovarannot\n";
		# }
		
		if ($annovarannot[2] !~ m/unknown/i) {
			$newfunc = translateGencodefxn($annovarannot[2]);
			if ($annovarannot[0] =~ 'splic') {
				$newfunc .= "splicing";
			}
			if ($annovarannot[3] =~ 'wholegene') {
				($newaa, $newproteinpos, $newcdnapos) = qw(wholegene wholegene wholegene);
			} else {
				my @varannotations = my @annotation = split(",", $annovarannot[3]);
				($newgene, $ensembl, $exon, $newcdnapos, $newproteinpos) = split(":", $varannotations[0]);
				$newaa = $newproteinpos;
			}
		} else {
			($newfunc, $newaa, $newproteinpos, $newcdnapos) = qw(unkn unkn unkn unkn);
		}
	} else {
		($newaa, $newproteinpos, $newcdnapos) = qw(NA NA NA);
	}
	# replace genelist, functionGVS, aminoAcids, proteinPosition, cDNAposition
	print OUT join("\t", @line[0..6]);
	print OUT "\t$newgene\t".join("\t", @line[8..9]);
	print OUT "\t$newfunc\t$newaa\t$newproteinpos\t$newcdnapos\t".join("\t", @line[14..$#line])."\n";
	$linecounter++;
}
close FILE;

close OUT;




sub translateGencodefxn {
	my $gencodefxn = $_[0];
	$gencodefxn =~ s/ SNV//;
	
	if ($gencodefxn eq 'nonsynonymous') {
		$gencodefxn = 'missense';
	}
	if ($gencodefxn eq 'stopgain') {
		$gencodefxn = 'stop-gained';
	}
	if ($gencodefxn eq 'stoploss') {
		$gencodefxn = 'stop-lost';
	}
	if ($gencodefxn eq 'UTR3') {
		$gencodefxn = 'utr-3';
	}
	if ($gencodefxn eq 'UTR5') {
		$gencodefxn = 'utr-5';
	}
	if ($gencodefxn eq 'intronic') {
		$gencodefxn = 'intron';
	}
	return $gencodefxn;
}



sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	die;
}