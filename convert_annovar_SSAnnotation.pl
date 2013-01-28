#!/usr/bin/env perl
#
# Description: Make an SSAnnotation-like file for analysis using my filter_model.pl script.
#
#
#
# Created by Jessica Chong on 2013-01-25.


use strict;
use warnings;
use Getopt::Long;


my ($annovarinput, $vcfinput, $outputfile);

GetOptions(
	'annovar=s' => \$annovarinput, 
	'vcf=s' => \$vcfinput, 
	'out=s' => \$outputfile,
);

if (!defined $vcfinput) {
	optionUsage("option --vcf not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} 

my %iupac = (
	'A/G' => 'R',
	'C/T' => 'Y',
	'G/C' => 'S',
	'A/T' => 'W',
	'G/T' => 'K',
	'A/C' => 'M',
	'G/A' => 'R',
	'T/C' => 'Y',
	'C/G' => 'S',
	'T/A' => 'W',
	'T/G' => 'K',
	'C/A' => 'M',
	'C/G/T' => 'B',
	'A/G/T' => 'D',
	'A/C/T' => 'H',
	'A/C/G' => 'V',
);


# perl ~/bin/annovar/convert2annovar.pl .vcf -format vcf4 > output.annovar

# my $gencodeannotate = `perl ~/bin/annovar/annotate_variation.pl $annovarinput ~/bin/annovar/humandb/ -dbtype wgEncodeGencodeBasicV14 --buildver hg19`;
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype avsift park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19 
# perl ~/bin/annovar/annotate_variation.pl park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ -filter -dbtype ljb_pp2  --buildver hg19 
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype gerp++gt2 park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype esp6500_all park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype 1000g2012apr_al park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19



my %annovar;
my $linecounter = 1;
open (FILE, "$annovarinput.variant_function") or die "Cannot read $annovarinput.variant_function file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	$annovar{$linecounter} = [@line[0..1]];		# exonic or not, gene
	$linecounter++;
}
close FILE;


open (FILE, "$annovarinput.exonic_variant_function") or die "Cannot read $$annovarinput.exonic_variant_function file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	my $linenum = $line[0];
	$linenum =~ s/line//;
	push(@{$annovar{$linenum}}, @line[1..2]);	# variant function, gene-related info
}
close FILE;


# my @annovarannot = @{$annovar{$linecounter}};
# my ($newgene, $newfunc, $newaa, $newproteinpos, $newcdnapos, $ensembl, $exon);
# $newgene = $annovarannot[1];
# 
# if ($annovarannot[0] =~ 'exonic' && $annovarannot[0] !~ 'ncRNA') {
# 	# if ($annovarannot[1] !~ 'synonymous') {
# 	# 	print "@annovarannot\n";
# 	# }
# 	
# 	if ($annovarannot[2] !~ m/unknown/i) {
# 		$newfunc = translateGencodefxn($annovarannot[2]);
# 		if ($annovarannot[0] =~ 'splic') {
# 			$newfunc .= "splicing";
# 		}
# 		if ($annovarannot[3] =~ 'wholegene') {
# 			($newaa, $newproteinpos, $newcdnapos) = qw(wholegene wholegene wholegene);
# 		} else {
# 			my @varannotations = my @annotation = split(",", $annovarannot[3]);
# 			($newgene, $ensembl, $exon, $newcdnapos, $newproteinpos) = split(":", $varannotations[0]);
# 			$newaa = $newproteinpos;
# 		}
# 	} else {
# 		($newfunc, $newaa, $newproteinpos, $newcdnapos) = qw(unkn unkn unkn unkn);
# 	}
# } else {
# 	($newfunc, $newaa, $newproteinpos, $newcdnapos) = qw(NA NA NA NA);
# }
# 
# 
# my %gene;
# open (GENE, "$annovarinput.variant_function");
# 	
# close GENE;



my %vcfdata;
my @subjectgenotypecols;
my $beginsubjectcols;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "$vcfinput") or die "Cannot read $vcfinput file: $!.\n";
while ( <FILE> ) {
	if ($_ =~ '##') { next; }
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	if ($_ =~ '#CHR') {
		my @header = split("\t", $_);
		for (my $i=0; $i<=$#header; $i++) {
			my $columnname = $header[$i];
			if ($columnname eq 'FORMAT') {
				$beginsubjectcols = $i+1;
				last;
			}
		}
		@subjectgenotypecols = $beginsubjectcols..$#header;
		print "Chr\tPos\tref\talt\tfilterFlagGATK\tgeneList\tfunctionGVS\tproteinchange\tcDNAchange\tconsScoreGERP\tpolyPhen2\tSIFT\t";
		print join("\t", @header[@subjectgenotypecols])."\t";
		print join("Depth\t", @header[@subjectgenotypecols])."Depth\t";
		print join("Qual\t", @header[@subjectgenotypecols])."Qual\t";
		print "PrctFreqinCMG\tPrctFreqinOutsidePop\n";
	} else {
		my ($chr, $pos, $rsid, $ref, $alt, $varqual, $varfilter, $varinfo, $varformat) = @line[0..($beginsubjectcols-1)];
		my (@subjectgts, @subjectgtquals, @subjectdepths);
		if ($varformat eq 'GT:AD:DP:GQ:PL') {
			foreach my $datastring (@line[@subjectgenotypecols]) {
				if ($datastring eq './.') {
					push(@subjectgts, 'N');
					push(@subjectgtquals, "NA");
					push(@subjectdepths, "NA");
				} else {
					my @metrics = split(":", $datastring);
					push(@subjectgts, vcfgeno2iupac($metrics[0], $ref, $alt));
					push(@subjectgtquals, $metrics[3]);
					push(@subjectdepths, $metrics[2]);	
				}
			}
		}
		my ($gene, $function, $protein, $cdna, $gerp, $polyphen, $sift);
		
		
		
		# $vcfdata{$chr}{$pos}{"$ref.$alt"} = [$varfilter, join("\t", @subjectgts), join("\t", @subjectgtquals), join("\t", @subjectdepths)];	
		if (length($ref)>12) {
			$ref = 0;
		}
		if (length($alt)>12) {
			$alt = 1;
		}
		print "$chr\t$pos\t$ref\t$alt\t$varfilter\tGene"."\t".join("\t", @subjectgts)."\t".join("\t", @subjectdepths)."\t".join("\t", @subjectgtquals)."\n";
		exit;
	}
}
close FILE;
close OUT;


sub vcfgeno2iupac {
	my ($vcfgeno, $ref, $alt) = @_;
	my $newgeno = $vcfgeno;
	if (length($ref)==1 && length($alt)==1) {
		if ($vcfgeno eq '0/1') {
			$newgeno = $iupac{$vcfgeno};
		} elsif ($vcfgeno eq '0/0') {
			$newgeno = $ref;
		} elsif ($vcfgeno eq '1/1') {
			$newgeno = $alt;
		}
	} elsif (length($ref)<=12 && length($alt)<=12) {
		if ($vcfgeno eq '0/1') {
			$newgeno = "$ref/$alt";
		} elsif ($vcfgeno eq '0/0') {
			$newgeno = "$ref/$ref";
		} elsif ($vcfgeno eq '1/1') {
			$newgeno = "$alt/$alt";
		}
	}
	return $newgeno;
}








sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	die;
}