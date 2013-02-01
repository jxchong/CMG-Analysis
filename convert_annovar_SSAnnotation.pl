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


# perl ~/bin/annovar/convert2annovar.pl park.remapped_and_pindel.sorted.reminvalid.vcf -format vcf4 > park.remapped_and_pindel.sorted.reminvalid.annovar
# perl ~/bin/annovar/annotate_variation.pl park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ -dbtype wgEncodeGencodeBasicV14 -geneanno --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype avsift park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ --buildver hg19 
# perl ~/bin/annovar/annotate_variation.pl park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ -filter -dbtype ljb_pp2  --buildver hg19 

# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype gerp++gt2 park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype esp6500_all park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype 1000g2012apr_al park.remapped_and_pindel.sorted.reminvalid.annovar ~/bin/annovar/humandb/ --buildver hg19



my %annovar;
my $linecounter = 1;
open (FILE, "$annovarinput.variant_function") or die "Cannot read $annovarinput.variant_function file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($component, $gene, $chr, $start, $stop, @line) = split ("\t", $_);
	if (defined $annovar{"$chr.$start"}) {
		print "$chr.$start already exists: $_\n";
	}
	$annovar{"$chr.$start.$stop."}{'variantcomponent'} = $component;					# exonic or not, gene
	$annovar{"$chr.$start.$stop"}{'gene'} = $gene;									# gene
	$linecounter++;
}
close FILE;


#  whole gene (frameshift insertion, annovar=+1 vs vcf)

open (FILE, "$annovarinput.exonic_variant_function") or die "Cannot read $annovarinput.exonic_variant_function file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($linenum, $mutationtype, $geneannot, $chr, $start, $stop, @line) = split ("\t", $_);	
	my ($newgene, $newfunc, $newaa, $newproteinpos, $newcdnapos, $ensembl, $exon);	
	if ($annovar{"$chr.$start.$stop"}{'variantcomponent'} =~ 'exonic' && $annovar{"$chr.$start"}{'variantcomponent'} !~ 'ncRNA') {	
		if ($mutationtype !~ m/unknown/i) {
			$newfunc = translateGencodefxn($mutationtype);
			if ($annovar{"$chr.$start.$stop"}{'variantcomponent'} =~ 'splic') {
				$newfunc .= "splicing";
			}
			if ($geneannot =~ 'wholegene') {
				($newaa, $newproteinpos, $newcdnapos) = qw(wholegene wholegene wholegene);
			} else {
				my @varannotations = split(",", $geneannot);
				($newgene, $ensembl, $exon, $newcdnapos, $newproteinpos) = split(":", $varannotations[0]);
				$newaa = $newproteinpos;
			}
		} else {
			($newfunc, $newaa, $newproteinpos, $newcdnapos) = qw(unkn unkn unkn unkn);
		}
	} else {
		($newfunc, $newaa, $newproteinpos, $newcdnapos) = qw(NA NA NA NA);
	}
	
	$annovar{"$chr.$start.$stop"}{'mutationfunction'} = $newfunc;
	$annovar{"$chr.$start.$stop"}{'mutationlocusinfo'} = [$newaa, $newcdnapos];
}
close FILE;



my @annotationsources = qw(gerp++gt2 avsift ljb_pph2);
foreach my $annotsource (@annotationsources) {
	open (FILE, "$annovarinput.hg19_".$annotsource."_dropped") or die "Cannot read $annovarinput.hg19_".$annotsource."_dropped file: $!.\n";
	while ( <FILE> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($annotsource, $thisannot, $chr, $start, $stop, @line) = split ("\t", $_);
		$annovar{"$chr.$start.$stop"}{$annotsource} = $thisannot;		
	}
	close FILE;
}


my %vcfdata;
my @subjectgenotypecols;
my $beginsubjectcols;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "zcat $vcfinput |") or die "Cannot read $vcfinput file: $!.\n";
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
		print join("Qual\t", @header[@subjectgenotypecols])."Qual\n";
		# print "PrctFreqinCMG\tPrctFreqinOutsidePop\n";
	} else {
		my ($chr, $pos, $rsid, $ref, $alt, $varqual, $varfilter, $varinfo, $varformat) = @line[0..($beginsubjectcols-1)];
		my (@subjectgts, @subjectgtquals, @subjectdepths);
		if ($varformat eq 'GT:AD:DP:GQ:PL') {
			foreach my $datastring (@line[@subjectgenotypecols]) {
				if ($datastring eq './.') {
					push(@subjectgts, 'N');
					push(@subjectgtquals, "0");
					push(@subjectdepths, "0");
				} else {
					my @metrics = split(":", $datastring);
					push(@subjectgts, vcfgeno2iupac($metrics[0], $ref, $alt));
					push(@subjectgtquals, $metrics[3]);
					push(@subjectdepths, $metrics[2]);	
				}
			}
		}
		
		my $variantlookup = "$chr.$pos.$pos";
		my $gene = $annovar{$variantlookup}{'gene'};
		my $component = $annovar{$variantlookup}{'variantcomponent'};
		my $function = $annovar{$variantlookup}{'mutationfunction'};
		my ($protein, $cdna) = @{$annovar{$variantlookup}{'mutationlocusinfo'}};
		my $gerp = $annovar{$variantlookup}{'gerp++gt2'};
		my $polyphen = $annovar{$variantlookup}{'ljb_pph2'};
		my $sift = $annovar{$variantlookup}{'avsift'};

		if (length($ref)>12) {
			$ref = 0;
		}
		if (length($alt)>12) {
			$alt = 1;
		}
		print "$chr\t$pos\t$ref\t$alt\t$varfilter\t$gene\t$function\t$protein\t$cdna\t$gerp\t$polyphen\t$sift\t";
		print join("\t", @subjectgts)."\t".join("\t", @subjectdepths)."\t".join("\t", @subjectgtquals)."\n";
		# print "\t$prctfrqcmg\t$prctfrqoutside\n";
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