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

# perl ~/bin/annovar/convert2annovar.pl .vcf -format vcf4 > output.annovar

# my $gencodeannotate = `perl ~/bin/annovar/annotate_variation.pl $annovarinput ~/bin/annovar/humandb/ -dbtype wgEncodeGencodeBasicV14 --buildver hg19`;
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype avsift park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19 
# perl ~/bin/annovar/annotate_variation.pl park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ -filter -dbtype ljb_pp2  --buildver hg19 
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype gerp++gt2 park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype esp6500_all park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19
# perl ~/bin/annovar/annotate_variation.pl -filter -dbtype 1000g2012apr_al park.remapped_and_pindel.sorted.annovar ~/bin/annovar/humandb/ --buildver hg19



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
		print "Chr\tPos\tref\talt\tfilterFlagGATK\tgeneList\t".join("\t", @header[@subjectgenotypecols])."\t";
		print "functionGVS\tproteinchange\tcDNAchange\tconsScoreGERP\tpolyPhen2\tSIFT\t";
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
					push(@subjectgts, $metrics[0]);
					push(@subjectgtquals, $metrics[3]);
					push(@subjectdepths, $metrics[2]);	
				}
			}
		}
		# $vcfdata{$chr}{$pos}{"$ref.$alt"} = [$varfilter, join("\t", @subjectgts), join("\t", @subjectgtquals), join("\t", @subjectdepths)];	
		# print "$_\n";
		print "$chr\t$pos\t$ref\t$alt\t$varfilter\tGene"."\t".join("\t", @subjectgts)."\t".join("\t", @subjectdepths)."\t".join("\t", @subjectgtquals)."\n";
		exit;
	}
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