#!/usr/bin/env perl
#
# Description: Combine annovar annotations in tab-delimited file.  Will load all annotations in memory so best for smallish numbers of variants
#	Keep vcf input allele style for matching later on by using:
#		convert2annovar.pl -format vcf4 --includeinfo --allallele input.vcf | cut -f1-10 | sort -k1,1 -k2,2n > output.annovar
#
#
# Created by Jessica Chong on 2013-04-03.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($outputfile, $annovarstem, $annovarfilter, $help);

GetOptions(
	'annovarstem=s' => \$annovarstem,
	'annovarfilters=s' => \$annovarfilter,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} elsif (!defined $annovarstem) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --annovarstem not defined\n");
} elsif (!defined $annovarfilter) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --annovarfilters not defined\n");
}




# <li>The VCF definition would allow for Mutation entries like &lt;DEL&gt; or comma separated lists of observed bases 
#  * (such as TA,TC,TG). Annovar throws out these entries, so there won't be any Annovar annotations for them.</li>
#  * <li>Annovar writes the DIPs in a different format. Here a comparison of how the 2 tools would write the mutations:
#  * <ul><li>Insertion: GATK:    Ref: T Obs: TACG</li>
#  * <li>               Annovar: Ref: - Obs:  ACG</li>
#  * <li>Deletion:  GATK:    Ref: TACG Obs: T</li>
#  * <li>           Annovar: Ref:  ACG Obs: -</li>
#  * </ul></li>
#  * <li>The start position of DIPs that are Deletions in annovar is +1 of the start position GATK would write. 
#  * So for deletions we have to add 1 to the start position when looking for the annovar entry.</li>
# http://www.icbi.at/svnsimplex/CommonBioCommands/tags/simplex-1.0/CommonBioCommands/src/at/tugraz/genome/biojava/cli/cmd/io/ann/MergeGATKAndAnnovarToAMTFormatCSCommand.java


my $annotationtable_ref = mergeGeneAnnot($annovarstem);
if ($annovarfilter ne 'none') {
	mergeFilterAnnot($annovarfilter, $annovarstem, $annotationtable_ref);
}

open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";

# print header line
print $output_handle "#CHROM\tPOS\tREF\tALT\tVarClass\tNearestGene\tAAchange\t";
my @filternames = split(",", $annovarfilter);
if ($annovarfilter ne 'none') {
	foreach (@filternames) {
		print $output_handle "$_\t";
	}
}
print $output_handle "\n";


# my %annotationtable = %{$annotationtable_ref};
foreach my $desiredchr ((1..22), "X", "Y" ,"M") {
	if (defined ${$annotationtable_ref}{$desiredchr}) {						
		foreach my $pos ( sort {$a <=> $b} keys %{ ${$annotationtable_ref}{$desiredchr} } ) {
			my %variants = %{ ${$annotationtable_ref}{$desiredchr}{$pos} };
			foreach my $alleles (keys %variants) {
				my ($ref,$alt) = split("/", $alleles);
				my $geneannotation = $variants{$alleles}{"gene"};
				print $output_handle "$desiredchr\t$pos\t$ref\t$alt\t$geneannotation";
				if ($annovarfilter ne 'none') {
					foreach my $filter (@filternames) {
						my $filterannotation;
						if (!defined $variants{$alleles}{$filter}) {
							$filterannotation = '.';
						} else {
							$filterannotation = $variants{$alleles}{$filter};
						}
						print $output_handle "\t$filterannotation";
					}
				}
				print $output_handle "\n";
			}
		}
	}
}
close $output_handle;







sub mergeFilterAnnot {
	my ($annovarfilter, $annovarstem, $annotationtable_ref) = @_;
	my @filternames = split(",", $annovarfilter);
	
	foreach my $filtername (@filternames) {
		my $filedropped = "$annovarstem.hg19_$filtername"."_dropped";
		my $filefiltered = "$annovarstem.hg19_$filtername"."_filtered";
		
		my $filewithfilterdata;
		if (checkFileContainsFilterVals($filedropped, $filtername) == 1)  {
			$filewithfilterdata = $filedropped;
		} elsif (checkFileContainsFilterVals($filefiltered, $filtername) == 1) {
			$filewithfilterdata = $filefiltered;
		} else {
			die "Neither $filedropped nor $filefiltered contains values for $filtername\n";
		}
		
		open (my $filterdata_handle, "$filewithfilterdata");
		while (<$filterdata_handle>) {
			$_ =~ s/\s+$//;					# Remove line endings		
			my ($matchfiltername, $filterval, $chr, $annovar_start, $annovar_stop, $annovar_ref, $annovar_alt, $annovar_length, $vcf_start, $vcf_stop, $vcf_ref, $vcf_alt) = split ("\t", $_);
			${$annotationtable_ref}{$chr}{$vcf_start}{"$vcf_ref/$vcf_alt"}{$filtername} = $filterval;
		}
		close $filterdata_handle;
	}
}

sub mergeGeneAnnot {
	my $annovarstem = $_[0];
	my $allvarfile = "$annovarstem.variant_function";
	my $exonicfile = "$annovarstem.exonic_variant_function";
	
	my %annotations;
	open (ALL, "$allvarfile") or die "Cannot open $allvarfile file.\n";
	open (EXON, "$exonicfile") or die "Cannot open $exonicfile file.\n";
	while ( <ALL> ){
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split ("\t", $_);
		my ($generalclass, $nearestgene, $chr, $annovar_start, $annovar_stop, $annovar_ref, $annovar_alt) = @line[0..6];
		my ($vcf_start, $vcf_ref, $vcf_alt) = ($line[8], $line[10], $line[11]);
		my ($exonclass, $aminochange, $exonchr, $exonpos) = (('.')x5);
		if ($generalclass =~ /^exonic[;,]*/) {
			my $newline = <EXON>;
			my @exonicinfo = split("\t", $newline);
			($exonclass, $aminochange) = ($exonicinfo[1], $exonicinfo[2]);
			$aminochange =~ s/,$//;
			my ($exonchr, $exonstart, $exonstop) = ($exonicinfo[3], $exonicinfo[4], $exonicinfo[5]);
			if ($exonstart != $annovar_start) {
				print STDERR "\nError: next line in exonic_variant_function doesn't match variant_function\n";
				print STDERR "variant_function: @line\n";
				print STDERR "exonic_variant_function: @exonicinfo\n";
				die;
			}
		} 
		if ($exonclass eq '.') {
			$exonclass = $generalclass;
		}
		if ($nearestgene =~ /\(/) {
			$nearestgene =~ s/\((.*)\)//;
			$aminochange = $1;
		}
		$annotations{$chr}{$vcf_start}{"$vcf_ref/$vcf_alt"}{"gene"} = "$exonclass\t$nearestgene\t$aminochange";	
	}
	close EXON;
	close ALL;
	
	return (\%annotations);
}

sub checkFileContainsFilterVals {
	my ($filename, $filtername) = @_;
	my $sampleline = `head -1 $filename`;
	my @samplevals = split("\t", $sampleline);
	my $return = 0;
	if ($samplevals[0] eq $filtername) {
		$return = 1;
	}
	return $return;
}



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


combine_annovar_annotations.pl - 


=head1 SYNOPSIS


perl B<combine_annovar_annotations.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--annovarstem> F<annovar filestem>	

	stem name of all annovar annotations

=item B<--annovarfilters> F<annovar filers>	

	comma-separated list of filter names annovar added to files, such as "gerp++gt2" or "avsift"

=item B<--out> F<output file>

	name of output file

=item B<--help> I<help>

	print documentation

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
