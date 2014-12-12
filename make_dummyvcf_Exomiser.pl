#!/usr/bin/env perl
#
# Description: Given sample ID and familyID (value in FamilieswHits), make a dummy single-sample vcf file from tab-delimited file produced by filter_model.pl.
#
#
#!/
# Created by Jessica Chong on 2014-12-10.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($tsvinputfile, $outputfile, $keepsampleID, $keepfamilyID, $help);

GetOptions(
	'tsv=s' => \$tsvinputfile, 
	'keepsample=s' => \$keepsampleID,
	'keepfamily=s' => \$keepfamilyID,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $tsvinputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} elsif (!defined $keepsampleID) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --keepsample not defined\n");
} elsif (!defined $keepfamilyID) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --keepfamily not defined\n");
}


# copy dummy vcf header to output
`cp /nfs/home/jxchong/references/dummyheader.vcf $outputfile`;


open (my $output_handle, ">>", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$keepsampleID\n";
open (my $input_handle, "$tsvinputfile") or die "Cannot read $tsvinputfile: $!.\n";
my $headerline = <$input_handle>;
$headerline =~ s/\s+$//;
my @header = split("\t", $headerline);
my $samplegenotypecol = 'NA';
for (my $i=0; $i<=$#header; $i++) {
	my $columnname = $header[$i];
	if ($keepsampleID eq $columnname) {
		$samplegenotypecol = $i;
	}
}
if ($samplegenotypecol eq 'NA') {
	die "Could not find sample with ID=$keepsampleID in input $tsvinputfile\n";
}

while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split("\t", $_);
	my ($chr, $start, $stop, $type, $ref, $alt, $filter, $genelist, $rsID) = @line[0..9];
	my $FamilieswHits = $line[$#line];
	my $genotype = call2vcfgeno($line[$samplegenotypecol], $ref, $alt);
	my $qual = 100;
	my $info = '.';
	my $format = 'GT';
	if ($keepfamilyID eq $FamilieswHits) {
			print $output_handle "$chr\t$start\t$rsID\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$genotype\n";
	}
}
close $input_handle;
close $output_handle;






sub call2vcfgeno {
	my ($genotype, $ref, $alt) = @_;
	my $vcfgenotype;
	my @allelearray = ($ref, split(",", $alt));
	my $allelecount = 0;
	my %allalleles = map { $_ => $allelecount++} @allelearray;
	if ($genotype eq 'N') {
		$vcfgenotype = './.';
	} else {
		my ($allele1, $allele2) = split("/", $genotype);
		$vcfgenotype = "$allalleles{$allele1}/$allalleles{$allele2}";
	}
	return $vcfgenotype;
}



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


xxx.pl - 


=head1 SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--tsv> F<input tsv file>	

	input tsv file created by filter_model.pl

=item B<--out> F<output file>

	name of output file

=item B<--keepfamily> F<string>

	string to extract from FamilieswHits column

=item B<--keepsample> F<string>

	sampleID column to keep from vcf

=item B<--help> I<help>

	print documentation

=back


=head1 FILES


xx


=head1 EXAMPLES


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
