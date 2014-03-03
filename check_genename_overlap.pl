#!/usr/bin/env perl
#
# Description: Check for overlap in gene names between a gene list file containing a candidate list and the results from filtering using filter_model.pl
#
#
#
# Created by Jessica Chong on 20xx-xx-xx.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($inputfile, $outputfile, $genelistfile, $help);

GetOptions(
	'in=s' => \$inputfile, 
	'genelist=s' => \$genelistfile,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $inputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} 


my %genelist;
open (my $genelist_handle, "$genelistfile") or die "Cannot read $genelistfile: $!.\n";
while ( <$genelist_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	$genelist{$_} = 1;
}
close $genelist_handle;


my $genecol;
open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (my $input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
my $headerline = <$input_handle>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
for (my $i=0;$i<=$#header;$i++) {
	if ($header[$i] eq 'geneList') { $genecol = $i; }
}
print $output_handle join("\t", @header)."\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split("\t", $_);
	my @genes = split(",", $line[$genecol]);
	my $include = 0;
	foreach my $gene (@genes) {
		if (exists $genelist{$gene}) { $include = 1; }
	}
	if ($include == 1) {
		print $output_handle join("\t", @line)."\n";
	}
}
close $input_handle;
close $output_handle;




################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


check_genename_overlap.pl - Check for overlap in gene names between a gene list file containing a candidate list and the results from filtering using filter_model.pl


=head1 SYNOPSIS


perl B<check_genename_overlap.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> F<input file>	

	input file (output from filter_model.pl script)

=item B<--genelist> F<gene list file file>

	just a list of genes in the first column

=item B<--out> F<output file>

	name of output file

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
