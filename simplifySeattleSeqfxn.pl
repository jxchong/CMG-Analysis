#!/usr/bin/env perl
#
# Description:
#
#
#
# Created by Jessica Chong on 20xx-xx-xx.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($inputfile, $outputfile, $fxncol, $help);

GetOptions(
	'in=s' => \$inputfile, 
	'desiredgenefile=s' => \$desiredgenefile,
	'desiredgenescol=i' => \$desiredgenecol,
	'fxncol=s' => \$fxncol,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $inputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} elsif (!defined $fxncol) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --fxncol not defined\n");
} elsif (!defined $desiredgenefile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --desiredgenefile not defined\n");
} elsif (!defined $desiredgenescol) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --desiredgenescol not defined\n");
}

my %gene2accession;
my $accessionfile = '~/references/locus_info.137.txt.gz';
open (my $input_handle, "zcat $accessionfile |") or die "Cannot read $accessionfile: $!.\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;
	my @line = split("\t", $_);
	# gene => accession
	push(@$gene2accession{$line[1]}, $line[8]);
close $input_handle;


my %desiredaccessions;
open (my $input_handle, "$desiredgenefile") or die "Cannot read $desiredgenefile: $!.\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;
	my @line = split("\t", $_);
	my $desiredgene = $line[$desiredgenecol-1];
	# accessions => gene
	$desiredaccessions{$gene2accession{$desiredgene}} = $desiredgene;
close $input_handle;



open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (my $input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
my $header = <$input_handle>;
$header =~ s/\s+$//;
my @headercols = split("\t", $header);
for (my $i=0; $i<$#headercols; $i++) {
	my $hcol = $headercols[$i];
	if ($hcol =~ /FG|INFO|functionGVS/i) {
		print $output_handle "$hcol\tFxnSeverity\tWorstFxn\t";
	} else {
		print $output_handle "$hcol\t";
	}
}
print $output_handle "$headercols[$#headercols]\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split("\t", $_);
	
	my $accessionstring = getfieldstring($line[$fxncol-1], "GM");			# only works for vcfs
	my $fxnstring = getfieldstring($line[$fxncol-1], "FG");
	
	my ($worstfxnval, $worstfxn, $worstgene) = getSSfxn($fxnstring, $accessionstring);
	print $output_handle join("\t", @line[0..($fxncol-1)])."\t";
	print $output_handle "$worstfxnval\t$worstfxn\t";
	print $output_handle join("\t", @line[($fxncol)..$#line])."\n";
}
close $input_handle;
close $output_handle;



sub getSSfxn {
	my $functionimpact = $_[0];
	my $accessionstring = $_[1];
	
	my %fxnvalues = map {$_ => 0} qw(intergenic near-gene-3 near-gene-5 none);
	@fxnvalues{qw(intron utr utr-3 utr-5)} = (1) x 4;
	@fxnvalues{qw(coding-synonymous)} = 2;
	@fxnvalues{qw(coding-synonymous-near-splice)} = 3;
	@fxnvalues{qw(missense missense-near-splice)} = (4) x 2;
	@fxnvalues{qw(stop-lost stop-lost-near-splice)} = (5) x 2;
	@fxnvalues{qw(splice-3 splice-5)} = (6) x 2;
	@fxnvalues{qw(coding)} = 7;
	@fxnvalues{qw(stop-gained stop-gained-near-splice)} = (8) x 2;
	@fxnvalues{qw(coding-near-splice)} = 9;
	@fxnvalues{qw(codingComplex coding-Complex-near-splice)} = (10) x2;
	@fxnvalues{qw(frameshift frameshift-near-splice)} = (11) x 2;
	
	my @filters = split(",", $functionimpact);
	my @accessions = split(",", $accessionstring);
	my $worstfxnval = -1;
	my $worstfxn;
	
	for (my $i=0; $i<=$#filters; $i++) {
		my $filter = $filters[$i]; 
		my $accession = $accessions[$i];
		# print "fxn=$eachfilter, value=$fxnvalues{$eachfilter}\n";
		if ($fxnvalues{$eachfilter} > $worstfxnval) {
			$worstfxnval = $fxnvalues{$eachfilter};
			$worstfxn = $eachfilter;
		}
	}

	return ($worstfxnval, $worstfxn, $worstgene);
}


sub getfieldstring {
	my $columndata = $_[0];
	my $desiredfield = $_[1];
	my $string;
	
	if ($columndata =~ /;$desiredfield=/) {		# must be a vcf-derived column, split out FG value
		my @columncontents = split(";", $columndata);
		foreach my $column (@columncontents) {
			my ($fieldname, $value) = split("=", $column);
			if (!defined $value) {
				$value = '.';
			}
			if ($fieldname eq 'FG') {
				$string = $value;
			}
		}		
	} else {							# assume column must be only the desired one, derived from tab-delimited SSAnnotation or SeattleSeqAnnotation file or other
		$string = $columndata;
	}
	return $string;
}






################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


simplifySeattleSeqfxn.pl - 


=head1 SYNOPSIS


perl B<simplifySeattleSeqfxn.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> F<input file>	

	input file, output from SeattleSeq, assumes tab-delimited

=item B<--fxncol> F<column number (starting from 1)>	

	column containing SeattleSeq function annotation (column number for either the INFO column from vcf or functionGVS column)

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
