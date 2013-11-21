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


my ($inputfile, $outputfile, $chrcol, $poscol, $refcol, $altcol, $help);

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $inputfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} 


my $caddfile = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/CADD/v1.0/whole_genome_SNVs.tsv.gz';

open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (my $input_handle, "$inputfile") or die "Cannot read $inputfile: $!.\n";
while ( <$input_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	print $output_handle "$_";
	
	if ($_ =~ /^#/) {
		if ($_ =~ /^#chr/i) {
			my @header = split("\t", $_);
			for (my $i=0;$i<=$#header;$i++) {
				if ($header[$i] =~ /#chr/i) { $chrcol = $i; }
				if ($header[$i] =~ /start/i) { $poscol = $i; }
				if ($header[$i] =~ /ref/i) { $refcol = $i; }
				if ($header[$i] =~ /alt/i) { $altcol = $i; }
				if ($i >= 6) {last;}
			}
			# print "Found chrcol=$chrcol, poscol=$poscol, refcol=$refcol, altcol=$altcol\n";
			print $output_handle "\tCADDphred";
		}
	} else {
		my @line = split("\t", $_);
		$line[0] =~ s/chr//;
		# print "Executing tabix $caddfile $line[$chrcol]:$line[$poscol]-$line[$poscol]\n";
		my @cadd_data = `tabix $caddfile $line[$chrcol]:$line[$poscol]-$line[$poscol]` or die "Cannot execute tabix\n";
		my $phred = "NA";
		
		foreach my $scoreline (@cadd_data) {
			$scoreline =~ s/\s+$//;
			my ($chr,$pos,$ref,$alt,$raw,$caddphred) = split("\t", $scoreline);
			if ($ref eq $line[$refcol] && $alt eq $line[$altcol]) {
				$phred = $caddphred;
			}
		}
		print $output_handle "\t$phred";
	}
	print $output_handle "\n";
}
close $input_handle;
close $output_handle;




################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


xxx.pl - 


=head1 SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> F<input file>	

	input file

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
