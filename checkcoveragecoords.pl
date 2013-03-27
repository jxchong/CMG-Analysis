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


my ($bedfile, $coveragefile, $outputfile, $outputdir, $help, $mindepthval);

GetOptions(
	'mindepth:f' => \$mindepthval,
	'bedfile=s' => \$bedfile,
	'coverage=s' => \$coveragefile, 
	'outdir=s' => \$outputdir,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $bedfile) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --bedfile not defined.\n")
} elsif (!defined $coveragefile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --coverage not defined\n");
} elsif (!defined $mindepthval) {
	$mindepthval = 6;
} elsif (!defined $outputdir) {
	$outputdir = '.';
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
}


my @splitname = split("/", $coveragefile);
my $coveragestem = $splitname[$#splitname];
$coveragestem =~ s/\.out//;

if ($coveragefile =~ /all\.coverage\.out$/ && ! -e "$outputdir/$coveragestem.tsv.gz") {
	my $convertedfile = "$coveragestem.tsv.gz";
	my $cmd = q/if ($F[0] eq "Locus") {print "Chr\tPos\t".join("\t", @F[1..$#F])."\n";} else { ($chr,$bp)=split(":", $F[0]); print "$chr\t$bp\t".join("\t", @F[1..$#F])."\n";}/;
	print "Converting $coveragefile to $outputdir/$convertedfile\n";
	if (!-e "$outputdir/$convertedfile") {
		`perl -ane \'$cmd\' $coveragefile | bgzip > $outputdir/$convertedfile`;
	}
	print "Tabix indexing $outputdir/$convertedfile\n";
	`tabix -s 1 -b 2 -e 2 $outputdir/$convertedfile`;
} else {
	open (my $output_handle, ">", "$outputdir/$outputfile") or die "Cannot write to $outputdir/$outputfile: $!.\n";	
	my $headerline = `zcat $coveragestem.tsv.gz | head -1`;
	$headerline =~ s/\s+$//;
	my @header = split("\t", $headerline);
	print $output_handle join("\t", @header[0..3])."\tName\tMin\tMax\t".join("\t", @header[4..$#header])."\n";
	
	open (my $input_handle, "$bedfile") or die "Cannot read $bedfile: $!.\n";
	while ( <$input_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $start, $stop, @bedline) = split(/\s+/, $_);
		my @coveragedata = `tabix $coveragestem.tsv.gz $chr:$start-$stop\n`;
		my ($countsites, $countbelowdepth) = (0,0);
		foreach my $line (@coveragedata) {
			$line =~ s/\s+$//;					# Remove line endings
			$countsites++;
			my @data = split("\t", $line);
			my ($min, $max) = ($data[4], 0);
			my $countsubjectsbelowdepth = 0;
			for (my $i=4; $i<=$#data; $i++) {
				if ($data[$i] < $min) { $min = $data[$i]; }
				if ($data[$i] > $max) { $max = $data[$i]; }
				if ($data[$i] < $mindepthval) {
					$countsubjectsbelowdepth++;
				}
			}
			if ($countsubjectsbelowdepth > 0) { $countbelowdepth++; }
			print $output_handle join("\t", @data[0..3]);
			if (defined $bedline[0]) {
				print $output_handle "\t$bedline[0]\t";
			} else {
				print $output_handle "\t.\t";
			}
			print $output_handle "$min\t$max\t".join("\t", @data[4..$#data])."\n";
		}
		print "At least one subject had DP<$mindepthval at $countbelowdepth out of $countsites sites in $chr:$start-$stop\n";
	}
	close $input_handle;
	close $output_handle;
}






################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


xxx.pl - 


=head1 SYNOPSIS


perl B<checkcoveragecoords.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--mindepth> I<minimum depth>	

	if depth less than this value in at least one subject, alert

=item B<--bed> F<coordinate file>	

	file of coordinates to fetch coverage for, in bed format

=item B<--coverage> F<coverage file>	

	coverage input file

=item B<--out> F<output file>

	name of output file

=item B<--outdir> I<output directory>

	path to output directory

=item B<--help> I<help>

	print documentation

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


Takes the standards UWCMG coverage.out file, optionally reformats it into a bgzip/tabix-indexable file, and given a UCSC BED-format file with coordinates, retrieves the per-subject depth for all positions between those coordinates.  Given a mindepth value, will count number of sites where at least one subject had a depth below that value.


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
