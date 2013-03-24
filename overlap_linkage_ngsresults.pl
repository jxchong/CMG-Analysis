#!/usr/bin/env perl
#
# Description: Get all candidate variants overlapping regions with linkage LOD > cutoff
#
#
#
# Created by Jessica Chong on 2013-03-23

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX;

my ($candidatefile, $linkageresults, $linkagemap, $outputfile, $minlod, $help);

GetOptions(
	'linkage=s' => \$linkageresults, 
	'map=s' => \$linkagemap,
	'candidates=s' => \$candidatefile,
	'out=s' => \$outputfile,
	'minlod=f' => \$minlod,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
	pod2usage(-verbose=>1, -exitval=>1) if $help;

if (!defined $linkageresults) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --linkage not defined.\n")
} elsif (!defined $candidatefile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --candidates not defined.\n")
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined.\n")
} elsif (!defined $minlod) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --minlod not defined.\n")
} elsif (!defined $linkagemap) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --map not defined.\n")
}


print "Reading map and linkage results\n";
my %snpname2pos;
open (FILE, "$linkagemap") or die "Cannot read $linkagemap file: $!.\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $label, $cm) = split ("\t", $_);
	my $pos = ($cm * 10**6);
	$snpname2pos{$label} = $pos;
	# print "adding $chr,$label,$cm => $pos to array\n";
}
close FILE;

my (%lookuplod, %lookuppos);
open (FILE, "$linkageresults") or die "Cannot read $linkageresults file: $!.\n";
<FILE>;
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $cm, $label, $model, $lod, @rest) = split ("\t", $_);
	my $pos = $snpname2pos{$label};
	push(@{$lookuppos{$chr}}, $pos);
	push(@{$lookuplod{$chr}}, $lod);    
}
close FILE;


my $filteredout = 0;

print "Checking candidate variants file\n";
open (my $outputfile_handle, ">", $outputfile) or die "Cannot write to $outputfile: $!.\n";
open (my $candidatefile_handle, "$candidatefile") or die "Cannot read $candidatefile file: $!.\n";
my $currchr = 0;
my (@lookupposchr, @lookuplodchr);
while ( <$candidatefile_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	if ($_ =~ /^#/ || $_ =~ /^chr\tpos/) { 
		print $outputfile_handle "$_\tLOD\tLinkageExclude\n";
		next;										# skip header line
	}
	my ($chr, $pos, @line) = split ("\t", $_);

	# print "\nVariant: $chr:$pos\n";
	
	if ($chr ne $currchr) {
		print "Loading data for chromosome $chr\n";
		$currchr = $chr;
		if ($chr ne 'X' && $chr ne 'Y') {
			@lookupposchr = @{$lookuppos{$chr}};
			@lookuplodchr = @{$lookuplod{$chr}};
		}
	}
		
	if ($chr ne 'X' && $chr ne 'Y') {
		my ($lodleft, $lodright) = getNearbyLOD($pos, \@lookupposchr, \@lookuplodchr);
		if ($lodleft >= $minlod || $lodright >= $minlod) {
			print $outputfile_handle "$chr\t$pos\t".join("\t", @line)."\t$lodleft,$lodright\tN\n";
		} else {
			print $outputfile_handle "$chr\t$pos\t".join("\t", @line)."\t$lodleft,$lodright\tY\n";
			$filteredout++;
		}
	} else {
		print $outputfile_handle "$chr\t$pos\t".join("\t", @line)."\tNA,NA\tNA\n";
	}
}
close $candidatefile_handle;
close $outputfile_handle;


print "Filtered out $filteredout with LOD < $minlod\n";




sub getNearbyLOD {
	my ($targetpos, $lookuppos_ref, $lookuplod_ref) = @_;
	my @posarray = @{$lookuppos_ref};
	
	my $startpoint = 0;
	my $endpoint = $#posarray;
	my $midpoint;
 	while ($endpoint-$startpoint > 1) {
		$midpoint = floor(($startpoint+$endpoint) / 2);
		# print "mid=$midpoint ($posarray[$midpoint] bp), start=$startpoint, end=$endpoint\n";
	
		if ($targetpos > $posarray[$midpoint]) {
			$startpoint = $midpoint;
		} elsif ($targetpos < $posarray[$midpoint]) {
			$endpoint = $midpoint;
		} elsif ($targetpos == $posarray[$midpoint]) {
			$startpoint = $midpoint;
			$endpoint = $midpoint;
		}
	}
	
	# print "STOP: mid=$midpoint, start=$startpoint, end=$endpoint\n";
	my ($lodleft, $lodright) = (${$lookuplod_ref}[$startpoint], ${$lookuplod_ref}[$endpoint]);	
	return ($lodleft, $lodright);
}





################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME




=head1 SYNOPSIS


perl B<.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--candidates> I<candidate variants input file>	

	file with first column=chr, second column=position

=item B<--out> I<output file>

	name of output file

=item B<--linkage> I<.tbl linkage analysis output from Merlni>

	.tbl file output by Merlin after running linkage analysis
	
=item B<--minlod> I<minimum LOD>

	variant must be in region with a minimum LOD of at least this value to be included in the output
	
=item B<--map> I<map file>

	map file used in running linkage analysis

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


This script will 


=head1 CAVEATS


xx

=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
