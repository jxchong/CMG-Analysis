#!/usr/bin/env perl
#
# Description: Look for mutations matching a particular model given an SSAnnotation file produced by SeattleSeq.
#
#
#
# Created by Jessica Chong on 2012-10-15.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, $inheritmodel, $mafcutoff, $excludeGVSfunction);
my $cmgfreqcutoff = 0.2;

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'subjectreq=s' => \$subjectdeffile,
	'minhits=i' => \$minhits,
	'GATKkeep=s' => \$filters,
	'N=s' => \$isNhit,
	'model:s' => \$inheritmodel,
	'mafcutoff:f' => \$mafcutoff,
	'excludefunction=s' => \$excludeGVSfunction,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $subjectdeffile) {
	optionUsage("option --subjectreq not defined\n");
} elsif (!defined $minhits) {
	optionUsage("option --minhits not defined\n");
} elsif (!defined $filters) {
	optionUsage("option --GATKkeep not defined\n");
} elsif (!defined $isNhit) {
	optionUsage("option --N not defined\n");
} elsif (!defined $mafcutoff) {
	$mafcutoff = 0.01;
} elsif (!defined $inheritmodel) {
	$inheritmodel = 'NA';
} elsif (!defined $excludeGVSfunction) {
	optionUsage("option --excludefunction not defined\n");
}


# if inputfile has "wfreqs" in filename, assume add_error_populationfreqs.pl has already been run and values added to last two columns
my $filehasfreqs = 0;
if ($inputfile =~ /wfreqs/i) {
	$filehasfreqs = 1;
	print STDOUT "$inputfile contains frequency of variants in CMG and outside populations\n";
} else {
	print STDOUT "add_error_populationfreqs.pl was not run on $inputfile, so filtering based on frequency of variant in other samples will be much slower\n";
}

my %allowedGATKfilters;
if ($filters eq 'all' || $filters eq 'any') {
	%allowedGATKfilters = map {$_ => 1} qw(QUALFilter QDFilter LowQual SBFilter PASS ABFilter LowQual HRunFilter SnpCluster QUALFilter);
} else {
	%allowedGATKfilters = map {$_ => 1} split(',', $filters);
}

my %GVStoexclude;
if ($excludeGVSfunction eq 'default') {
	%GVStoexclude = map {$_ => 1} qw(intron intergenic coding-synonymous utr-3 utr-5 near-gene-3 near-gene-5);
} else {
	%GVStoexclude = map {$_ => 1} split(',', $excludeGVSfunction);
}



my $logfile = "$outputfile.log";
open (LOG, ">$logfile") or die "Cannot write to $logfile: $!.\n";
print LOG "input=$inputfile\n";
print LOG "output=$outputfile\n";
print LOG "subjectreq=$subjectdeffile\n";

my %countuniquefamilies_hash;
my @orderedsubjects;
my %subjects;
open (SUBJECTS, "$subjectdeffile") or die "Cannot read $subjectdeffile: $!.\n";
while (<SUBJECTS>) {
	$_ =~ s/\s+$//;					# Remove line endings
	print LOG "$_\n";
	my ($familyid, $subjectid, $father, $mother, $relation, $desiredgeno) = split("\t", $_);
	$subjects{$subjectid} = [$familyid, $father, $mother, $relation, $desiredgeno];
	push(@orderedsubjects, $subjectid);
	if ($subjectid !~ '#') {
		$countuniquefamilies_hash{$familyid} = 1;
	} else {
		print LOG "$subjectid is being skipped in this analysis\n";
	}
}
close SUBJECTS;
my $countuniquefamilies = scalar(keys %countuniquefamilies_hash);
print LOG "\n";


my %genehits;
my $countinputvariants = 0;
my $printparams = 0;
open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
if ($inputfile =~ m/.vcf/ && $inputfile =~ m/.gz/) {
	open (FILE, "zcat $inputfile |") or die "Cannot read $inputfile file: $!.\n";
} else {
	open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
}
my $headerline = <FILE>;
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
my @subjectcolumns;
my ($polyphencol, $phastconscol, $gerpcol, $gatkfiltercol);
my $isannotated = 0;
for (my $i=0; $i<=$#header; $i++) {
	my $columnname = $header[$i];
	if ($columnname =~ /Gtype/i) {
		$columnname =~ s/Gtype//; 
	}
	if (defined $subjects{$columnname} || defined $subjects{"#$columnname"}) {
		push(@subjectcolumns, $i);
	}
	if ($columnname =~ /Freqin/i) {
		$isannotated = 1;
	} elsif ($columnname =~ /polyPhen/i) {
		$polyphencol = $i;
	} elsif ($columnname =~ /GERP/i) {
		$gerpcol = $i;
	} elsif ($columnname =~ /PhastCons/i) {
		$phastconscol = $i;
	} elsif ($columnname =~ /filterFlagGATK/i) {
		$gatkfiltercol = $i;
	}
}

if ($isannotated == 1) {
	print LOG "File contains frequencies in outside populations and frequency observed in other CMG subjects\n";
} else {
	print LOG "!!!!!! File doesn't contain frequencies in outside populations and frequency observed in other CMG subjects\n";
	print LOG "Run: perl ~/bin/add_error_populationfreqs.pl --in $inputfile --out <outputfile> --capture <capture array for systematic error filtering>\n";
	die;
}

print OUT "$headerline\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	$countinputvariants++;
	
	if ($printparams == 0) {
		print LOG "Requiring hits in gene in at least $minhits subjects/families\n";
		print LOG "Missing genotypes/no calls are counted as: $isNhit\n";
		print LOG "Excluding all variants with annotations: ".join(" ", keys %GVStoexclude)."\n";
		print LOG "Only allow variants with GATK filter: ".join(" ", keys %allowedGATKfilters)."\n";
		print LOG "Excluding variants with MAF>$mafcutoff in ESP and/or 1000 Genomes\n";
		print LOG "Excluding variants with frequency>=$cmgfreqcutoff in CMG subjects (likely systematic error)\n";
		print LOG "Compound het analysis? $inheritmodel\n";
		$printparams = 1;
	}
	
	my @line = split ("\t", $_);
	my $countmatches = 0;
	my ($chr, $pos, $vartype, $ref, $alt);
	my @subjectgenotypes;
	my ($gene, $functionimpact, $gerp, $polyphen, $phastcons);
	# my $indbsnp = 'NA';
	# my $inUWexomes = 'NA';
	# my $UWexomescovered = 'NA';
	my $filterset;
	
	if ($inputfile =~ m/SeattleSeqAnnotation134/) {
		if ($inputfile =~ m/snps/) {
			my $sampleallelestring = $line[5];
			my @samplealleles = split("/", $sampleallelestring);
			if (!defined $samplealleles[1]) {
				$samplealleles[1] = 'N';
				next;
			}
			$gene = $line[20];
			$filterset = 'NA';
			$functionimpact = $line[8];
			# $indbsnp = $line[0];
			$gerp = $line[$gerpcol];
			$phastcons = $line[$phastconscol];
			$polyphen = $line[$polyphencol];		
			$vartype = 'SNP';
			($ref, $alt) = selectRefAlt($line[3], $samplealleles[0], $samplealleles[1]);
			@subjectgenotypes = split(",", $line[4]);
		}
		if ($inputfile =~ m/indels/) {
			# ($vartype, $ref, $alt) = ('indel', $line[3], $line[4]);		# complicated parsing; not implemented
			die "Have not implemented parsing of this file format yet\n";
		}
	} elsif ($inputfile =~ m/SSAnnotation/) {
		($chr, $pos, $vartype, $ref, $alt) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		$filterset = $line[$gatkfiltercol];
		$gene = $line[7];
		@subjectgenotypes = @line[@subjectcolumns];
		# $indbsnp = $line[8];
		$functionimpact = $line[10];
		$phastcons = $line[$phastconscol];
		$gerp = $line[$gerpcol];
		$polyphen = $line[$polyphencol];
		# $inUWexomes = $line[16];
		# $UWexomescovered = $line[17];
	} else {
		die "Input file ($inputfile) isn't an SSAnnotation or SeattleSeqAnnotation134 file\n";
	} 
	
	# Check that user input corresponds to input data files
	if (scalar(@subjectgenotypes) != scalar(@orderedsubjects)) {
		die "Your input subject definition file ($subjectdeffile) lists a different number of subjects (".scalar(@orderedsubjects).") than are contained (".scalar(@subjectgenotypes).") in the input data file ($inputfile)\n";
	}
	
	# if (isNovel($indbsnp, $inUWexomes, $UWexomescovered) && shouldfunctionFilter(\%GVStoexclude, $functionimpact)==0 && ($filterset eq 'NA' || checkGATKfilters(\%allowedGATKfilters, $filterset))) {
	if (shouldfunctionFilter(\%GVStoexclude, $functionimpact)==0 && ($filterset eq 'NA' || checkGATKfilters(\%allowedGATKfilters, $filterset))) {
		my %checkfamilies;
		for (my $i=0; $i<=$#subjectgenotypes; $i++) {
			my $genotype = $subjectgenotypes[$i];
			my $subjectid = $orderedsubjects[$i];
			if ($subjectid !~ '#') {
				my ($familyid, $father, $mother, $relation, $desiredgeno) = @{$subjects{$subjectid}};
				if (!exists $checkfamilies{$familyid}) {
					$checkfamilies{$familyid}{'father'} = -1;
					$checkfamilies{$familyid}{'mother'} = -1;
					$checkfamilies{$familyid}{'child'} = -1;
				}
				my $ismatch += checkGenoMatch($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit);
				
				if ($relation eq 'father') {
					$checkfamilies{$familyid}{'father'} = $ismatch; 
				}
				if ($relation eq 'mother') {
					$checkfamilies{$familyid}{'mother'} = $ismatch; 
				}
				if ($relation eq 'child') {
					$checkfamilies{$familyid}{'child'} = $ismatch;
				}
			} 
		}

		my $countfamiliesmatch = 0;
		my %matchingfamilies;
		my %matchingtrios;
		while (my ($familyid, $thisfamily_ref) = each %checkfamilies) {
			my %thisfamily = %{$thisfamily_ref};

			my $thisfamilymatch = 0;
			my $familysize = grep { $thisfamily{$_} != -1 } keys %thisfamily;	
			# print "family $familyid is size=$familysize\n";		
			if ($inheritmodel eq 'compoundhet' && $familysize >= 3) {							# assumes trios; doesn't handle more than one affected child per family (yet)
				if (($thisfamily{'father'}+$thisfamily{'mother'}) == 1 && $thisfamily{'child'} == 1) {
					$countfamiliesmatch += 1;
					my $matchsubj;
					if ($thisfamily{'father'} == 1) {
						$matchsubj = 'father';
					}
					if ($thisfamily{'mother'} == 1) {
						$matchsubj = 'mother';
					}
					$matchingtrios{$familyid} = $matchsubj;	
					# print "$gene, chr$chr:$pos = hit for trio $familyid\n";
				}
			} else {
				while (my ($relation, $thismatch) = each %thisfamily) {
					if ($thismatch != -1) {
						$thisfamilymatch += $thismatch;
					}
				}

				if ($thisfamilymatch == $familysize) {
					$countfamiliesmatch += 1;
					$matchingfamilies{$familyid} = 1;
					# print "$gene, chr$chr:$pos = hit for family $familyid\n";
				}
			}
		}

		if ($countfamiliesmatch > 0) {
			my $data = join("\t", @line);
			push(@{$genehits{$gene}}, [$data, \%matchingfamilies, \%matchingtrios]);							# store this as a hit
		}
	}
	# if ($chr == 3) {
	# 	last;
	# }
}
close FILE;



# Check putative hit list for all genes and count number of valid hits
print "\nChecking all genes for desired number of hits in desired number of families\n";
my $countgeneswhits = 0;
my $countgenesrejectedhits = 0;
my $countoutputvariants = 0;
my $countkeptvariants = 0;
my $counterrorvariants = 0;
my $countcommonvariants = 0;
my $countexcludedvariants = 0;
my $genecounter = 0;
my $totalgenes = scalar(keys %genehits);
while (my ($gene, $results_ref) = each %genehits) {
	$genecounter++;
	# print "\n===========================================================================\n";
	if ($genecounter % 1000 == 0) {
		print "Checking gene #$genecounter (out of $totalgenes) for hits\n";
	}
	my @hitdata = @{$results_ref};
	
	# account for possible compound het model, then count hits in each family for this gene
	# do before filtering of common variants and systematic errors for increased efficiency (data access, even by tabix, is slower)
	my $resultsFamiliesvsModel_ref = checkFamiliesvsModel(\@hitdata, $inheritmodel);
	# review all families for required number of hits in this gene
	my $enoughfamilieshavehits = checkFamiliesforHits($resultsFamiliesvsModel_ref, $inheritmodel, $countuniquefamilies, $minhits);													
	
	if ($enoughfamilieshavehits == 1) {
		# filter out common variants and systematic errors
		my @filteredoutput;
		foreach my $hit (@hitdata) {
			my $hitvarinfo = ${$hit}[0];
				my @thishit = split("\t", $hitvarinfo);
				my ($chr,$pos,$vartype,$ref,$alt) = @thishit[0..4];
			my %matchingfamilies = %{${$hit}[1]};
			my %matchingtrios = %{${$hit}[2]};
			my $iserror = 0;
			my $iscommon = 0;
			if ($filehasfreqs == 1) {
				my $lastcol = $#thishit;
				my $freqinCMG = $thishit[($lastcol-1)]/100;							# storing allele freq as percentage
				my $freqinOutside = $thishit[$lastcol]/100;							# storing allele freq as percentage
				if ($freqinCMG >= $cmgfreqcutoff) {
					$iserror = 1;
				}
				if ($freqinOutside > $mafcutoff) {
					$iscommon = 1;
				}				
			} 
			if ($iserror==1) {
				$counterrorvariants++;
			}
			if ($iscommon==1) {
				$countcommonvariants++;
			}
			if ($iserror==0 && $iscommon==0) {
				push(@filteredoutput, $hit);
			} else {
				$countexcludedvariants++;
			}
		}
		
		# after filtering, recheck to make sure enough families still have enough hits in this gene, then output results
		my $resultsFamiliesvsModel_ref_postfilter = checkFamiliesvsModel(\@filteredoutput, $inheritmodel);
		my $enoughfamilieshavehits_postfilter = checkFamiliesforHits($resultsFamiliesvsModel_ref_postfilter, $inheritmodel, $countuniquefamilies, $minhits);													
		if ($enoughfamilieshavehits_postfilter == 1) {
			$countgeneswhits++;
			# print "... Printing hits for this gene ($gene)\n";
		  	foreach my $hit (@filteredoutput) {
		  		my $hitvarinfo = ${$hit}[0];
		  		$countoutputvariants++;
		  		print OUT "$hitvarinfo\n";
		  	}
		} else {
			$countgenesrejectedhits++;
			# print "... After filtering, no hits left for this gene ($gene)\n";
		}
	}
}

print LOG "\nResults summary:\n";
print LOG "In $countgeneswhits gene(s), matched $countoutputvariants variants\n";
print LOG "Total $countinputvariants variants examined from input\n";
print LOG "Total $countexcludedvariants additional possible hits in $countgenesrejectedhits gene(s) were excluded\n";
print LOG "N=$counterrorvariants are systematic errors and N=$countcommonvariants are common variants\n";

close LOG;



sub checkFamiliesvsModel {
	# account for possible compound het model, then recount matching families
	my ($inputhitdata_ref, $inheritmodel) = @_;
	my %familieswHitsinGene;
	foreach my $hit (@{$inputhitdata_ref}) {
		my $hitvarinfo = ${$hit}[0];
			my @thishit = split("\t", $hitvarinfo);
		my %matchingfamilies = %{${$hit}[1]};
		my %matchingtrios = %{${$hit}[2]};

		while (my($familyid, $ismatch) = each %matchingfamilies) {
			$familieswHitsinGene{$familyid}++;																										# DEBUG
		}

		if ($inheritmodel eq 'compoundhet') {
			while (my($familyid, $matchsubj) = each %matchingtrios) {
				if (!defined $familieswHitsinGene{$familyid}) {
					$familieswHitsinGene{$familyid}{'father'} = 0;
					$familieswHitsinGene{$familyid}{'mother'} = 0;
					$familieswHitsinGene{$familyid}{'child'} = 0;
				}
				$familieswHitsinGene{$familyid}{$matchsubj}++;
				$familieswHitsinGene{$familyid}{'child'}++;
			}
		}
	}
	return \%familieswHitsinGene;
}


sub checkFamiliesforHits {
	my($resultsFamiliesvsModel_ref, $inheritmodel, $countuniquefamilies, $minhits) = @_;
	my %familieswHitsinGene = %{$resultsFamiliesvsModel_ref};
	
	my $countfamiliesmatch = 0;
	while (my($familyid, $familyhits) = each %familieswHitsinGene) {																										# DEBUG
		if ($inheritmodel eq 'compoundhet') {
			if (ref($familyhits) eq 'HASH') {
				if (${$familyhits}{'child'} >= 2 && ${$familyhits}{'father'} >= 1 && ${$familyhits}{'mother'} >= 1) {
					$countfamiliesmatch++;
				}
			} elsif (!ref($familyhits)) {																										# DEBUG
				if ($familyhits >= 2) {
					$countfamiliesmatch++;
				}
			}
		} else {
			if ($familyhits >= 1) {
				$countfamiliesmatch++;
			}
		}
	}
	
	# if desired number of families have hits in gene
	if ($countfamiliesmatch >= $minhits) {
		return 1;
	} else {
		return 0;
	}
}

sub checkGATKfilters {
	my ($allowedfilters_ref, $filterset) = @_;
	my @thesefilters = split(';', $filterset);
	my $keep = 1;
	foreach my $filter (@thesefilters) {
		if (!defined ${$allowedfilters_ref}{$filter}) {
			$keep = 0;
		} 
	}
	return $keep;
}

# sub isNovel {
# 	my ($indbsnp, $inUWexomes, $UWexomescovered) = @_;
# 	my $score = 0;
# 	# if ($indbsnp eq 'none' || $indbsnp eq 'NA' || $indbsnp eq '.' || $indbsnp eq '0') {
# 	# 	$score += 1;
# 	# }
# 	if ($inUWexomes == 0 || $UWexomescovered == 0 || $inUWexomes eq 'NA') {
# 		$score += 1;
# 	}
# 	if ($score > 0) {
# 		return 1;
# 	} else {
# 		return 0;
# 	}
# }

sub shouldfunctionFilter {
	my ($GVStoexclude_ref, $functionimpact) = @_;
	my $score = 0;
	if (defined ${$GVStoexclude_ref}{$functionimpact}) {
		$score += 1;
	}
	if ($score > 0) {
		return 1;
	} else {
		return 0;
	}
}

sub checkGenoMatch {
	my ($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit) = @_;
		
	my @alleles;
	my $ismatch = 0;
	
	if (($genotype eq 'N' && $isNhit eq 'hit') || $desiredgeno eq 'any') {
		$ismatch = 1;
	} elsif ($genotype eq 'N' && $isNhit eq 'nohit') {
		$ismatch = 0;
	} else {
		if ($vartype eq 'SNP') {
			if ($genotype eq $ref) {
				@alleles = ($ref, $ref);
			} elsif ($genotype eq $alt) {
				@alleles = ($alt, $alt);
			} else {
				@alleles = ($ref, $alt);
			}
		} elsif ($vartype eq 'indel') {
			@alleles = split("/", $genotype);
		}
	
		# print "$vartype\t$genotype\t$desiredgeno\t$ref\t$alt\t@alleles\n";				# DEBUG
	
		if ($desiredgeno eq 'ref') {
			if ($alleles[0] eq $ref && $alleles[1] eq $ref) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'alt') {
			if ($alleles[0] eq $alt && $alleles[1] eq $alt) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'het') {
			if ($alleles[0] eq $ref && $alleles[1] eq $alt) {
				$ismatch = 1;
			}
		} 
	}
	return $ismatch;
}


sub selectRefAlt {
	my ($origref, $allele1, $allele2) = @_;
	my ($ref, $alt);
	if ($origref eq $allele1) {
		$ref = $allele1;
		$alt = $allele2;
	} else {
		$ref = $allele2;
		$alt = $allele1;
	}
	return $ref, $alt;
}




sub optionUsage {
	my $errorString = $_[0];
	print "$errorString";
	print "perl $0 \n";
	print "\t--in\tinput file\n";
	print "\t--out\toutput file\n";
	print "\t--subjectreq\tpedigree-like file listing families and required subject genotypes\n";
	print "\t--minhits\tminimum number of families/individuals with hits\n";
	print "\t--GATKkeep\tGATK quality filters that should be kept (comma-delimited with no spaces, or keep 'all')\n";
	print "\t--N\thit or nothit (how should we count missing genotypes)\n";
	print "\t--excludefunction\tcomma separated list of GVS variant function classes to be excluded, or 'default'\n";
		print "\t\tdefault=intron,intergenic,coding-synonymous,utr-3,utr-5,near-gene-3,near-gene-5\n";
	print "\t--model\toptional: ('compoundhet' if compound het model desired, otherwise don't specify this option at all)\n";
	print "\t--mafcutoff\toptional (cutoff MAF for filtering out common variants using 1000 Genomes and/or ESP; any var with freq > cutoff is excluded; default 0.01)\n";
	die;
}

