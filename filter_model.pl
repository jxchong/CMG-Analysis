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



my ($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, $inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minqual, $maxmissesperfamily);
my $debugmode;
my %genehits;
my ($countinputvariants, $printparams, $workingchr) = (0, 0, 0);

my ($count_dpexclude, $count_qualexclude, $count_notunique, $count_examined_variants, $count_variants_excluded_function, $count_variants_excluded_gatk) = ((0) x 6);
my $countvariantsmatchmodel = 0;
my $counterrorvariants = 0;
my $countcommonvariants = 0;
my $countexcludedvariants = 0;

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'subjectreq=s' => \$subjectdeffile,
	'minhits=i' => \$minhits,
	'maxmissesperfamily=i' => \$maxmissesperfamily,
	'GATKkeep=s' => \$filters,
	'N=s' => \$isNhit,
	'excludefunction=s' => \$excludeGVSfunction,
	'model:s' => \$inheritmodel,
	'mafcutoff:f' => \$mafcutoff,
	'errorcutoff:f' => \$cmgfreqcutoff,
	'dp:f' => \$mindp,
	'gq:f' => \$minqual,
	'debug:i' => \$debugmode,
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
} elsif (!defined $maxmissesperfamily) {
	optionUsage("option --maxmissesperfamily not defined\n");
}
if (!defined $mindp) {
	$mindp = 10;
}
if (!defined $minqual) {
	$minqual = 30;
}
if (!defined $mafcutoff) {
	$mafcutoff = 0.01;
}
if (!defined $cmgfreqcutoff) {
	$cmgfreqcutoff = 0.2;
} 
if (!defined $inheritmodel) {
	$inheritmodel = 'NA';
}

if (!defined $excludeGVSfunction) {
	optionUsage("option --excludefunction not defined\n");
}
if (!defined $debugmode) {
	$debugmode = 0;
} else {
	print STDOUT "Debug mode (#$debugmode) is on\n";
}


my %allowedGATKfilters;
if ($filters eq 'all' || $filters eq 'any') {
	%allowedGATKfilters = map {$_ => 1} qw(QUALFilter QDFilter LowQual SBFilter PASS ABFilter LowQual HRunFilter SnpCluster QUALFilter);
} elsif ($filters eq 'default') {
	%allowedGATKfilters = map {$_ => 1} qw(PASS);
} else {
	%allowedGATKfilters = map {$_ => 1} split(',', $filters);
} 

my %GVStoexclude;
if ($excludeGVSfunction eq 'default') {
	%GVStoexclude = map {$_ => 1} qw(intron intergenic coding-synonymous utr-3 utr-5 near-gene-3 near-gene-5 none);
} elsif ($excludeGVSfunction eq 'none')  {
	%GVStoexclude = map {$_ => 1} qw(NA);
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
########## Order of subjects in this file must correspond to order of genotype columns
while (<SUBJECTS>) {
	$_ =~ s/\s+$//;					# Remove line endings
	print LOG "$_\n";
	my ($familyid, $subjectid, $father, $mother, $sex, $relation, $desiredgeno) = split("\t", $_);
	$subjects{$subjectid} = [$familyid, $father, $mother, $sex, $relation, $desiredgeno];
	push(@orderedsubjects, $subjectid);
	if ($subjectid !~ '#') {
		push(@{$countuniquefamilies_hash{$familyid}}, $relation);
	} else {
		print LOG "$subjectid is being skipped in this analysis\n";
	}
}
close SUBJECTS;
my $countuniquefamilies = scalar(keys %countuniquefamilies_hash);
print LOG "\n";



############## BEGIN READING INPUT and WRITING LOG FILE ############## 


open (OUT, ">$outputfile") or die "Cannot write to $outputfile: $!.\n";
if ($inputfile =~ /\.gz$/) {
	open (FILE, "zcat $inputfile |") or die "Cannot read $inputfile file: $!.\n";
} else {
	open (FILE, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
}
print STDOUT "Reading in genotypes from $inputfile\n";


###################### BEGIN parse the header #############################
my ($polyphencol, $phastconscol, $gerpcol, $gatkfiltercol, $freqinCMGcol, $freqinOutsidecol, @genotypecolumns, @qualcolumns, @dpcolumns, @keepcolumns);
my $filehasfreqs = 0;
my ($vcfGQcol, $vcfDPcol, $vcfGTcol);
my $headerline;

if ($inputfile =~ m/\.vcf/) {						# skip all the metadata lines at the top of the file
	while (<FILE>) {
		next if ($_ !~ "#CHROM");
		$headerline = <FILE>;
		$headerline =~ s/\s+$//;					# Remove line endings
		my @header = split("\t", $headerline);
		last;
	}
} else {
	$headerline = <FILE>;
}
$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("\t", $headerline);
for (my $i=0; $i<=$#header; $i++) {
	my $columnname = $header[$i];
	if ($columnname =~ /Qual/i) {
		push(@qualcolumns, $i);
	} elsif ($columnname =~ /Depth/i && $columnname ne 'averageDepth') {
		push(@dpcolumns, $i);
	} elsif ($columnname =~ /Gtype/i) {
		$columnname =~ s/Gtype//; 
	}
	if ($columnname !~ /Qual/i && $columnname !~ /Depth/i) {
		push(@keepcolumns, $i);
	}
	if (defined $subjects{$columnname} || defined $subjects{"#$columnname"}) {
		push(@genotypecolumns, $i);
	}
	if ($columnname =~ /Freqin/i) {
		$filehasfreqs = 1;
		if ($columnname =~ /PrctFreqinCMG/i) {
			$freqinCMGcol = $i;
		} elsif ($columnname =~ /PrctFreqinOutsidePop/i) {
			$freqinOutsidecol = $i;
		}		
	} elsif ($columnname =~ /polyPhen/i) {
		$polyphencol = $i;
	} elsif ($columnname =~ /GERP/i) {
		$gerpcol = $i;
	} elsif ($columnname =~ /PhastCons/i) {
		$phastconscol = $i;
	} elsif ($columnname =~ /filterFlagGATK/i) {
		$gatkfiltercol = $i;
	} elsif ($columnname =~ /FILTER/i) {
		$gatkfiltercol = $i;
	}
}

###################### END parse the header #############################

if ($filehasfreqs == 1) {
	print LOG "File contains frequencies in outside populations and frequency observed in other CMG subjects\n";
} else {
	print STDERR "!!!!!! File doesn't contain frequencies in outside populations and frequency observed in other CMG subjects\n";
	print LOG "!!!!!! File doesn't contain frequencies in outside populations and frequency observed in other CMG subjects\n";
	print LOG "Run: perl ~/bin/add_error_populationfreqs.pl --in $inputfile --out <outputfile> --capture <capture array for systematic error filtering>\n";
	die;
}

# print OUT "$headerline\tFamilieswHits\n";
print OUT join("\t", @header[@keepcolumns])."\tFamilieswHits\n";

while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	$countinputvariants++;
	
	if ($printparams == 0) {
		print LOG "Requiring hits in gene in at least $minhits subjects/families\n";
		print LOG "Allowing up to $maxmissesperfamily subjects per family to not have the correct genotype\n";
		if (@dpcolumns && @qualcolumns) {
			print LOG "Genotypes require at least $mindp depth and at least $minqual qual, otherwise genotype set to be missing\n";
		} else {
			print LOG "Per-sample depth and quality information not available\n";
		}
		print LOG "Missing genotypes/no calls are counted as: $isNhit\n";
		print LOG "Excluding all variants with annotations: ".join(" ", keys %GVStoexclude)."\n";
		print LOG "Only allow variants with GATK filter: ".join(" ", keys %allowedGATKfilters)."\n";
		print LOG "Excluding variants with MAF>$mafcutoff in ESP and/or 1000 Genomes\n";
		print LOG "Excluding variants with frequency>=$cmgfreqcutoff in CMG subjects (likely systematic error)\n";
		print LOG "Special analysis? $inheritmodel\n";
		$printparams = 1;
	}
	
	my @line = split ("\t", $_);
	my $countmatches = 0;
	my ($chr, $pos, $vartype, $ref, $alt);
	my (@subjectgenotypes, @subjectquals, @subjectdps);
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
		@subjectgenotypes = @line[@genotypecolumns];
		if (@dpcolumns) {
			@subjectdps = @line[@dpcolumns];
		}
		if (@qualcolumns) {
			@subjectquals = @line[@qualcolumns];			
		}
		$functionimpact = $line[10];
		$phastcons = $line[$phastconscol];
		$gerp = $line[$gerpcol];
		$polyphen = $line[$polyphencol];
		# $inUWexomes = $line[16];
		# $UWexomescovered = $line[17];
	} 
	# elsif ($inputfile =~ m/\.vcf/) {
	# 	($chr, $pos, $ref, $alt) = ($line[0], $line[1], $line[3], $line[4]);
	# 	if ($inputfile =~ m/snps/i) {
	# 		$vartype = 'SNP';
	# 	} elsif ($inputfile =~ m/indels/i) {
	# 		$vartype = 'indel';
	# 	}
	# 	$filterset = $line[$gatkfiltercol];
	# 	@subjectgenotypes = @line[@genotypecolumns];
	# 	if (@dpcolumns) {
	# 		@subjectdps = @line[@dpcolumns];
	# 	}
	# 	if (@qualcolumns) {
	# 		@subjectquals = @line[@qualcolumns];			
	# 	}
		
		# $gene = $line[7];
		# $functionimpact = $line[10];
		# $phastcons = $line[$phastconscol];
		# $gerp = $line[$gerpcol];
		# $polyphen = $line[$polyphencol];
	# } 
	else {
		die "Input file ($inputfile) isn't an SSAnnotation or SeattleSeqAnnotation134 file\n";
	} 
	
	# Check that user input corresponds to input data files
	if (scalar(@subjectgenotypes) != scalar(@orderedsubjects)) {
		die "Your input subject definition file ($subjectdeffile) lists a different number of subjects (".scalar(@orderedsubjects).") than are contained (".scalar(@subjectgenotypes).") in the input data file ($inputfile)\n";
	}
	# Print current chromosome being processed (for user knowledge)
	if ($workingchr ne $chr) {
		print STDOUT "Reading variants on chr $chr\n";
		$workingchr = $chr;
	}
	
	if ($debugmode >= 2) { print STDOUT "looking at $chr:$pos $vartype $ref/$alt\n"; } 			## DEBUG
	
	my ($iserror, $iscommon) = (0,0);
	my $freqinCMG = $line[$freqinCMGcol]/100;							# storing allele freq as percentage
	my $freqinOutside = $line[$freqinOutsidecol]/100;					# storing allele freq as percentage
	if ($freqinCMG >= $cmgfreqcutoff) {
		$iserror = 1;
		$counterrorvariants++;
		if ($debugmode >= 2) { print STDOUT "rejecting b/c freqinCMG $freqinCMG >= $cmgfreqcutoff\n"; } 			## DEBUG
	}
	if ($freqinOutside > $mafcutoff) {
		$iscommon = 1;
		$countcommonvariants++;
		if ($debugmode >= 2) { print STDOUT "rejecting b/c freqinOutside $freqinOutside >= $mafcutoff\n"; }			## DEBUG
	}				
	if (shouldfunctionFilter(\%GVStoexclude, $functionimpact) == 1) {
		$count_variants_excluded_function++;
		if ($debugmode >= 2) { print STDOUT "rejecting b/c of function filter\n"; }			## DEBUG
	}
	if (passGATKfilters(\%allowedGATKfilters, $filterset) == 0) {
		$count_variants_excluded_gatk++;
		if ($debugmode >= 2) { print STDOUT "rejecting b/c of GATK filter\n"; }			## DEBUG
	}

	if ($iserror==0 && $iscommon==0 && shouldfunctionFilter(\%GVStoexclude, $functionimpact)==0 && ($filterset eq 'NA' || passGATKfilters(\%allowedGATKfilters, $filterset))) {
		if ($debugmode >= 2) { print STDOUT "variant at $pos passes initial filters\n"; }			## DEBUG
		$count_examined_variants++;
		my %checkfamilies;
		my %qualityflags;
		my $countcarriers = 0;
		my $ishet_excluded = 0;
		for (my $i=0; $i<=$#subjectgenotypes; $i++) {										# determine genotype for each subject
			my ($qual, $dp) = ('NA', 'NA');
			if (@dpcolumns) {
				$dp = $subjectdps[$i];
			}
			if (@qualcolumns) {
				$qual = $subjectquals[$i];			
			}
			my $genotype = $subjectgenotypes[$i];
			my $subjectid = $orderedsubjects[$i];
			if ($subjectid !~ '#') {
				my ($familyid, $father, $mother, $sex, $relation, $desiredgeno) = @{$subjects{$subjectid}};
				if (!exists $checkfamilies{$familyid}) {									# initialize family data
					my @familymembers = @{$countuniquefamilies_hash{$familyid}};
					foreach my $member (@familymembers) {
						$checkfamilies{$familyid}{$member} = -1;
						$qualityflags{$familyid}{$member} = -1;
					}
				}
				
				my $thissubjflag = "";
				if ($dp ne 'NA' && $qual ne 'NA') {															# if depth or qual information available
					if ($dp < $mindp) {
						$thissubjflag .= "DP";
					}
					if ($qual < $minqual) {
						$thissubjflag .= "GQ";
					}
					# for this variant, determine if this subject has at least one copy of the alt allele, if so, count as a carrier (for determining if this is a unique de novo)
					if ($dp >= $mindp || $qual >= $minqual) {	
						if (checkGenoMatch($vartype, $ref, $alt, $genotype, 'alt', $isNhit, $sex, $chr) || checkGenoMatch($vartype, $ref, $alt, $genotype, 'het', $isNhit, $sex, $chr)) {
							$countcarriers++;
						}
					}
				} elsif (checkGenoMatch($vartype, $ref, $alt, $genotype, 'alt', $isNhit, $sex, $chr) || checkGenoMatch($vartype, $ref, $alt, $genotype, 'het', $isNhit, $sex, $chr)) {
					# for this variant, determine if this subject has at least one copy of the alt allele, if so, count as a carrier (for determining if this is a unique de novo)
					$countcarriers++;
				}
				

				my $ismatch = checkGenoMatch($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr);	
				$checkfamilies{$familyid}{$relation} = $ismatch;
				$qualityflags{$familyid}{$relation} = $thissubjflag;
				if ($debugmode >= 3) { my $matchtext = 'is'; if ($ismatch==0) {$matchtext='is not';} print STDOUT "$familyid-$relation $subjectid genotype $matchtext a match to $desiredgeno: with GQ/DP flag=($thissubjflag)\n"; }			## DEBUG
			} 
		}
		
		my @countfamiliesrejectqual = (0,0);
		my $countfamiliesmatchmodel = 0;
		# my $countfamiliesmatch = 0;
		my %matchingfamilies;
		my %matchingfamilyunits;
		while (my ($familyid, $thisfamily_ref) = each %checkfamilies) {						# for each family, check if genotype for each family member matches model
			my %thisfamily = %{$thisfamily_ref};
			my $thisfamilymatch = 0;
			my $familysize = grep { $thisfamily{$_} != -1 } keys %thisfamily;	
			my @rejectquality = (0,0);

			if ($inheritmodel =~ 'compoundhet' && $familysize >= 3) {							
				if ( ($thisfamily{'father'}+$thisfamily{'mother'} == 1) || ($inheritmodel eq 'compoundhetmosaic' && $thisfamily{'father'}+$thisfamily{'mother'} == 2) ) {										# only one parent can be het
					if ($debugmode >= 3) { print STDOUT "family=$familyid parent matches = $thisfamily{'father'} + $thisfamily{'mother'}\n"; }
					my @familymembers = @{$countuniquefamilies_hash{$familyid}};
					foreach my $member (@familymembers) {
						if ($member ne 'mother' && $member ne 'father' && $thisfamily{$member} == 1) {
							$thisfamilymatch += 1;
							if ($qualityflags{$familyid}{$member} =~ m/DP/i) {
								$rejectquality[0] = 1;															# store whether someone failed the DP filter
							} elsif ($qualityflags{$familyid}{$member} =~ m/GQ/i) {
								$rejectquality[1] = 1;															# store whether someone failed the GQ filter
							}
						}
					}
					my $matchsubj;
					if ($thisfamily{'father'} == 1) {
						$matchsubj = 'father';
					}
					if ($thisfamily{'mother'} == 1) {
						$matchsubj = 'mother';
					}
					if ($thisfamilymatch == ($familysize-2) || ($familysize-2-$thisfamilymatch <= $maxmissesperfamily)) {					# if all affected children have at least one hit
						$countfamiliesmatchmodel++;
						if (($rejectquality[0]+$rejectquality[1]) == 0) {
							if ($debugmode >= 3) { print STDOUT "family=$familyid has a hit at $gene:$pos, the variant comes from $matchsubj and all members pass the GQ/DP check\n"; }
							$matchingfamilyunits{$familyid}{'source'} = $matchsubj;
							foreach my $member (@familymembers) {
								$matchingfamilyunits{$familyid}{$member} = $thisfamily{$member};
							}
						} else {
							$countfamiliesrejectqual[0] += $rejectquality[0];
							$countfamiliesrejectqual[1] += $rejectquality[1];
							if ($debugmode >= 3) { print STDOUT "family=$familyid has a hit in $gene:$pos, the variant comes from $matchsubj but all members do NOT pass the GQ/DP check\n"; }
						}
					}
				} 
			} else {
				while (my ($relation, $thismatch) = each %thisfamily) {
					if ($thismatch != -1) {
						$thisfamilymatch += $thismatch;
						if ($qualityflags{$familyid}{$relation} =~ m/DP/i) {
							$rejectquality[0] = 1;
						} elsif ($qualityflags{$familyid}{$relation} =~ m/GQ/i) {
							$rejectquality[1] = 1;
						}
					}
				}
				
				if ($debugmode >= 3) { print STDOUT "family=$familyid has $thisfamilymatch individuals out of $familysize matching the desired genotype, allowing only $maxmissesperfamily\n"; }
				
				if ($thisfamilymatch == $familysize || ($familysize>=2 && $maxmissesperfamily>0 && ($familysize-$thisfamilymatch) <= $maxmissesperfamily)) {
					$countfamiliesmatchmodel++;
					if (($rejectquality[0]+$rejectquality[1]) == 0) {
						$countfamiliesmatchmodel += 1;
						$matchingfamilies{$familyid} = 1;
					} else {
						$countfamiliesrejectqual[0] += $rejectquality[0];
						$countfamiliesrejectqual[1] += $rejectquality[1];
						if ($debugmode >= 2) { print STDOUT "rejecting b/c of GQ/DP\n"; }			## DEBUG
					}
				}
			}
		}


		if ($countfamiliesmatchmodel > 0) {
			$countvariantsmatchmodel++;
		}
		$count_dpexclude += $countfamiliesrejectqual[0];
		$count_qualexclude += $countfamiliesrejectqual[1];

		if ($debugmode >= 3) { print STDOUT scalar(keys %matchingfamilyunits)." family units has a hit at $gene:$pos\n"; }			## DEBUG
		if ($countfamiliesmatchmodel > 0) {
			my $data = join("\t", @line);
			if ($inheritmodel eq 'unique') {
				if ($countcarriers <= 1) {																				# if 0(only the original subject if DP/GQ not available) or 1 carriers
					push(@{$genehits{$gene}}, [$data, \%matchingfamilies, \%matchingfamilyunits]);						# under a unique de novo model, only store this as a hit if a single individual in a single family has the mutation
				} else {
					$count_notunique++;
					if ($debugmode >= 2) { print STDOUT "rejecting b/c of not unique\n"; }			## DEBUG
				}
			} else {
				push(@{$genehits{$gene}}, [$data, \%matchingfamilies, \%matchingfamilyunits]);							# store this as a hit
			}
		}
	} else {
		$countexcludedvariants++;
	}
	
	if ($debugmode >= 2) { print STDOUT "\n"; }			## DEBUG
	
	# if ($chr == 12) {
	# 	print STDOUT "stopping at chromosome $chr\n";
	# 	last;
	# }
}
close FILE;



# Check putative hit list for all genes and count number of valid hits across all families
print "\nChecking all genes for desired number of hits in desired number of families\n";
my $countgeneswhits = 0;
my $countgenesrejectedhits = 0;
my $countoutputvariants = 0;
my $genecounter = 0;
my $totalgenes = scalar(keys %genehits);
while (my ($gene, $results_ref) = each %genehits) {
	$genecounter++;
	if ($debugmode >= 1) { print "\n===========================================================================\n"; }
	if ($genecounter % 1000 == 0) {
		print "Checking gene #$genecounter (out of $totalgenes) for hits\n";
	}
	my @hitdata = @{$results_ref};
	
	# account for possible compound het model, then count hits in each family for this gene
	if ($debugmode >= 1) { print "in $gene:\n"; }		
	my $resultsFamiliesvsModel_ref = checkFamiliesvsModel(\@hitdata, $inheritmodel);
	# review all families for required number of hits in this gene
	my $enoughfamilieshavehits = checkFamiliesforHits($resultsFamiliesvsModel_ref, $inheritmodel, $countuniquefamilies, $minhits);													
	
	if ($enoughfamilieshavehits == 1) {
		$countgeneswhits++;
	  	foreach my $hit (@hitdata) {
	  		my $hitvarinfo = ${$hit}[0];
			my %matchingfamilies = %{${$hit}[1]};
			my %matchingfamilyunits = %{${$hit}[2]};
			if ($debugmode >= 1) { print "Enough families, total=".scalar(keys %matchingfamilies)."+".scalar(keys %matchingfamilyunits)." have hits\n"; }
			
			my @familyids = ((keys %matchingfamilies), (keys %matchingfamilyunits));
	  		$countoutputvariants++;
			my @hitvarinfoarray = split("\t", $hitvarinfo);
	  		# print OUT "$hitvarinfo\t".join(";", @familyids)."\n";
			print OUT join("\t", @hitvarinfoarray[@keepcolumns])."\t".join(";", @familyids)."\n";
	  	}
	} else {
		$countgenesrejectedhits++;
		if ($debugmode >= 1) { print "Rejecting $gene b/c not enough families/subjects within family have hits\n"; }
	}
}

print LOG "\nResults summary:\n";
print LOG "In $countgeneswhits gene(s), identified $countoutputvariants variants\n";
print LOG "Total $countinputvariants variants in input\n";
print LOG "Total $count_examined_variants variants from input examined after excluding $countexcludedvariants:\n";
print LOG "\tN=$counterrorvariants are systematic errors\n";
print LOG "\tN=$countcommonvariants are common variants\n";
print LOG "\tN=$count_variants_excluded_gatk based on GATK filter\n";
print LOG "\tN=$count_variants_excluded_function based on functional annotation\n";


print LOG "Filtered out candidate hits from list of $countvariantsmatchmodel variants matching inheritance model:\n";
print LOG "\tN=$count_dpexclude variants excluded based on DP; N=$count_qualexclude variants based on GQ\n";
print LOG "\tN=$count_notunique excluded because they were not unique within this dataset\n";
close LOG;
 


print STDOUT "In $countgeneswhits gene(s), matched $countoutputvariants variants\n";


sub checkFamiliesvsModel {																										# sum up hits in each family
	# account for possible compound het model (count number of hits per affected child in family)
	my ($inputhitdata_ref, $inheritmodel) = @_;
	my %familieswHitsinGene;
	foreach my $hit (@{$inputhitdata_ref}) {
		my $hitvarinfo = ${$hit}[0];
			my @thishit = split("\t", $hitvarinfo);
		my %matchingfamilies = %{${$hit}[1]};
		my %matchingfamilyunits = %{${$hit}[2]};

		while (my($familyid, $ismatch) = each %matchingfamilies) {
			$familieswHitsinGene{$familyid}++;	
			if ($debugmode >= 2) { print "family=$familyid has a match\n"; }			
		}

		# total up hits per family member
		if ($inheritmodel =~ 'compoundhet') {
			if ($debugmode >= 1) { print "For this hit, implementing compoundhet model (add to total hits per individual in total ".scalar(keys %matchingfamilyunits)." family units)\n"; }	
			while (my($familyid, $thisfamilydata_ref) = each %matchingfamilyunits) {
				if ($debugmode >= 1) { print "Counting hits in individuals in family=$familyid\n"; }	
				if (!defined $familieswHitsinGene{$familyid}) {
					$familieswHitsinGene{$familyid}{'father'} = 0;
					$familieswHitsinGene{$familyid}{'mother'} = 0;
				}
				my $matchsubj = ${$thisfamilydata_ref}{'source'};
				$familieswHitsinGene{$familyid}{$matchsubj}++;											# hit comes from xxx parent
				if ($debugmode >= 1) { print "a hit in family=$familyid comes from $matchsubj\n"; }	
				
				my @familymembers = @{$countuniquefamilies_hash{$familyid}};
				foreach my $member (@familymembers) {
					if ($member ne 'mother' && $member ne 'father') {							# FIX!!!! (need to accomodate allowing for some kids with no hits)
						$familieswHitsinGene{$familyid}{$member} += ${$thisfamilydata_ref}{$member};								# add one to hits in children
						if ($debugmode >= 2) { print "family=$familyid $member added ${$thisfamilydata_ref}{$member} hits making a total $familieswHitsinGene{$familyid}{$member} hits\n"; }			
					}
				}				
			}
		}
	}
	return (\%familieswHitsinGene);
}


sub checkFamiliesforHits {																				# make sure required number of hits in each family/individual
	my($resultsFamiliesvsModel_ref, $inheritmodel, $countuniquefamilies, $minhits) = @_;
	my %familieswHitsinGene = %{$resultsFamiliesvsModel_ref};
	
	my $countfamiliesmatch = 0;
	my @trackfamilyidswHits;
	while (my($familyid, $familyhits) = each %familieswHitsinGene) {																		
		if ($inheritmodel =~ 'compoundhet') {
			if (ref($familyhits) eq 'HASH') {
				if ($debugmode >= 1) { print "Counting hits in each individual in $familyid\n"; }	
				my $counthitsinfamily = 0;
				my @familymembers = keys %{$familyhits};
				foreach my $member (@familymembers) {
					if ($debugmode >= 1) { print "${$familyhits}{$member} hit(s) in $member\n"; }	
					if ($member ne 'mother' && $member ne 'father' && ${$familyhits}{$member} >= 2) {
						$counthitsinfamily++;
					}
				}

				if ($debugmode >= 1) { print "$counthitsinfamily hits in family vs ".scalar(@familymembers)." family members\n"; }	
				
				if (${$familyhits}{'mother'} >= 1 && ${$familyhits}{'father'} >= 1) {
					# $counthitsinfamily += 2;
					if ($counthitsinfamily+$maxmissesperfamily+2 >= scalar(@familymembers)) {
						$countfamiliesmatch++;
					}
				}
			} elsif (!ref($familyhits)) {																									
				if ($familyhits >= 2) {
					$countfamiliesmatch++;
				}
			}
		} else {
			if ($debugmode >= 1) { print "family=$familyid has hits? $familyhits\n"; }
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

sub passGATKfilters {
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
	my ($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr) = @_;
		
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
			if ($sex == 1 && ($chr eq "X" || $chr eq "Y") && ($alleles[0] eq $alt && $alleles[1] eq $alt)) {
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
	print "\t--maxmissesperfamily\tmax number of individuals with genotypes not matching model in each family unit\n";
	print "\t--GATKkeep\tGATK quality filters that should be kept (comma-delimited with no spaces, or keep 'all' or 'default')\n";
	print "\t\tdefault=PASS\n";
	print "\t--N\thit or nothit (how should we count missing genotypes)\n";
	print "\t--excludefunction\tcomma separated list of GVS variant function classes to be excluded, or 'default'\n";
		print "\t\tdefault=intron,intergenic,coding-synonymous,utr-3,utr-5,near-gene-3,near-gene-5\n";
	print "\t--model\toptional: ('compoundhet' if compound het model desired, 'unique' if filtering for de novo unique; otherwise don't specify this option at all)\n";
	print "\t--mafcutoff\toptional: (cutoff MAF for filtering out common variants using 1000 Genomes and/or ESP; any var with freq > cutoff is excluded; default 0.01)\n";
	print "\t--errorcutoff\toptional: (cutoff frequency for filtering out systematic errors based on frequency in all CMG exomes to date)\n";
	die;
}

