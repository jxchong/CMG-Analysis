#!/usr/bin/env perl
#
# Description: Look for mutations matching a particular model given an SSAnnotation file produced by SeattleSeq.
#	Note compound het mosaic model doesn't work with parents lacking data
#
#
# Created by Jessica Chong on 2012-10-15.


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Moose;

my $help = 0;
my ($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, $inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minGQ, $maxmismatchesperfamily, $debugmode, $logfile);
my ($allowedGATKfilters_ref, $GVStoexclude_ref);
GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'subjectreq=s' => \$subjectdeffile,
	'minhits=i' => \$minhits,
	'maxmissesperfamily=i' => \$maxmismatchesperfamily,
	'GATKkeep=s' => \$filters,
	'N=s' => \$isNhit,
	'excludefunction=s' => \$excludeGVSfunction,
	'model:s' => \$inheritmodel,
	'mafcutoff:f' => \$mafcutoff,
	'errorcutoff:f' => \$cmgfreqcutoff,
	'dp:f' => \$mindp,
	'gq:f' => \$minGQ,
	'debug:i' => \$debugmode,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;

($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, 
	$inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minGQ, $maxmismatchesperfamily, 
	$debugmode, $allowedGATKfilters_ref, $GVStoexclude_ref, $logfile) = checkandstoreOptions($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, $inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minGQ, $maxmismatchesperfamily, $debugmode);
my %allowedGATKfilters = %{$allowedGATKfilters_ref};
my %GVStoexclude = %{$GVStoexclude_ref};

open (my $log_filehandle, ">", $logfile) or die "Cannot write to $logfile: $!.\n";
print $log_filehandle "input=$inputfile\noutput=$outputfile\nsubjectreq=$subjectdeffile\n";
printParamstoLog($log_filehandle, $minhits, $maxmismatchesperfamily, $mindp, $minGQ, $isNhit, join(" ", keys %GVStoexclude), join(" ", keys %allowedGATKfilters), $mafcutoff, $cmgfreqcutoff, $inheritmodel);

###################### READ PEDIGREE/SUBJECT DATA ###################### 
my ($countuniquefamilies_hash_ref, $orderedsubjects_ref, $subjects_ref, $countuniquefamilies) = readPedigree($subjectdeffile, $log_filehandle);
my %countuniquefamilies_hash = %{$countuniquefamilies_hash_ref};
my @orderedsubjects = @{$orderedsubjects_ref};
my %subjects = %{$subjects_ref};
	
###################### BEGIN READING INPUT and HEADER ###################### 
open (my $output_filehandle, ">", $outputfile) or die "Cannot write to $outputfile: $!.\n";
my $input_filehandle;
if ($inputfile =~ /\.gz$/) {
	open ($input_filehandle, "zcat $inputfile |") or die "Cannot read $inputfile file: $!.\n";
} else {
	open ($input_filehandle, "$inputfile") or die "Cannot read $inputfile file: $!.\n";
}
print STDOUT "Reading in genotypes from $inputfile\n";

my $countskipheaderlines = 0;
my ($vcfGQcol, $vcfDPcol, $vcfGTcol);
my $headerline;
if ($inputfile =~ m/\.vcf/) {						# skip all the metadata lines at the top of the file
	while (<$input_filehandle>) {
		$headerline = $_;
		$countskipheaderlines++;
		next if ($headerline !~ "#CHROM");
		$headerline =~ s/\s+$//;					# Remove line endings
		my @header = split("	", $headerline);
		last;
	}
} else {
	$headerline = <$input_filehandle>;
}

if ($debugmode >= 1) {
	print STDOUT "$countskipheaderlines lines skipped in header\n";
}

$headerline =~ s/\s+$//;					# Remove line endings
my @header = split("	", $headerline);
my $inputfiletype = determineInputType($headerline, $inputfile);					# SSAnnotation, vcf, or SeattleSeqAnnotation
print "Input file type is $inputfiletype\n";
my ($polyphencol, $phastconscol, $gerpcol, $gatkfiltercol, $freqinCMGcol, $freqinOutsidecol, $genotypecolumns_ref, $qualcolumns_ref, $dpcolumns_ref, $keepcolumns_ref) = parseHeader(@header);
my @genotypecolumns = @{$genotypecolumns_ref};
my @qualcolumns = @{$qualcolumns_ref};
my @dpcolumns = @{$dpcolumns_ref};
my @keepcolumns = @{$keepcolumns_ref};

# Check that user input corresponds to input data files
if (scalar(@genotypecolumns) != scalar(@orderedsubjects)) {
	die "Your input subject definition file ($subjectdeffile) lists a different number of subjects (".scalar(@orderedsubjects).") than are contained (".scalar(@genotypecolumns).") in the input data file ($inputfile)\n";
}

###################### PRINT HEADER FOR OUTPUT #############################
if ($inputfiletype ne 'vcf') {
	print $output_filehandle "#".join("	", @header[@keepcolumns])."	FamilieswHits\n";
} else {
	print $output_filehandle "#chr	start	end	type	ref	alt	GATKflag	geneList	rsID	functionGVS	aminoacids	proteinPos	cDNAPos	PhastCons	GERP	Polyphen	AltFreqinCMG	AltFreqOutside	clinAssoc	KEGG	";
	print $output_filehandle join("Gtype	", @orderedsubjects)."Gtype	";
	print $output_filehandle join("DP	", @orderedsubjects)."DP	";
	print $output_filehandle join("GQ	", @orderedsubjects)."GQ	FamilieswHits\n";
}


###################### start checking variants in input #############################
my ($countinputvariants, $printparams, $workingchr, $count_dpexclude, $count_qualexclude, $count_notunique, $count_examined_variants, $count_variants_excluded_function, $count_variants_excluded_gatk) = ((0) x 9);
my ($countvariantsmatchmodel, $counterrorvariants, $countcommonvariants, $countexcludedvariants) = ((0) x 4);
my %genehits;
while ( <$input_filehandle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	$countinputvariants++;
	
	my @line = split ("	", $_);
	my $countmatches = 0;
	my ($chr, $pos, $vartype, $ref, $alt, $filterset);
	my ($subjgeno_ref, $subjdps_ref, $subjquals_ref, @subjectgenotypes, @subjectquals, @subjectdps);
	my ($gene, $functionimpact, $gerp, $polyphen, $phastcons, $rsid, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside, $clinical, $kegg);
	
	if ($debugmode >= 3) {
		print STDOUT "\n";
	}
	
	if ($inputfiletype eq 'SeattleSeqAnnotation') {
		# parse_SeattleSeqAnnotation134_byline();
	} elsif ($inputfiletype eq 'SSAnnotation') {
		($chr, $pos, $vartype, $ref, $alt, $filterset, $gene, $subjgeno_ref, $subjdps_ref, $subjquals_ref, $functionimpact) = parse_SSAnnotation134_byline(\@genotypecolumns, \@dpcolumns, \@qualcolumns, @line);
		$freqinCMG = $line[$freqinCMGcol];		
		$freqinOutside = $line[$freqinOutsidecol];
		if ($debugmode >= 3) {
			print STDOUT "Looking at: $chr, $pos, $vartype, $ref, $alt, $filterset, $gene, $functionimpact, $freqinCMG, $freqinOutside\n";
		}
	} elsif ($inputfiletype eq 'vcf') {
		($chr, $pos, $vartype, $ref, $alt, $filterset, $gene, $subjgeno_ref, $subjdps_ref, $subjquals_ref, $functionimpact, $gerp, $polyphen, $phastcons, $rsid, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside, $clinical, $kegg) = parse_vcf_byline(\@genotypecolumns, @line);
		# if ($debugmode >= 4) {
		# 	if ($chr < 2) {
		# 		next;
		# 	} elsif ($chr > 2) {
		# 		last;
		# 	}			
		# }
		if ($debugmode >= 3) {
			print STDOUT "Looking at: $chr, $pos, $vartype, $ref, $alt, $filterset, $gene, $functionimpact, $gerp, $polyphen, $phastcons, $rsid, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside\n";
		}
	} else {
		die "Input file ($inputfile) isn't an SSAnnotation or SeattleSeqAnnotation134 file\n";
	} 
	@subjectgenotypes = @{$subjgeno_ref};
	@subjectdps = @{$subjdps_ref};
	@subjectquals = @{$subjquals_ref};
	
	if ($debugmode >= 4) {
		print STDOUT "genotypes=@subjectgenotypes\n";
		print STDOUT "dps=@subjectdps\n";
		print STDOUT "GQs=@subjectquals\n";
	}
	
	if ($workingchr ne $chr && ($chr eq '1' || $chr eq '6' || $chr eq '12' || $chr eq '22')) {
		print STDOUT "Reading variants on chr $chr\n";			# Print current chromosome being processed (for user knowledge)
		$workingchr = $chr;
	}
			
	my ($iserror, $iscommon) = (0,0);
	$freqinCMG = $freqinCMG/100;							# store allele freq as percentage, so convert back to frequency
	$freqinOutside = $freqinOutside/100;					# store allele freq as percentage, so convert back to frequency
	if ($freqinCMG >= $cmgfreqcutoff) {
		$iserror = 1;
		$counterrorvariants++;
		if ($debugmode >= 4) { print STDOUT "rejecting b/c freqinCMG $freqinCMG >= $cmgfreqcutoff\n"; } 			## DEBUG
	}
	if ($freqinOutside > $mafcutoff) {
		$iscommon = 1;
		$countcommonvariants++;
		if ($debugmode >= 4) { print STDOUT "rejecting b/c freqinOutside $freqinOutside >= $mafcutoff\n"; }			## DEBUG
	}
	my $resultfunctionfilter = shouldfunctionFilter(\%GVStoexclude, $functionimpact);
	my $resultGATKfilter = passGATKfilters(\%allowedGATKfilters, $filterset);
	if ($resultfunctionfilter == 1) {
		$count_variants_excluded_function++;
		if ($debugmode >= 4) { print STDOUT "rejecting b/c of function filter\n"; }			## DEBUG
	}
	if ($resultGATKfilter == 0) {
		$count_variants_excluded_gatk++;
		if ($debugmode >= 4) { print STDOUT "rejecting b/c of GATK filter\n"; }			## DEBUG
	}

	if ($iserror==0 && $iscommon==0 && $resultfunctionfilter==0 && ($filterset eq 'NA' || $resultGATKfilter==1)) {
		if ($debugmode >= 2) { print STDOUT "variant at $chr:$pos in $gene passes initial filters\n"; }			## DEBUG
		$count_examined_variants++;
		my %checkfamilies;
		my %typeofgeno_all;
		my %compoundhetcarriers;
		my %qualityflags;
		my $countcarriers = 0;
		my $ishet_excluded = 0;
		for (my $i=0; $i<=$#subjectgenotypes; $i++) {										# determine genotype for each subject
			my ($qual, $dp) = ('NA', 'NA');
			if (@subjectdps) {
				$dp = $subjectdps[$i];
			}
			if (@subjectquals) {
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
					if ($qual < $minGQ) {
						$thissubjflag .= "GQ";
					}
					# for this variant, determine if this subject has at least one copy of the alt allele, if so, count as a carrier (for determining if this is a unique de novo)
					if ($dp >= $mindp || $qual >= $minGQ) {	
						if (checkGenoMatch($vartype, $ref, $alt, $genotype, 'alt', $isNhit, $sex, $chr, $debugmode) || checkGenoMatch($vartype, $ref, $alt, $genotype, 'het', $isNhit, $sex, $chr, $debugmode)) {
							$countcarriers++;
						}
					}
				} elsif (checkGenoMatch($vartype, $ref, $alt, $genotype, 'alt', $isNhit, $sex, $chr, $debugmode) || checkGenoMatch($vartype, $ref, $alt, $genotype, 'het', $isNhit, $sex, $chr, $debugmode)) {
					# for this variant, determine if this subject has at least one copy of the alt allele, if so, count as a carrier (for determining if this is a unique de novo)
					$countcarriers++;
				}
				
				my $ismatch = checkGenoMatch($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr, $debugmode);	
				$typeofgeno_all{$familyid}{$relation} = determineGenoType($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr);
				$checkfamilies{$familyid}{$relation} = $ismatch;
				$qualityflags{$familyid}{$relation} = $thissubjflag;
				if ($inheritmodel eq 'compoundhetmosaic') {											# for compound het mosaic model, which parent should be a carrier?
					if (!defined $compoundhetcarriers{$familyid}) {
						$compoundhetcarriers{$familyid} = 'NA';
					} elsif ($desiredgeno eq 'het' && ($relation eq 'father' || $relation eq 'mother')) {
						$compoundhetcarriers{$familyid} = $relation;								
					} elsif ($desiredgeno eq 'ref' && ($relation eq 'father' || $relation eq 'mother')) {
						if ($relation eq 'father') {
							$compoundhetcarriers{$familyid} = 'mother';	
						} else {
							$compoundhetcarriers{$familyid} = 'father';	
						}
					}
				}
				if ($debugmode >= 3) { my $matchtext = 'is'; if ($ismatch==0) {$matchtext='is not';} print STDOUT "$familyid-$relation (sex=$sex) $subjectid genotype ($genotype) $matchtext a match to $desiredgeno: with GQ/DP flag=($thissubjflag) based on DP=$dp and GQ=$qual\n"; }			## DEBUG
			} 
		}
		
		my @countfamiliesrejectqual = (0,0);
		my $countfamiliesmatchmodel = 0;
		# my $countfamiliesmatch = 0;
		my %matchingfamilies;
		my %matchingfamilyunits;
		while (my ($familyid, $thisfamily_ref) = each %checkfamilies) {						# for each family, check if genotype for each family member matches model
			my %thisfamily = %{$thisfamily_ref};											# storage for whether genotype matches expectations (either 1 or 0 for match or no match)
			my ($thisfamilykidsmatch, $thisfamilyparentsmatch, $countNgenotypes) = (0,0,0);
			my $familysize = grep { $thisfamily{$_} != -1 } keys %thisfamily;	
			my $nparentswithdata = (defined $thisfamily{'father'}) + (defined $thisfamily{'mother'});
			my $nkidswithdata = $familysize-$nparentswithdata;
			if ($inheritmodel =~ 'compoundhet') {		
				if (!defined $thisfamily{'father'}) { 
					$thisfamily{'father'} = 'NA';
				}
				if (!defined $thisfamily{'mother'}) {
					$thisfamily{'mother'} = 'NA';
				}
	
				if ($debugmode >= 3) { print STDOUT "family=$familyid parent matches = $thisfamily{'father'} + $thisfamily{'mother'}\n"; }
				my @familymembers = @{$countuniquefamilies_hash{$familyid}};
				foreach my $member (@familymembers) {									
					my $countasmatch = 0;
					my $failsquality = failsCallQuality($qualityflags{$familyid}{$member});
					if ($failsquality) {
						$countNgenotypes += 1;
					}
					if ($isNhit eq 'hit' && $failsquality || ($isNhit eq 'hit' && !$failsquality && $thisfamily{$member}==1) || ($isNhit eq 'nohit' && !$failsquality && $thisfamily{$member}==1)) {
						$countasmatch = 1;
					} 
					if ($countasmatch == 1) {
						if ($member ne 'mother' && $member ne 'father') {
							$thisfamilykidsmatch += 1;																	# counts matches in the kids	
						} else {
							$thisfamilyparentsmatch += 1;																# counts matches in the parents
						}
						if ($debugmode >= 3) { print STDOUT "$member in $familyid has the correct genotype\n"; }
					}
				}
				
				my $matchsubj = determineMatchSubject($thisfamily{'father'}, $thisfamily{'mother'}, $inheritmodel, $compoundhetcarriers{$familyid}, $typeofgeno_all{$familyid}{'father'}, $typeofgeno_all{$familyid}{'mother'});
				if ($matchsubj eq 'reject' || $thisfamilykidsmatch+$thisfamilyparentsmatch==0) {
					if ($debugmode >=4) { print STDOUT "reject\n"; }
					next;							# skip to evaluation of next family; doesn't match inheritance model
				} else {
					if ($debugmode >= 3) { print STDOUT "in family=$familyid, there are $thisfamilykidsmatch+$thisfamilyparentsmatch matches out of $familysize genotypes; $thisfamilyparentsmatch/$nparentswithdata parents match; $countNgenotypes genotypes in family have low DP/GQ\n"; }
					## Note to self: under compoundhet model, if parents are available, one parent will not match at the site					
					if ( ($thisfamilyparentsmatch==0 && $thisfamilykidsmatch+$maxmismatchesperfamily>=$familysize) 
					     || ($inheritmodel =~ 'compoundhet' && $thisfamilykidsmatch+$thisfamilyparentsmatch+1+$maxmismatchesperfamily>=$familysize)
						 || ($inheritmodel eq 'compoundhetmosaic' && $thisfamilykidsmatch+$thisfamilyparentsmatch+$maxmismatchesperfamily>=$familysize) 
						 ) {		
						$countfamiliesmatchmodel++;
						if ($debugmode >= 3) { print STDOUT "family=$familyid has a hit at $gene:$pos, variant comes from $matchsubj; enough members ($thisfamilykidsmatch+$thisfamilyparentsmatch, $maxmismatchesperfamily allowed misses) have hits\n"; }
						$matchingfamilyunits{$familyid}{'source'} = $matchsubj;
						foreach my $member (@familymembers) {
							$matchingfamilyunits{$familyid}{$member} = $thisfamily{$member};
						}
					} else {
						if ($debugmode >=4) { print STDOUT "no hits, reject\n"; }
					}
				}
			} else {
				while (my ($member, $thismatch) = each %thisfamily) {
					if ($thismatch != -1) {
						my $countasmatch = 0;
						my $failsquality = failsCallQuality($qualityflags{$familyid}{$member});
						if ($failsquality) {
							$countNgenotypes += 1;
						}
						if ($isNhit eq 'hit' && $failsquality) {
							$countasmatch = 1;
						} elsif ($isNhit eq 'hit' && !$failsquality && $thismatch==1) {
							$countasmatch = 1;
						} elsif ($isNhit eq 'nohit' && !$failsquality && $thismatch==1) {
							$countasmatch = 1;
						}
						if ($countasmatch == 1) {
							if ($member ne 'mother' && $member ne 'father') {
								$thisfamilykidsmatch += 1;																	# counts matches in the kids	
							} else {
								$thisfamilyparentsmatch += 1;																# counts matches in the parents
							}
							if ($debugmode >= 3) { print STDOUT "$member in $familyid has the correct genotype\n"; }
						}
					}
				}
				if ($countNgenotypes == $familysize) {
					if ($debugmode >= 2) { print STDOUT "rejecting b/c everyone in family had missing/poor quality genotypes\n"; }			## DEBUG
					next;
				} else {
					if ($debugmode >= 3) { print STDOUT "family=$familyid has $thisfamilykidsmatch+$thisfamilyparentsmatch matching individuals out of $familysize matching the desired genotype, allowing $maxmismatchesperfamily misses; $countNgenotypes genotypes in family have low DP/GQ\n"; }
					if (($thisfamilyparentsmatch+$thisfamilykidsmatch) == $familysize || ($familysize>=2 && ($maxmismatchesperfamily+$thisfamilyparentsmatch+$thisfamilykidsmatch) >= $familysize)) {
						$countfamiliesmatchmodel++;
						$matchingfamilies{$familyid} = 1;
					} else {
						if ($debugmode >= 2) { print STDOUT "rejecting b/c not enough matches in family\n"; }			## DEBUG
					}
				}
			}
		}


		if ($countfamiliesmatchmodel > 0) {
			$countvariantsmatchmodel++;
		}
		# $count_dpexclude += $countfamiliesrejectqual[0];
		# $count_qualexclude += $countfamiliesrejectqual[1];

		if ($debugmode >= 3) { print STDOUT scalar(keys %matchingfamilyunits)." family units and ".scalar(keys %matchingfamilies)." families out of ".scalar(keys %checkfamilies)." families of some type have a hit at $gene:$pos\n"; }			## DEBUG
		if ((scalar(keys %matchingfamilyunits)+scalar(keys %matchingfamilies)) > 0) {
			my $data;
			if ($inputfiletype ne 'vcf') {
				$data = join("	", @line);
			} else {
				$data = "$chr	$pos	$pos	$vartype	$ref	$alt	$filterset	$gene	$rsid	$functionimpact	$aminoacids	$proteinpos	$cdnapos	$phastcons	$gerp	$polyphen	$freqinCMG	$freqinOutside	$clinical	$kegg	".join("	", (@subjectgenotypes, @subjectdps, @subjectquals));	
			}
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
	
	# if ($debugmode >= 2) { print STDOUT "\n"; }			## DEBUG
	
	# if ($chr == 12) {
	# 	print STDOUT "stopping at chromosome $chr\n";
	# 	last;
	# }
	# exit;												## DEBUG
}
close $input_filehandle;



# Check putative hit list for all genes and count number of valid hits across all families
print STDOUT "\nChecking all genes for desired number of hits in desired number of families\n";
my %geneswhits;
my $countgenesrejectedhits = 0;
my $countoutputvariants = 0;
my $genecounter = 0;
my $totalgenes = scalar(keys %genehits);
while (my ($gene, $results_ref) = each %genehits) {
	$genecounter++;
	if ($debugmode >= 1) { print STDOUT "\n===========================================================================\n"; }
	if ($genecounter % 1000 == 0) {
		print STDOUT "Checking gene #$genecounter (out of $totalgenes) for hits\n";
	}
	my @hitdata = @{$results_ref};
	
	# account for possible compound het model, then count hits in each family for this gene
	if ($debugmode >= 1) { print STDOUT "in $gene:\n"; }		
	my $resultsFamiliesvsModel_ref = checkFamiliesvsModel(\@hitdata, $inheritmodel);
	# review all families for required number of hits in this gene
	my $enoughfamilieshavehits = checkFamiliesforHits($resultsFamiliesvsModel_ref, $inheritmodel, $countuniquefamilies, $minhits);													
	
	if ($enoughfamilieshavehits == 1) {
		$geneswhits{$gene}++;
	  	foreach my $hit (@hitdata) {
	  		my $hitvarinfo = ${$hit}[0];
			my %matchingfamilies = %{${$hit}[1]};
			my %matchingfamilyunits = %{${$hit}[2]};
			if ($debugmode >= 1) { print STDOUT "Enough families, total=".scalar(keys %matchingfamilies)."+".scalar(keys %matchingfamilyunits)." have hits\n"; }
			
			my @familyids = ((keys %matchingfamilies), (keys %matchingfamilyunits));
	  		$countoutputvariants++;
			my @hitvarinfoarray = split("	", $hitvarinfo);
	  		# print $output_filehandle "$hitvarinfo	".join(";", @familyids)."\n";
			if ($inputfiletype ne 'vcf') {
				print $output_filehandle join("	", @hitvarinfoarray[@keepcolumns])."	".join(";", @familyids)."\n";
			} else {
				print $output_filehandle join("	", @hitvarinfoarray)."	".join(";", @familyids)."\n";
			}
	  	}
	} else {
		$countgenesrejectedhits++;
		if ($debugmode >= 1) { print STDOUT "Rejecting $gene b/c not enough families/subjects within family have hits\n"; }
	}
}

my @genes = keys %geneswhits;
my $countgeneswhits = scalar(@genes);
print $log_filehandle "\nResults summary:\n";
print $log_filehandle "In $countgeneswhits gene(s), identified $countoutputvariants variants.\n";
print $log_filehandle "Genes: ".join(", ", @genes)."\n";
print $log_filehandle "Total $countinputvariants variants in input\n";
print $log_filehandle "Total $count_examined_variants variants from input examined after excluding $countexcludedvariants:\n";
print $log_filehandle "	N=$counterrorvariants are systematic errors\n";
print $log_filehandle "	N=$countcommonvariants are common variants\n";
print $log_filehandle "	N=$count_variants_excluded_gatk based on GATK filter\n";
print $log_filehandle "	N=$count_variants_excluded_function based on functional annotation\n";

print $log_filehandle "Filtered out candidate hits from list of $countvariantsmatchmodel variants matching inheritance model:\n";
# print $log_filehandle "	N=$count_dpexclude variants excluded based on DP; N=$count_qualexclude variants based on GQ\n";
print $log_filehandle "	N=$count_notunique excluded because they were not unique within this dataset\n";
close $log_filehandle;
 


print STDOUT "In $countgeneswhits gene(s), matched $countoutputvariants variants.\n";
print STDOUT "Genes: ".join(", ", @genes)."\n";


sub checkFamiliesvsModel {																										# sum up hits in each family
	# account for possible compound het model (count number of hits per affected child in family)
	my ($inputhitdata_ref, $inheritmodel) = @_;
	my %familieswHitsinGene;
	foreach my $hit (@{$inputhitdata_ref}) {
		my $hitvarinfo = ${$hit}[0];
			my @thishit = split("	", $hitvarinfo);
		my %matchingfamilies = %{${$hit}[1]};
		my %matchingfamilyunits = %{${$hit}[2]};

		while (my($familyid, $ismatch) = each %matchingfamilies) {
			$familieswHitsinGene{$familyid}++;	
			if ($debugmode >= 2) { print STDOUT "family=$familyid has a match\n"; }			
		}

		# total up hits per family member
		if ($inheritmodel =~ 'compoundhet') {
			if ($debugmode >= 1) { print STDOUT "For this hit, implementing compoundhet model (add to total hits per individual in total ".scalar(keys %matchingfamilyunits)." family units)\n"; }	
			while (my($familyid, $thisfamilydata_ref) = each %matchingfamilyunits) {
				if ($debugmode >= 1) { print STDOUT "Counting hits in individuals in family=$familyid\n"; }	
				if (!defined $familieswHitsinGene{$familyid}) {
					$familieswHitsinGene{$familyid}{'father'} = 0;
					$familieswHitsinGene{$familyid}{'mother'} = 0;
				}
				my $matchsubj = ${$thisfamilydata_ref}{'source'};
				if ($matchsubj eq 'either') {
					$familieswHitsinGene{$familyid}{'father'}++;											# hit comes from xxx parent
					$familieswHitsinGene{$familyid}{'mother'}++;											# hit comes from xxx parent
				} else {
					$familieswHitsinGene{$familyid}{$matchsubj}++;											# hit comes from xxx parent
				}
				if ($debugmode >= 1) { print STDOUT "a hit in family=$familyid comes from $matchsubj\n"; }	
				
				my @familymembers = @{$countuniquefamilies_hash{$familyid}};
				foreach my $member (@familymembers) {
					if ($member ne 'mother' && $member ne 'father') {				
						$familieswHitsinGene{$familyid}{$member} += ${$thisfamilydata_ref}{$member};								# add one to hits in children
						if ($debugmode >= 2) { print STDOUT "family=$familyid $member added ${$thisfamilydata_ref}{$member} hits making a total $familieswHitsinGene{$familyid}{$member} hits\n"; }			
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
				if ($debugmode >= 1) { print STDOUT "Counting hits in each individual in $familyid\n"; }	
				my ($counthitsinfamily, $countparentswithdata) = 0;
				my @familymembers = keys %{$familyhits};
				foreach my $member (@familymembers) {
					if ($debugmode >= 1) { print STDOUT "${$familyhits}{$member} hit(s) in $member\n"; }	
					if ($member ne 'mother' && $member ne 'father' && ${$familyhits}{$member} >= 2) {
						$counthitsinfamily++;
					} elsif (($member eq 'mother' || $member eq 'father') && ${$familyhits}{$member} >= 1) {		# is this ok in a compoundhetmosaic model?
						$counthitsinfamily++;
					}
				}
				if ($debugmode >= 1) { print STDOUT "$counthitsinfamily hits in family vs ".scalar(@familymembers)." family members\n"; }		
				if ($debugmode >= 3) { print STDOUT "Mother has ${$familyhits}{'mother'} hits and father has ${$familyhits}{'father'} hits in family $familyid\n"; }	
				if ($counthitsinfamily+$maxmismatchesperfamily >= scalar(@familymembers)) {
					$countfamiliesmatch++;
				}
			} elsif (!ref($familyhits)) {																									
				if ($familyhits >= 2) {
					$countfamiliesmatch++;
				}
			}
		} else {
			if ($debugmode >= 1) { print STDOUT "family=$familyid has hits? $familyhits\n"; }
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


sub checkandstoreOptions {
	my ($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, 
		$inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minGQ, $maxmismatchesperfamily, 
		$debugmode) = @_;
	my (%allowedGATKfilters, %GVStoexclude);
	
	if (!defined $inputfile) {
		 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --in not defined.\n")
	} elsif (!defined $outputfile) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
	} elsif (!defined $subjectdeffile) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --subjectreq not defined\n");
	} elsif (!defined $minhits) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --minhits not defined\n");
	} elsif (!defined $filters) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --GATKkeep not defined\n");
	} elsif (!defined $isNhit) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --N not defined\n");
	} elsif (!defined $maxmismatchesperfamily) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --maxmismatchesperfamily not defined\n");
	}
	if ($isNhit ne 'hit' && $isNhit ne 'nohit') {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --N defined but not valid (should be 'hit' or 'nohit')\n");
	}
	# if ($isNhit eq 'hit') {
	# 	$mindp = 0;
	# 	$minGQ = 0;
	# } else {
		if (!defined $mindp) {
			$mindp = 8;
		}
		if (!defined $minGQ) {
			$minGQ = 20;
		}
	# }

	if (!defined $mafcutoff) {
		$mafcutoff = 0.01;
	}
	if (!defined $cmgfreqcutoff) {
		$cmgfreqcutoff = 0.2;
	} 
	if (!defined $inheritmodel) {
		$inheritmodel = 'NA';
	} elsif ($inheritmodel ne 'compoundhet' && $inheritmodel ne 'compoundhetmosaic' && $inheritmodel ne 'unique') {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --model defined but not valid (should be compoundhet or compoundhetmosaic or unique)\n");
	}
	if (!defined $excludeGVSfunction) {
		pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --excludefunction not defined\n");
	}
	if (!defined $debugmode) {
		$debugmode = 0;
	} else {
		print STDOUT "Debug mode (#$debugmode) is on\n";
	}
	if ($filters eq 'all' || $filters eq 'any') {
		%allowedGATKfilters = map {$_ => 1} qw(QUALFilter QDFilter LowQual SBFilter PASS ABFilter LowQual HRunFilter SnpCluster QUALFilter);
	} elsif ($filters eq 'default') {
		%allowedGATKfilters = map {$_ => 1} qw(PASS SBFilter ABFilter);
	} else {
		%allowedGATKfilters = map {$_ => 1} split(',', $filters);
	} 
	if ($excludeGVSfunction eq 'default') {
		%GVStoexclude = map {$_ => 1} qw(intron intergenic utr utr-3 utr-5 near-gene-3 near-gene-5 none coding-synonymous synonymous);
	} elsif ($excludeGVSfunction eq 'NA')  {
		%GVStoexclude = map {$_ => 1} qw(NA);
	} else {
		%GVStoexclude = map {$_ => 1} split(',', $excludeGVSfunction);
	}
	
	my $logfile = "$outputfile.log";
	return ($inputfile, $outputfile, $subjectdeffile, $minhits, $filters, $isNhit, 
		$inheritmodel, $mafcutoff, $excludeGVSfunction, $cmgfreqcutoff, $mindp, $minGQ, $maxmismatchesperfamily, 
		$debugmode, \%allowedGATKfilters, \%GVStoexclude, $logfile);
}

sub printParamstoLog {
	my ($log_filehandle, $minhits, $maxmismatchesperfamily, $mindp, $minGQ, $isNhit, $gvstoexclude_string, $allowedgatk_string, $mafcutoff, $cmgfreqcutoff, $inheritmodel) = @_;

	print $log_filehandle "Requiring hits in gene in at least $minhits subjects/families\n";
	print $log_filehandle "Allowing up to $maxmismatchesperfamily subjects per family to not have the correct genotype\n";
	# if (@dpcolumns && @qualcolumns) {
		print $log_filehandle "Genotypes require at least $mindp depth and at least $minGQ qual, otherwise genotype set to be missing\n";
	# } else {
		# print $log_filehandle "Per-sample depth and quality information not available\n";
	# }
	print $log_filehandle "Missing genotypes/no calls are counted as: $isNhit\n";
	print $log_filehandle "Excluding all variants with annotations: $gvstoexclude_string\n";
	print $log_filehandle "Only allow variants with GATK filter: $allowedgatk_string\n";
	print $log_filehandle "Excluding variants with MAF>$mafcutoff in ESP and/or 1000 Genomes\n";
	print $log_filehandle "Excluding variants with frequency>=$cmgfreqcutoff in CMG subjects (likely systematic error)\n";
	print $log_filehandle "Special analysis? $inheritmodel\n";
}

sub readPedigree {
	my ($subjectdeffile, $log_filehandle) = @_;
	my (%countuniquefamilies_hash, @orderedsubjects, %subjects);
	my %validvalues = (1 => 1, 2 => 1, -9 => 1, 0 => 1);
	
	########## Order of subjects in this file must correspond to order of genotype columns
	open (SUBJECTS, "$subjectdeffile") or die "Cannot read $subjectdeffile: $!.\n";
	while (<SUBJECTS>) {
		$_ =~ s/\s+$//;					# Remove line endings
		if ($_ =~ /^#/) {
			next;						# skip header lines
		}
		print $log_filehandle "$_\n";
		my ($familyid, $subjectid, $father, $mother, $sex, $phenotype, $relation, $desiredgeno, @other) = split(/[ \t]+/, $_);
		if ($desiredgeno ne 'het' && $desiredgeno ne 'ref' && $desiredgeno ne 'alt' && $desiredgeno ne 'any') {
			print $log_filehandle "$familyid	$subjectid has an illegal desired genotype ($desiredgeno)\n";
			die "$familyid	$subjectid has an illegal desired genotype ($desiredgeno)\n";
		}
		if (!defined $validvalues{$phenotype}) {
			print $log_filehandle "$familyid	$subjectid has an illegal phenotype value ($phenotype)\n";
			die "$familyid	$subjectid has an illegal phenotype value ($phenotype)\n";
		}
		if (!defined $validvalues{$sex}) {
			print $log_filehandle "$familyid	$subjectid has an illegal sex value ($sex)\n";
			die "$familyid	$subjectid has an illegal sex value ($sex)\n";
		}
		$subjects{$subjectid} = [$familyid, $father, $mother, $sex, $relation, $desiredgeno];
		push(@orderedsubjects, $subjectid);
		if ($subjectid !~ '#') {
			push(@{$countuniquefamilies_hash{$familyid}}, $relation);
		} else {
			print $log_filehandle "$subjectid is being skipped in this analysis\n";
		}
	}
	close SUBJECTS;
	my $countuniquefamilies = scalar(keys %countuniquefamilies_hash);
	print $log_filehandle "\n";
	
	if ($debugmode >= 1) {
		print STDOUT "Read in $countuniquefamilies families from $subjectdeffile: subject data order: @orderedsubjects\n";
	}
	return (\%countuniquefamilies_hash, \@orderedsubjects, \%subjects, $countuniquefamilies);
}

sub parseHeader {
	##
	# add code to check for undefined columns
	##
	my @header = @_;
	my ($polyphencol, $phastconscol, $gerpcol, $gatkfiltercol, $freqinCMGcol, $freqinOutsidecol, @genotypecolumns, @qualcolumns, @dpcolumns, @keepcolumns);
	for (my $i=0; $i<=$#header; $i++) {
		my $columnname = $header[$i];
		if ($columnname =~ /Qual/i) {
			push(@qualcolumns, $i);
		} elsif ($columnname =~ /Depth/i && $columnname ne 'averageDepth') {
			push(@dpcolumns, $i);
		} elsif ($columnname =~ /Gtype/i) {
			$columnname =~ s/Gtype//; 
		}
		# if ($columnname !~ /Qual/i && $columnname !~ /Depth/i) {
			push(@keepcolumns, $i);
		# }
		if (defined $subjects{$columnname} || defined $subjects{"#$columnname"}) {
			push(@genotypecolumns, $i);
		}
		if ($columnname =~ /Freqin/i) {
			if ($columnname =~ /FreqinCMG/i) {
				$freqinCMGcol = $i;
			} elsif ($columnname =~ /FreqinOutsidePop/i) {
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
	return ($polyphencol, $phastconscol, $gerpcol, $gatkfiltercol, $freqinCMGcol, $freqinOutsidecol, \@genotypecolumns, \@qualcolumns, \@dpcolumns, \@keepcolumns);
}

sub passGATKfilters {
	my ($allowedfilters_ref, $filterset) = @_;
	my $keep = 1;
	if ($filterset ne '.') {
		my @thesefilters = split(';', $filterset);
		foreach my $filter (@thesefilters) {
			if (!defined ${$allowedfilters_ref}{$filter}) {
				$keep = 0;
			} 
		}
	}
	return $keep;
}

sub shouldfunctionFilter {
	my ($GVStoexclude_ref, $functionimpact) = @_;
	my $score = 0;
	my @filters = split(",", $functionimpact);
	foreach my $eachfilter (@filters) {
		if (defined ${$GVStoexclude_ref}{$eachfilter}) {
			$score += 1;
		}
	}

	if ($score == scalar(@filters)) {
		return 1;
	} else {
		return 0;
	}
}

sub checkGenoMatch {
	my ($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr, $debugmode) = @_;		
	my @alleles;
	my $ismatch = 0;
	
	if (($genotype eq 'N' && $isNhit eq 'hit') || $desiredgeno eq 'any') {
		$ismatch = 1;
	} elsif ($genotype eq 'N' && $isNhit eq 'nohit') {
		$ismatch = 0;
	} else {
		my %altalleles = map {$_ => 1} split(",", $alt);	

		if ($inputfiletype eq 'SSAnnotation') {
			if ($genotype eq $ref) {
				@alleles = ($ref, $ref);
			} elsif ($genotype eq $alt) {
				@alleles = ($alt, $alt);
			} else {
				@alleles = ($ref, $alt);
			}
			if ($vartype eq 'indel') {
				@alleles = split("/", $genotype);
			}
		} else {
			@alleles = split("/", $genotype);
		}
	
		if ($desiredgeno eq 'ref') {
			if ($alleles[0] eq $ref && $alleles[1] eq $ref) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'alt') {
			if (defined $altalleles{$alleles[0]} && defined $altalleles{$alleles[1]}) {
				$ismatch = 1;
			}
		} elsif ($desiredgeno eq 'het') {
			if (($alleles[0] eq $ref || $alleles[1] eq $ref) && (defined $altalleles{$alleles[1]} || defined $altalleles{$alleles[0]})) {
				$ismatch = 1;
			}
			if ($sex == 1 && ($chr eq "X" || $chr eq "Y") && (defined $altalleles{$alleles[0]} && defined $altalleles{$alleles[1]})) {
				$ismatch = 1;
			}
		} 
	}
	# if ($debugmode >= 3) {
		# print STDOUT "$vartype	geno=$genotype	desiredgeno=$desiredgeno	ref=$ref	alt=$alt	genoalleles=@alleles	ismatch=$ismatch\n";				# DEBUG
	# }
	return $ismatch;
}

sub failsCallQuality {
	my $qualityvalues = $_[0];
	
	my $return = 0;
	if ($qualityvalues =~ m/DP/i || $qualityvalues =~ m/GQ/i) {
		$return = 1;
	}
	return $return;
}

sub determineGenoType {
	my ($vartype, $ref, $alt, $genotype, $desiredgeno, $isNhit, $sex, $chr) = @_;		
	my @alleles;
	my $typeofgeno;
	my %altalleles = map {$_ => 1} split(",", $alt);	
	if ($inputfiletype eq 'SSAnnotation') {
		if ($genotype eq $ref) {
			@alleles = ($ref, $ref);
		} elsif ($genotype eq $alt) {
			@alleles = ($alt, $alt);
		} else {
			@alleles = ($ref, $alt);
		}
		if ($vartype eq 'indel') {
			@alleles = split("/", $genotype);
		}
	} else {
		@alleles = split("/", $genotype);
	}
	
	if ($genotype eq 'N' || $genotype eq 'NA') {
		$typeofgeno = 'N';
	} elsif ($alleles[0] eq $ref && $alleles[1] eq $ref) {
		$typeofgeno = 'homref';
	} elsif (defined $altalleles{$alleles[0]} && defined $altalleles{$alleles[1]}) {
		$typeofgeno = 'homalt';
	} elsif (($alleles[0] eq $ref || $alleles[1] eq $ref) && (defined $altalleles{$alleles[1]} || defined $altalleles{$alleles[0]})) {
		$typeofgeno = 'het';
	}
	return $typeofgeno;
}


sub vcfgeno2calls {
	my ($genotype, $ref, $alt) = @_;
	my $newgenotype;
	my @allalleles = ($ref, split(",", $alt));
	if ($genotype eq '.' || $genotype eq './.') {
		$newgenotype = 'N';
	} else {
		my @oldgeno = split("/", $genotype);
		$newgenotype = "$allalleles[$oldgeno[0]]/$allalleles[$oldgeno[1]]";
	}
	return $newgenotype;
}

sub selectRefAlt {
	### not compatible with multi-allelic SNPs
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

sub determineMatchSubject {
	my ($fatherismatch, $motherismatch, $inheritmodel, $requiredcarrier, $fathertypeofgeno, $mothertypeofgeno) = @_;
	my $matchsubj = 'reject';
	if (!defined $fathertypeofgeno) { $fathertypeofgeno = 'N'; }
	if (!defined $mothertypeofgeno) { $mothertypeofgeno = 'N'; }
	if (!defined $requiredcarrier) { $requiredcarrier = 'N'; }
	if ($debugmode >= 4) {
		print STDOUT "determineMatchSubject: $fatherismatch, $motherismatch, $inheritmodel, $requiredcarrier, $fathertypeofgeno, $mothertypeofgeno\n";
	}
	if ("$fatherismatch$motherismatch" !~ 'NA' && "$fathertypeofgeno$mothertypeofgeno" !~ 'N') {											# if both parents have genotypes
		### if ($fatherismatch+$motherismatch) == 2, then both parents would be carriers for the same rare mutation?!, in which case it would be homozygous recessive not compound het; just very very unlikely
		if (($fatherismatch+$motherismatch) == 1 || ($inheritmodel eq 'compoundhetmosaic' && ($fatherismatch+$motherismatch) == 2)) {
			if ($inheritmodel eq 'compoundhetmosaic') {
				if ($requiredcarrier eq 'father') {
					if ($fatherismatch==1) {
						$matchsubj = 'father';
					} else {
						$matchsubj = 'mother';
					}
				} elsif ($requiredcarrier eq 'mother') {
					if ($motherismatch==1) {
						$matchsubj = 'mother';
					} else {
						$matchsubj = 'father';
					}
				}
			} elsif ($inheritmodel eq 'compoundhet') {
				if ($debugmode >=3) { print STDOUT "father is $fathertypeofgeno, mother is $mothertypeofgeno: "; }
				if (($fathertypeofgeno eq 'het' || $fathertypeofgeno eq 'homalt') xor ($mothertypeofgeno eq 'het' || $mothertypeofgeno eq 'homalt')) {
					if ($debugmode >=3) { print STDOUT "pass\n"; }
					if ($fatherismatch == 1) {
						$matchsubj = 'father';
					} elsif ($motherismatch == 1) {
						$matchsubj = 'mother';
					}	
				} else {
					if ($debugmode >=3) { print STDOUT "fail\n" }
				}
			}
		}
	} elsif ($fathertypeofgeno =~ 'N' && $mothertypeofgeno =~ 'N') {							# if both parents do not have a genotype
		$matchsubj = 'either';
	} elsif ($fathertypeofgeno =~ 'N') {														# if father does not have genotype but mother does
		if ($motherismatch == 1) {
			$matchsubj = 'mother';
		} elsif ($motherismatch == 0) {
			$matchsubj = 'father';
		}
	} elsif ($mothertypeofgeno =~ 'N') {														# if mother does not have genotype but father does
		if ($fatherismatch == 1) {
			$matchsubj = 'father';
		} elsif ($fatherismatch == 0) {
			$matchsubj = 'mother';
		}
	}
	
	# need to think about this more but I think a compound het + mosaic model only matters if both parents were sequenced
	# if one or neither parent was sequenced, can filter with other inheritance models that look identical
	if ($debugmode >= 4) {
		print STDOUT "determineMatchSubject: matchsubj=$matchsubj\n";
	}
	return $matchsubj;
}
sub determinevartype {
	my ($ref, $alt) = @_;
	my $vartype = 'SNP';
	my @altalleles = split(",", $alt);
	if (length($ref)>1 || length($altalleles[0])>1) {
		$vartype = 'indel';
	}
	return $vartype;
}

sub determineInputType {
	my ($firstline, $filename) = @_;
	my $inputtype = 'SSAnnotation';
	if ($filename =~ '.vcf' || $firstline =~ 'VCFv4') {
		$inputtype = 'vcf';
	} elsif ($filename =~ 'SeattleSeqAnnotation' || $firstline =~ /^# inDBSNPOrNot/) {
		$inputtype = 'SeattleSeqAnnotation';
	}
	return $inputtype;
}

sub retrieveVCFfields {
	my ($vcfcolumn, $vcfcoltype) = @_;
	my @columncontents;
	my %vcffields;
	foreach my $name (qw(GL FG CP CG PP AAC PP CDP AFCMG AFPOP CA KP PH)) {
		$vcffields{$name} = '.';
	}
	
	if ($vcfcoltype eq 'genotype') {
		@columncontents = split(":", $vcfcolumn);
	} elsif ($vcfcoltype eq 'info') {
		@columncontents = split(";", $vcfcolumn);
	} else {
		@columncontents = split(/[;:]/, $vcfcolumn);
	}
	foreach my $column (@columncontents) {
		my ($fieldname, $value) = split("=", $column);
		if (!defined $value) {
			$value = '.';
		}
		$vcffields{$fieldname} = $value;
	}
	return \%vcffields;
}

# parse input files line by line
sub parse_SSAnnotation134_byline {
	my ($subj_gtcolumns_ref, $dpcols_ref, $qualcols_ref, @line) = @_;
	my (@subjectquals, @subjectdps);
	my @genotypecolumns = @{$subj_gtcolumns_ref};
	my @dpcolumns = @{$dpcols_ref};
	my @qualcolumns = @{$qualcols_ref};
	
	my ($chr, $pos, $vartype, $ref, $alt) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
	my $filterset = $line[$gatkfiltercol];
	my $gene = $line[7];
	my $functionimpact = $line[10];
	# my $phastcons = $line[$phastconscol];
	# my $gerp = $line[$gerpcol];
	# my $polyphen = $line[$polyphencol];
	# $inUWexomes = $line[16];
	# $UWexomescovered = $line[17];
	
	my @subjectgenotypes = @line[@genotypecolumns];
	if (@dpcolumns) {
		@subjectdps = @line[@dpcolumns];
	}
	if (@qualcolumns) {
		@subjectquals = @line[@qualcolumns];			
	}

	return ($chr, $pos, $vartype, $ref, $alt, $filterset, $gene, \@subjectgenotypes, \@subjectdps, \@subjectquals, $functionimpact);
}

sub parse_vcf_byline {
	my ($subj_columns_ref, @line) = @_;
	my (@subjectquals, @subjectdps, @subjectgenotypes);
	my ($gene, $functionimpact, $gerp, $polyphen, $phastcons, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside, $clinical, $kegg);
			
	my ($chr, $pos, $rsid, $ref, $alt) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
	my $vartype = determinevartype($ref, $alt);
	my $filterset = $line[6];
	
	my $infofields_ref = retrieveVCFfields($line[7], 'info');
	$gene = ${$infofields_ref}{'GL'};
	$functionimpact = ${$infofields_ref}{'FG'};
	$phastcons = ${$infofields_ref}{'CP'};
	$gerp = ${$infofields_ref}{'CG'};
	$polyphen = ${$infofields_ref}{'PH'};
	$aminoacids = ${$infofields_ref}{'AAC'};
	$proteinpos = ${$infofields_ref}{'PP'};
	$cdnapos = ${$infofields_ref}{'CDP'};
	$freqinCMG = ${$infofields_ref}{'AFCMG'};
	$freqinOutside = ${$infofields_ref}{'AFPOP'};
	$clinical = ${$infofields_ref}{'CA'};
	$kegg = ${$infofields_ref}{'KP'};
	# $inUWexomes = $line[16];
	# $UWexomescovered = $line[17];
	
	my @format = split(":", $line[8]);
	my ($gtcolnum, $dpcolnum, $gqcolnum, $pindel_rd, $pindel_ad);
	my $altdp = 0;
	for (my $i=0; $i<=$#format; $i++) {
		if ($format[$i] eq 'GT') {
			$gtcolnum = $i;
		} elsif ($format[$i] eq 'DP') {
			$dpcolnum = $i;
		} elsif ($format[$i] eq 'GQ') {
			$gqcolnum = $i;
		} elsif ($format[$i] eq 'RD') {
			$pindel_rd = $i;
		} elsif ($format[$i] eq 'AD') {
			$pindel_ad = $i;
		} elsif ($format[$i] eq 'DP4') {
			# samtools singlesample calls will put depth info (DP4) in INFO column; need to move to FORMAT column to create multisample vcf for this script to read
			$dpcolnum = $i;
			$altdp = 1;
		}
	}
	
	# in UW CMG files, the only time DP,GQ are missing is if the call came from Pindel, so this assumes there is info from Pindel to retrieve
	foreach my $subject (@line[@{$subj_columns_ref}]) {
		my @subjectdata = split(":", $subject);
		push(@subjectgenotypes, vcfgeno2calls($subjectdata[$gtcolnum], $ref, $alt));
		if ($subject eq './.' || $subject eq '.') {
			push(@subjectdps, 0);
			push(@subjectquals, 0);
		} else {
			if (defined $dpcolnum) {
				if ($altdp == 0) {
					push(@subjectdps, $subjectdata[$dpcolnum]);
				} else {
					my @dp4 = split(",", $subjectdata[$dpcolnum]);
					push(@subjectdps, ($dp4[0]+$dp4[1]+$dp4[2]+$dp4[3]));
				}
			} else {
				push(@subjectdps, ($subjectdata[$pindel_rd]+$subjectdata[$pindel_ad]));
			}
			if (defined $gqcolnum) {
				push(@subjectquals, $subjectdata[$gqcolnum]);
			} else {
				push(@subjectquals, 99);
			}
		}
	}
	# if ($debugmode >= 4) {
	# 	print STDOUT "$chr, $pos, $vartype, $ref, $alt, $filterset, $gene, @subjectgenotypes, @subjectdps, @subjectquals, $functionimpact, $gerp, $polyphen, $phastcons, $rsid, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside\n";			## DEBUG
	# }
	return ($chr, $pos, $vartype, $ref, $alt, $filterset, $gene, \@subjectgenotypes, \@subjectdps, \@subjectquals, $functionimpact, $gerp, $polyphen, $phastcons, $rsid, $aminoacids, $proteinpos, $cdnapos, $freqinCMG, $freqinOutside, $clinical, $kegg);
}

sub parse_SeattleSeqAnnotation134_byline {
	# if ($inputfile =~ m/snps/) {
	# 	my $sampleallelestring = $line[5];
	# 	my @samplealleles = split("/", $sampleallelestring);
	# 	if (!defined $samplealleles[1]) {
	# 		$samplealleles[1] = 'N';
	# 		next;
	# 	}
	# 	$gene = $line[20];
	# 	$filterset = 'NA';
	# 	$functionimpact = $line[8];
	# 	# $indbsnp = $line[0];
	# 	$gerp = $line[$gerpcol];
	# 	$phastcons = $line[$phastconscol];
	# 	$polyphen = $line[$polyphencol];		
	# 	$vartype = 'SNP';
	# 	($ref, $alt) = selectRefAlt($line[3], $samplealleles[0], $samplealleles[1]);
	# 	@subjectgenotypes = split(",", $line[4]);
	# }
	# if ($inputfile =~ m/indels/) {
	# 	# ($vartype, $ref, $alt) = ('indel', $line[3], $line[4]);		# complicated parsing; not implemented
	# 	die "Have not implemented parsing of this file format yet\n";
	# }
}






################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


filter_model.pl - filter annotated vcfs or SeattleSeq annotation files based on various Mendelian inheritance models


=head1 SYNOPSIS


perl B<filter_model.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> I<input file>	

	input file (either vcf or SSAnnotation)

=item B<--out> I<output file>

	name of output file

=item B<--subjectreq> I<pedigree file>

	ped file for samples IN ORDER of appearance in input file
	first 6 columns: family subject father mother sex phenotype
	7th column (temporary until I improve pedigree code): person's relationship in family (father, mother, child, etc).  Children must have distinct names (proband1, proband2, etc)
	8th column: desired genotype given inheritance model

=item B<--minhits> I<number>
	
	minimum number of families/individuals with hits
	
=item B<--N> I<hit/nohit>

	how to treat missing genotypes
		
=item B<--maxmissesperfamily> I<number>

	max number of subjects per family with genotypes not matching model

=item B<--GATKkeep> I<list>

	GATK quality filters that should be kept as possible hits
	Possible values: comma-delimited list with no spaces OR 'all' OR 'default' (default=PASS,SBFilter,ABFilter)
		
=item B<--excludefunction> I<list>
	
	GVS variant function classes to be excluded
	Possible values: comma separated list with no spaces OR 'default' (default=intron,intergenic,utr,utr-3,utr-5,near-gene-3,near-gene-5,none) OR 'NA'
		
=item B<--mafcutoff> I<value>

	cutoff MAF for filtering out common variants using 1000 Genomes and/or ESP; any var with freq > cutoff is excluded
	optional, default value = 0.01
		
=item B<--errorcutoff> I<value>

	cutoff frequency for filtering out systematic errors based on frequency in all CMG exomes to date
	optional, default value = 0.2

=item B<--model> I<inheritance model>

	optional, not necessary for homozygous recessive, dominant, de novo
	Possible values: 'compoundhet' OR 'compoundhetmosaic' OR 'unique' if filtering for unique de novo mutations (never seen before and not recurrent)>

=item B<--dp> I<number>

	minimum DP for subject genotype; otherwise considered missing
	optional, default=8
		
=item B<--gq> I<number>

	minimum GQ for subject genotype; otherwise considered missing
	optional, default=20

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


This script will filter genotypes for subjects from a vcf or SSAnnotation file and filter for variants matching a desired Mendelian inheritance model.  The inheritance model is specified in an extended pedigree-like file and with the --model argument.


=head1 CAVEATS


When allowing missing (which includes N as well as low DP/GQ) genotypes to be hits, think hard about allowing misses to be matches (can end up with a lot of variants).


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut

