#!/usr/bin/env perl
#
# Description: Look for mutations matching a particular model
#
#
#
# Created by Jessica Chong on 2012-10-15.


use strict;
use warnings;
use Getopt::Long;


my ($inputfile, $outputfile, $subjectdeffile, $allowedmisses, $filters, $isNhit, $inheritmodel);
my $capturearray = 'NA';


GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'subjectreq=s' => \$subjectdeffile,
	'misses=i' => \$allowedmisses,
	'GATKkeep=s' => \$filters,
	'N=s' => \$isNhit,
	'capture:s' => \$capturearray,
	'model:s' => \$inheritmodel,
);

if (!defined $inputfile) {
	optionUsage("option --in not defined\n");
} elsif (!defined $outputfile) {
	optionUsage("option --out not defined\n");
} elsif (!defined $subjectdeffile) {
	optionUsage("option --subjectreq not defined\n");
} elsif (!defined $allowedmisses) {
	optionUsage("option --misses not defined\n");
} elsif (!defined $filters) {
	optionUsage("option --GATKkeep not defined\n");
} elsif (!defined $isNhit) {
	optionUsage("option --N not defined\n");
}

if (!defined $inheritmodel) {
	$inheritmodel = 'NA';
}

my %allowedGATKfilters;
if ($filters eq 'all' || $filters eq 'any') {
	%allowedGATKfilters = map {$_ => 1} qw(QUALFilter QDFilter LowQual SBFilter PASS ABFilter LowQual HRunFilter SnpCluster QUALFilter);
} else {
	%allowedGATKfilters = map {$_ => 1} split(',', $filters);
}


my %GVStoexclude = (
	'intron' => 1,
	'intergenic' => 1,
	# 'coding-synonymous' => 1,
	# 'coding-synonymous-near-splice' => 1,
	# 'utr-3' => 1,
	# 'utr-5' => 1,
	'near-gene-3' => 1,
	'near-gene-5' => 1,
);



my %countuniquefamilies_hash;
my @orderedsubjects;
my %subjects;
open (SUBJECTS, "$subjectdeffile") or die "Cannot read $subjectdeffile: $!.\n";
while (<SUBJECTS>) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($familyid, $subjectid, $father, $mother, $relation, $desiredgeno) = split("\t", $_);
	$subjects{$subjectid} = [$familyid, $father, $mother, $relation, $desiredgeno];
	push(@orderedsubjects, $subjectid);
	if ($subjectid !~ '#') {
		$countuniquefamilies_hash{$familyid} = 1;
	}
}
close SUBJECTS;
my $countuniquefamilies = scalar(keys %countuniquefamilies_hash);


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
for (my $i=0; $i<=$#header; $i++) {
	if (defined $subjects{$header[$i]} || defined $subjects{"#$header[$i]"}) {
		push(@subjectcolumns, $i);
	}
}
print OUT "$headerline\n";
while ( <FILE> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	if ($_ =~ '#') {
		if ($_ =~ 'UnifiedGenotyper="analysis_type=UnifiedGenotyper') {
			if ($_ =~ 'nimblegen_solution_bigexome_2011') {
				$capturearray = 'bigexome';
			} elsif ($_ =~ 'v2') {
				$capturearray = 'v2';
			}
		} else {
			next;	
		}
	}
	$countinputvariants++;
	
	if ($printparams == 0) {
		if ($capturearray eq 'NA') {
			$capturearray = 'bigexome';
		}
		print "Excluding all variants with annotations: ".join(" ", keys %GVStoexclude)."\n";
		print "Excluding variants in systematic error file for $capturearray\n";
		print "Only allow variants with GATK filter: ".join(" ", keys %allowedGATKfilters)."\n";
		$printparams = 1;
	}
	
	my @line = split ("\t", $_);
	my $countmatches = 0;
	my ($chr, $pos, $vartype, $ref, $alt);
	my @subjectgenotypes;
	my ($gene, $functionimpact, $gerp, $polyphen, $phastcons);
	my $indbsnp = 'NA';
	my $inUWexomes = 'NA';
	my $UWexomescovered = 'NA';
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
			$indbsnp = $line[0];
			$gerp = $line[17];
			$phastcons = $line[16];
			$polyphen = $line[14];		
			$vartype = 'SNP';
			($ref, $alt) = selectRefAlt($line[3], $samplealleles[0], $samplealleles[1]);
			@subjectgenotypes = split(",", $line[4]);
		}
		if ($inputfile =~ m/indels/) {
			# ($vartype, $ref, $alt) = ('indel', $line[3], $line[4]);		# complicated parsing; not implemented
		}
	} elsif ($inputfile =~ m/SSAnnotation/) {
		($chr, $pos, $vartype, $ref, $alt) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		$filterset = $line[5];
		$gene = $line[7];
		@subjectgenotypes = @line[@subjectcolumns];
		$indbsnp = $line[8];
		$functionimpact = $line[10];
		$phastcons = $line[13];
		$gerp = $line[14];
		$polyphen = $line[15];
		$inUWexomes = $line[16];
		$UWexomescovered = $line[17];
	} 
	
	if (noveltyfunctionFilter(\%GVStoexclude, $functionimpact, $indbsnp, $inUWexomes, $UWexomescovered) && ($filterset eq 'NA' || checkGATKfilters(\%allowedGATKfilters, $filterset))) {
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
		my @matchingfamilyids;
		while (my ($familyid, $thisfamily_ref) = each %checkfamilies) {
			my %thisfamily = %{$thisfamily_ref};

			my $thisfamilymatch = 0;
			my $familysize = scalar(keys %thisfamily);

			if ($inheritmodel eq 'compoundhet' && $familysize >= 3) {							# assumes trios; doesn't handle more than one affected child per family (yet)
				if (($thisfamily{'father'}+$thisfamily{'mother'}) == 1 && $thisfamily{'child'} == 1) {
					$countfamiliesmatch += 1;
					${$genehits{$gene}{'families'}{$familyid}}[0] += $thisfamily{'father'};		#
					${$genehits{$gene}{'families'}{$familyid}}[1] += $thisfamily{'mother'};		#
					${$genehits{$gene}{'families'}{$familyid}}[2] += $thisfamily{'child'};		#
				}
			} else {
				while (my ($relation, $thismatch) = each %thisfamily) {
					if ($thismatch != -1) {
						$thisfamilymatch += $thismatch;
					}
				}

				if ($thisfamilymatch == $familysize) {
					$countfamiliesmatch += 1;
					push(@matchingfamilyids, $familyid);
					$genehits{$gene}{'families'}{$familyid} += 1;
				}
			}
		}

		if ($countfamiliesmatch > 0) {
			my $data = join("\t", @line);
			my $familiesmatching = join(',', @matchingfamilyids)
			push(@{$genehits{$gene}{'hits'}}, "$data\t$familiesmatching");
		}
	}
}
close FILE;



# Check putative hit list for all genes and count number of valid hits

my $countoutputvariants = 0;
my $countkeptvariants = 0;
my $counterrorvariants = 0;
while (my ($gene, $results_ref) = each %genehits) {
	my %hitdata = %{$results_ref};
	
	if (defined $hitdata{'hits'}) {
		# filter out common variants and systematic errors
		my @filteredoutput;
		my @output = @{$hitdata{'hits'}};
		foreach my $hit (@output) {
			my @thishit = split("\t", $hit);
			my ($chr,$pos,$vartype,$ref,$alt) = @thishit[0..4];
			my $familiesmatching = pop(@thishit);
			my %matchingfamilyids = map {$_ => 1} split(/,/, $familiesmatching);
			my $iserror = isSystematicError($chr,$pos,$vartype,$ref,$alt,$capturearray);
			my $iscommon = isCommonVar($chr,$pos,$vartype,$ref,$alt);
		
			if ($iserror==1 || $iscommon==1) {
				while (my($familyid, $familyhits_ref) = each %{$hitdata{'families'}}) {
					if (defined $matchingfamilyids{$familyid}) {
						if (ref($familyhits_ref) eq 'SCALAR') {
							$hitdata{'families'}{$familyid}--;
						} elsif (ref($familyhits_ref) eq 'ARRAY') {				### EXCLUDE VARIANT from trios/families under compound het model
							my @familyhits = @{$familyhits_ref};
							## subtract hit from the correct people in this family
						}
					}
				}
			} else {
				push(@filteredoutput, @thishit);
			}
		}
		
		# after filtering for common/error variants, account for possible compound het model, then recount matching families
		my $countfamiliesmatch = 0;
		while (my($familyid, $familyhits_ref) = each %{$hitdata{'families'}}) {		
			# compound het model only
			if (ref($familyhits_ref) eq 'ARRAY') {													# if at least a trio
				my @familyhits = @{$familyhits_ref};
				if ($familyhits[0] >= 1 && $familyhits[1] >= 1 && $familyhits[2] >= 2) {			# if at least one variant each from mom and dad
					$countfamiliesmatch++;
				}
			}	# end compound het model
			if (ref($familyhits_ref) eq 'SCALAR') {
				if ($familyhits_ref > 0) {
					$countfamiliesmatch++;
				}
			}
		}
	
		if ($countfamiliesmatch >= ($countuniquefamilies-$allowedmisses)) {
			foreach my $hit (@filteredoutput) {
				$countoutputvariants++;
				print OUT "$hit\n";
			}
		}
	}
}


print "Matched $countoutputvariants out of $countinputvariants\n";



sub isCommon {
	# remove common variation (ESP and 1kG)
	return 0;
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


sub noveltyfunctionFilter {
	my ($GVStoexclude_ref, $functionimpact, $indbsnp, $inUWexomes, $UWexomescovered) = @_;
	my $score = 0;
	if (!defined ${$GVStoexclude_ref}{$functionimpact}) {
		$score += 1;
	}
	if ($indbsnp eq 'none' || $indbsnp eq 'NA' || $indbsnp eq '.' || $indbsnp eq '0') {
		$score += 1;
	}
	if ($inUWexomes == 0 || $UWexomescovered == 0 || $inUWexomes eq 'NA') {
		$score += 1;
	}
	if ($score == 3) {
		return 1;
	} else {
		return 0;
	}
}


sub isSystematicError {
	my ($targetchr, $targetpos, $targettype, $targetref, $targetalt, $capturearray) = @_;
		
	my $errorpath = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references/systematic_error/2012_oct';
	my ($snperrorsfile, $indelerrorsfile);
	if ($capturearray eq 'bigexome') {
		$snperrorsfile = "$errorpath/snv.bigexome.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.bigexome.vcf.gz";	
	} elsif ($capturearray eq 'v2') {
		$snperrorsfile = "$errorpath/snv.v2.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.v2.vcf.gz";
	} else {
		$snperrorsfile = "$errorpath/snv.$capturearray.vcf.gz";	
		$indelerrorsfile = "$errorpath/indels.$capturearray.vcf.gz";
	}

	my @errors;
	if ($targettype =~ m/snp/i) {
		@errors = `tabix $snperrorsfile $targetchr:$targetpos-$targetpos`;
	} elsif ($targettype =~ m/indel/i) {
		@errors = `tabix $indelerrorsfile $targetchr:$targetpos-$targetpos`;
	}
	
	my $iserror = 0;
	foreach my $error (@errors) {
		my ($errorchr, $errorpos, $errortype, $errorref, $erroralt, @errordata) = split("\t", $error);
		if ($errorpos == $targetpos && $errorref eq $targetref && $erroralt eq $targetalt) {
			$iserror = 1;
		}
	}
	return $iserror;
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
	print "\t--subjectreq\toutput file\n";
	print "\t--misses\toutput file\n";
	print "\t--GATKkeep\tGATK quality filters that should be kept, comma-delimited with no spaces\n";
	print "\t--N\thit or nothit (how should we count missing genotypes)\n";
	print "\t--capture\tcapture array (optional: bigexome or v2)\n";
	print "\t--model\toptional: compoundhet)\n";
	die;
}

