perl filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq dominant.subjectreq --minhits 1 --maxmissesperfamily 0 --GATKkeep all --N hit --excludefunction intergenic --out testout.dominant.tsv --mafcutoff .01 --errorcutoff .2
perl filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq homozygrec.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.homozygrec.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 0
perl filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.compoundhet.tsv --mafcutoff .01 --errorcutoff .2 --model compoundhet --maxmissesperfamily 0
perl filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet-mosaic.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.compoundhet-mosaic-father.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 0 --model compoundhetmosaic
