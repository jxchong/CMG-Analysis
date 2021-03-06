perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq dominant.subjectreq --minhits 1 --maxmissesperfamily 0 --GATKkeep all --N hit --excludefunction intergenic --out testout.Nhit.dominant.tsv --mafcutoff .01 --errorcutoff .2
perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq homozygrec.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.Nhit.homozygrec.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 0



perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.Nhit.compoundhet.tsv --mafcutoff .01 --errorcutoff .2 --model compoundhet --maxmissesperfamily 0
perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.Nhit.1miss.compoundhet.tsv --mafcutoff .01 --errorcutoff .2 --model compoundhet --maxmissesperfamily 1 --debug 4
perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet.subjectreq --minhits 1 --GATKkeep all --N nohit --excludefunction intergenic --out testout.Nnohit.compoundhet.tsv --mafcutoff .01 --errorcutoff .2 --model compoundhet --maxmissesperfamily 0


perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet-mosaic.subjectreq --minhits 1 --GATKkeep all --N hit --excludefunction intergenic --out testout.Nhit.1miss.compoundhet-mosaic-father.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 1 --model compoundhetmosaic



perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet-mosaic.subjectreq --minhits 1 --GATKkeep all --N nohit --excludefunction intergenic --out testout.Nnohit.compoundhet-mosaic-father.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 0 --model compoundhetmosaic
perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet.subjectreq --minhits 1 --GATKkeep all --N nohit --excludefunction intergenic --out testout.Nnohit.1miss.compoundhet.tsv --mafcutoff .01 --errorcutoff .2 --model compoundhet --maxmissesperfamily 1
perl ~/bin/filter_model.pl --in TESTDATA.SSAnnotation.tsv --subjectreq compoundhet-mosaic.subjectreq --minhits 1 --GATKkeep all --N nohit --excludefunction intergenic --out testout.Nnohit.1miss.compoundhet-mosaic-father.tsv --mafcutoff .01 --errorcutoff .2 --maxmissesperfamily 1 --model compoundhetmosaic

