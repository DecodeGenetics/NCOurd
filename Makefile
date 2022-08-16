SHELL   := /bin/bash

# The python used might need to be changed
python  := python3
NBINOM  := $(python) bin/EMmodel.py --mink=0 --dclass=nbinom --detfile=
GEOM    := $(python) bin/EMmodel.py --mink=0 --dclass=geom --detfile=

penetrance := 0.9
F45   := $(shell for n in 1 2 3; do echo run/nbinom.F45_$${n}; done)

.SECONDARY:

all: 

######
# Data that needs to be downloaded and placed in the directory extdata
# F2_genotypes.vcf.gz
# F45_genotypes.vcf.gz
# GeneConversionsRanLi2019.csv
extdata/%:
	echo '    To download dataset from Ran, Li et al. (2019) =    go to https://www.nature.com/articles/s41467-019-11675-y =    and create the file $@' | sed 's/=/\n/g'

# Used to calculate penetrance in the RanLi19 dataset 
MarkersInBetween: extdata/GeneConversionsRanLi2019.csv
	$(python) bin/inBeweeners.py extdata/GeneConversionsRanLi2019.csv --markerdir=data > $@

printPenetrance: MarkersInBetween
	awk 'BEGIN {gcBetween=0; allBetween=0; print "gcBetween","allBetween","penetrance"} NR>1 {gcBetween+=$$8; allBetween+=$$9} END {print gcBetween,allBetween, gcBetween/allBetween}' $<

# Build informative marker files
data/%.markerset: extdata/%_genotypes.vcf.gz
	echo "Chr pos" | sed 's/ /\t/' > $@
	zcat $< | grep -v "^#" | awk '{print $$1"\t"$$2}' >> $@

# Build tract files
data/RanLi.tracts: data/F2.markerset data/F45.markerset extdata/GeneConversionsRanLi2019.csv
	$(python) bin/buildTracts.py --markerdir=data extdata/GeneConversionsRanLi2019.csv > $@

# Tract functions calculated from tract files
data/RanLi.tracts.p$(penetrance).twdet.gz: data/RanLi.tracts
	$(python) bin/tracts2twdet.py $< --penetrance=$(penetrance)
	touch $@

# Detection function calculated from informative markers
data/%.markerset.p$(penetrance).detdf: data/%.markerset
	$(python) bin/Markers2DF.py --penetrance=$(penetrance) $< 

# Run EMmodel to estimate the parameters for the length distribution of NCO events as a mixture of negative binomial distributions
run/nbinom.F45_%: data/F45.markerset.p$(penetrance).detdf data/RanLi.tracts.p$(penetrance).twdet.gz
	mkdir -p run
	$(NBINOM)$^ --ngroups=$* > $@

# Run EMmodel to estimate the parameter for the length distribution of NCO events as a single geometric distributions
run/geom.F45_1: data/F45.markerset.p$(penetrance).detdf data/RanLi.tracts.p$(penetrance).twdet.gz
	mkdir -p run
	$(GEOM)$^ --ngroups=1 > $@

run/EM.nbinom.done: $(F45)
	touch $@

run/EM.geom.done: run/geom.F45_1
	touch $@

# Summary
GENOME:=2500000000

run/nbinom.F45_2.summary: run/nbinom.F45_2
	python  bin/results2NCOs.py --detdf=data/F45.markerset.p0.9.detdf $< --genome=$(GENOME) > $@

# likelihood ratio test to determine which model best fits the data out for 1 gemetric and 1,2 and 3 negative binomial distributions
likelihoodratio: run/EM.geom.done run/EM.nbinom.done
	$(python) bin/likelihoodRatio.py run/geom.F45_1 $(F45) > $@ 


### Use gene conversion regions instead of whole genome ###

# Define gene conversion regions
data/RanLi.regions:
	cat data/RanLi.tracts | awk 'NR>1 {OFS="\t";print $$1,$$2-100000,$$3+100000}' | sed -e '1iChr\tlbd\tubd' > $@

# Detection function calculated from informative markers
data/%.markerset.regions.p$(penetrance).detdf: data/%.markerset data/RanLi.regions
	$(python) bin/Markers2DF.py --penetrance=$(penetrance) --regions=data/RanLi.regions $< --output=$@

# Run EMmodel for with region detection function
run/nbinom.F45.regions_2: data/F45.markerset.regions.p$(penetrance).detdf data/RanLi.tracts.p$(penetrance).twdet.gz
	$(python) bin/EMmodel.py --detfile=data/F45.markerset.regions.p$(penetrance).detdf data/RanLi.tracts.p$(penetrance).twdet.gz --ngroups=2 --mink=0 > $@

# Summary for regions
run/nbinom.F45.regions_2.summary: run/nbinom.F45.regions_2
	python  bin/results2NCOs.py --detdf=data/F45.markerset.regions.p0.9.detdf $< --genome=$$(awk 'BEGIN {genome=0}  NR>1 {genome=genome+$$3-$$2} END {print genome}' data/RanLi.regions) > $@

