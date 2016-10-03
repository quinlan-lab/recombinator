#rm -rf k
#python -m recombinator recombinator --ped K34175.ped --vcf K34175.subset.snps.phased.vcf.gz --prefix k
#python -m recombinator denovos K34175.ped K34175.subset.snps.phased.vcf.gz > k.denovos.bed

<<NO_XOS
python -m recombinator enrichment \
	--figure k.enrich.png \
	$XOS \
	k.denovos.bed
NO_XOS
