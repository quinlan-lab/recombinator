zgrep -Pv "^#|random|^GL" ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz \
	| sed -e 's/^chr//'  \
	| awk '($7 == "PASS"){split($8,a,";"); split(a[1],b,"="); ac=b[2]; if(ac > 30){ print $1":"$2"-"$2+length($4) }}' \
	| gzip -9c > omni-2.5.sites.gz
	
