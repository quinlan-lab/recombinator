<<DOC
After running an unphased VCF rough recombinator, we get a list of all sites
that are informative and their states. We are interested in sites where there is
a state-change. We want to phase the sites that bound that change so we can know
which sib inherited it.

We get all the unique sites that bound a change, combine them with the omni 2.5
sites, extract the union of those sites in a VCF that can then be sent to shapeit
for phasing.

This requires:
+ shapeit (v2.r837 or greater)
+ plink
+ bcftools
+ gargs
+ tabix 

DOC

set -euo pipefail

export prefix=$1 # this is the --prefix sent to recombinator.py
export vcf=$2 # this is the --vcf set to recombinator.py
export ped=$3

extract() {
	# find sites where state(i) != state(i + 1) and report both of those locations.
	zcat $1 | awk 'NR == 2 { state=$6 } (NR > 1 && $6 != state) { print last; print $1":"$2-1"-"$2; state=$6 }(NR > 1) { last=$1":"$2-1"-"$2;}' | uniq
}
export -f extract

ls $prefix/*/*/*{dad,mom}.bed.gz | gargs -s -p 32 "extract {}" | sort -t':' -u -k1,1 -k2,2n > $prefix/informative.sites

(cat $prefix/informative.sites; zcat data/omni-2.5.sites.gz ) | sort -ut':' -k1,1 -k2,2n > $prefix/all-informative.sites

(bcftools view -h $vcf;
cat $prefix/all-informative.sites \
	| gargs -s -n 100 -p 20 "bcftools view -H $vcf {}" \
	| awk '$7 == "PASS" || $7 == "." && $6 > 10' \
	| awk 'BEGIN{x=0;} $0 ~/^#/{ if(x==0) {print;} next}{x=1; print $0 | "LC_ALL=C sort -u --buffer-size 3G -k1,1 -k2,2n"}' \
	) \
	| bcftools view -m2 -M2 -c1 -o $prefix/sites.tophase.vcf.gz -O z

# see: https://github.com/quinlan-lab/recombinator/blob/bdf6c28f5a938683a076b9eddebd768864bde34c/duohmm.sh
mkdir -p $prefix/shapeit/

plinkify() {
	chrom=$1
	plink \
		--keep-allele-order \
		--chr $chrom --geno 0.05 \
		--mind 0.05 --vcf-half-call m \
		--biallelic-only --vcf $prefix/sites.tophase.vcf.gz \
		--make-bed \
		--out $prefix/shapeit/$chrom

}

run_shapeit(){
	chrom=$1
	shapeit \
		--aligned \
		--duohmm \
		-B $prefix/shapeit/$chrom \
		-T 2 \
		-M $MAP \
		-W 2 \
		-O $prefix/shapeit/$chrom.duohmm

	# convert to phased vcf.
	shapeit -convert --input-haps \
		$prefix/shapeit/$chrom.duohmm \
		--output-vcf $prefix/shapeit/$chrom.phased.vcf \
	bgzip -f $prefix/shapeit/$chrom.phased.vcf
	tabix $prefix/shapeit/$chrom.phased.vcf.gz

	python recombinator.py \
		--min-gq 20 \
		--min-depth 16 \
		--region $chrom \
		--vcf $prefix/shapeit/$chrom.phased.vcf.gz \
		--ped $ped \
		--prefix $prefix/phased/
}

export -f plinkify
export -f run_shapeit

seq 1 22 \
	| gargs -s -p 22 "plinkify {} && run_shapeit {}"
