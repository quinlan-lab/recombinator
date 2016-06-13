VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=4fams
export MAP=/scratch/ucgd/lustre/u6000771/Projects/src/shapeit/example/genetic_map.txt
PED=~u6000771/Data/ssc_519.ped
export VCF
export prefix

#seq 1 22 | xargs -P 22 -I {} sh -c "bcftools view -m2 -M2 -a -S samps.txt $VCF -c 1:nref -O z -o simons.chr{}.$prefix.vcf.gz {}"

set -eo nounset
mkdir -p results/

module load shapeit/v2.r837
module load duohmm/v0.1.7

export chrom=2
plink --recode --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --me 1 1 --vcf simons.chr$chrom.$prefix.vcf.gz --make-bed --out chr$chrom.$prefix

cut -f 1-6 fam.ped | grep -v ^# > chr$chrom.$prefix.fam

shapeit \
	-B chr$chrom.$prefix \
	--duohmm \
	-M $MAP \
	-W 5 \
	--output-max results/duohmm-$chrom-$prefix \
	--output-graph results/duohmm-$chrom-$prefix.graph \



seq 1 10 | xargs -P 10 -I{} sh -c  "shapeit -convert --input-graph results/duohmm-$chrom-$prefix.graph --output-sample results/sim$chrom.${prefix}.{} --seed {} && duohmm -H results/sim$chrom.${prefix}.{} -M $MAP -R results/sim$chrom.${prefix}.{}.rec"

mapavg.py results/sim${chrom}.${prefix}.*.rec > duohmm-${chrom}-${prefix}-recombinations.txt
