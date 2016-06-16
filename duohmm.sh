VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=4fams
PED=~u6000771/Data/ssc_519.ped
export VCF
export prefix

#seq 1 22 | xargs -P 22 -I {} sh -c "bcftools view -m2 -M2 -a -S samps.txt $VCF -c 1:nref -O z -o simons.chr{}.$prefix.vcf.gz {}"

set -eo nounset
mkdir -p results/

module load shapeit/v2.r837
module load duohmm/v0.1.7

# wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
export chrom=2
export MAP=/uufs/chpc.utah.edu/common/home/u6000771/Projects/src/shapeit/example/genetic_map_b37/genetic_map_chr${chrom}_combined_b37.txt
# the hapmap ones have an extra column
#cut -f 2- /scratch/ucgd/lustre/u6000771/Projects/src/shapeit/example/genetic_map_GRCh37_chr${chrom}.txt > $MAP

export name=chr$chrom.$prefix
<<DONE
plink --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf simons.$name.vcf.gz --make-bed --out $name
cut -f 1-6 fam.ped | grep -v ^# > $name.fam

# remove mendelian errors
plink --noweb \
     --bfile $name \
     --me 1 1 \
     --set-me-missing \
     --make-bed \
     --out $name-nomendel

# phase without using ped info.
shapeit \
	--noped \
	-B $name-nomendel \
	-T 1 \
	-M $MAP \
	-W 5 \
	--output-max results/duohmm-$name \
	--output-graph results/duohmm-$name.graph
set -xo nounset

duohmm \
	-H results/duohmm-$name \
	-M $MAP \
	-O results/duohmm-$name-corrected

DONE

#seq 1 10 | xargs -P 10 -I{} sh -c  "shapeit -convert --input-graph results/duohmm-$chrom-$prefix.graph --output-sample results/sim$chrom.${prefix}.{} --seed {} && duohmm -H results/sim$chrom.${prefix}.{} -M $MAP -R results/sim$chrom.${prefix}.{}.rec"
set -x
seq 1 10 | xargs -P 10 -I{} sh -c  "shapeit -convert --input-haps results/duohmm-$name-corrected --output-haps results/duohmm-$name-haps-{} --output-sample results/sim$name.{} --seed {} && duohmm -H results/duohmm-$name-haps-{} -M $MAP -R results/sim$name.{}.rec"

mapavg.py results/sim${name}.*.rec > duohmm-${name}-recombinations.txt
