VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=simons
PED=~u6000771/Data/ssc_519.ordered.ped
export VCF
export prefix
mkdir -p data/

#seq 1 22 | xargs -P 22 -I {} sh -c "bcftools view -m2 -M2 $VCF -c 1:nref -O z -o data/$prefix.chr{}.vcf.gz {}"
#exit

set -exo nounset
mkdir -p results/

module load shapeit/v2.r837
module load duohmm/v0.1.7

# wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
export chrom=22
export MAP=/uufs/chpc.utah.edu/common/home/u6000771/Projects/src/shapeit/example/genetic_map_b37/genetic_map_chr${chrom}_combined_b37.txt
# the hapmap ones have an extra column
#cut -f 2- /scratch/ucgd/lustre/u6000771/Projects/src/shapeit/example/genetic_map_GRCh37_chr${chrom}.txt > $MAP

export name=chr$chrom.$prefix
vcf=data/$prefix.chr${chrom}.vcf.gz

<<DONE
plink --chr $chrom --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf simons.$name.vcf.gz --make-bed --out $name

plink --chr $chrom --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf $vcf --make-bed --out $name

cut -f 1-6 $PED | grep -v ^# > ${name}.fam

# remove mendelian errors
plink --noweb \
     --bfile $name \
     --me 1 1 \
     --set-me-missing \
     --make-bed \
     --out $name-nomendel

DONE

# phase without using ped info.
shapeit \
	--noped \
	-B $name-nomendel \
	-T 1 \
	-M $MAP \
	-W 5 \
	--output-max results/duohmm-$name \
	--output-graph results/duohmm-$name.graph

# correct haplotypes with duohmm
#duohmm \
#	-H results/duohmm-$name \
#	-M $MAP \
#	-O results/duohmm-$name-corrected

seq 1 20 | xargs -P 21 -I{} bash -c  "shapeit -convert --input-graph results/duohmm-$name.graph --output-sample results/sim$name.{} --seed {} && duohmm -H results/sim-$name-{} -M $MAP -R results/sim$name.{}.rec"

mapavg.py results/sim${name}.*.rec > duohmm-${name}-recombinations.txt
