VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=simons
PED=~u6000771/Data/ssc_519.ordered.ped
export VCF
export prefix
mkdir -p data/

set -e -o pipefail
mkdir -p results/

module load shapeit/v2.r837
module load duohmm/v0.1.7

# wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
export chrom=$1
export MAP=/uufs/chpc.utah.edu/common/home/u6000771/Projects/src/shapeit/example/genetic_map_b37/genetic_map_chr${chrom}_combined_b37.txt
# the hapmap ones have an extra column
#cut -f 2- /scratch/ucgd/lustre/u6000771/Projects/src/shapeit/example/genetic_map_GRCh37_chr${chrom}.txt > $MAP

vcf=data/$prefix.chr${chrom}.vcf.gz

vcf=data/$prefix.sites.chrom.vcf.gz

<<XXX
bcftools view  \
	-e  'AVG(FMT/DP)<20 | MIN(FMT/DP)<6 | AVG(FMT/DP)>400 | AVG(FMT/GQ)<40 | MAF[*] < 0.03 | QUAL < 60'  \
	-O z -o $vcf $bcf $chrom
XXX

vcf=illumina-omni2.5M-sites.${chrom}.vcf.gz


export name=chr${chrom}.${prefix}
D=results/2016_08_12-shapeit-LCR/
D=results/2016_08_29-shapeit-2.5M/
mkdir -p $D

plink --keep-allele-order --chr $chrom --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf $vcf --make-bed --out $D/$name

cut -f 1-6 $PED | grep -v ^# > $D/${name}.fam

shapeit \
	--aligned \
	--duohmm \
	--run 3 \
	-B $D/$name \
	-T 33 \
	-M $MAP \
	-W 2 \
	-O $D/duohmm-$name

shapeit -convert --input-haps $D/duohmm-$name --output-vcf $D/$name.phased.vcf

bgzip -f $D/$name.phased.vcf
tabix $D/$name.phased.vcf.gz

python recombinator.py \
    --min-gq 20 \
    --min-depth 20 \
    --region $chrom \
    --vcf $D/$name.phased.vcf.gz \
    --ped $PED \
    --prefix $D/recomb
