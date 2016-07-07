VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=simons
PED=~u6000771/Data/ssc_519.ordered.ped
export VCF
export prefix
mkdir -p data/

#seq 1 22 | xargs -P 22 -I {} sh -c "bcftools view -m2 -M2 $VCF -c 1:nref -O z -o data/$prefix.chr{}.vcf.gz {}"
#exit



set -e -o pipefail
mkdir -p results/

module load shapeit/v2.r837
module load duohmm/v0.1.7

# wget http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
export chrom=22
export MAP=/uufs/chpc.utah.edu/common/home/u6000771/Projects/src/shapeit/example/genetic_map_b37/genetic_map_chr${chrom}_combined_b37.txt
# the hapmap ones have an extra column
#cut -f 2- /scratch/ucgd/lustre/u6000771/Projects/src/shapeit/example/genetic_map_GRCh37_chr${chrom}.txt > $MAP

export name=chr${chrom}.${prefix}.sites
vcf=data/$prefix.chr${chrom}.vcf.gz
vcf=illumina-sites.chr${chrom}.vcf.gz

(bcftools view -h $VCF;
zgrep chr22 illumina-1M-duo.bed.gz \
	| sort -k2,2n \
	| bedtools merge -i stdin \
	| sed -s 's/^chr//' | awk '{ print $1":"$2-1"-"$2+1 }' \
	| gargs -p 20 "bcftools view -H $VCF {}" \
	| awk 'length($5) == 1 && length($4) == 1' \
	| awk 'BEGIN{x=0;} $0 ~/^#/{ if(x==0) {print;} next}{x=1; print $0 | "LC_ALL=C sort -u --compress-program gzip --buffer-size 1G -k1,1 -k2,2n"}'
	) | bcftools view -m2 -M2 -c1 | bgzip -c > $vcf


echo "LINES:"
zgrep -cv ^# $vcf


plink --real-ref-alleles --chr $chrom --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf $vcf --make-bed --out $name

cut -f 1-6 $PED | grep -v ^# > ${name}.fam

rm -f $name-nomendel*
# remove mendelian errors
plink \
 --real-ref-alleles \
 --bfile $name \
 --me 1 1 \
 --set-me-missing \
 --make-bed \
 --out $name-nomendel

mkdir -p results/$name/
rm -f results/$name/*
mkdir -p results/$name/
# phase without using ped info.
shapeit \
	--aligned \
	--duohmm \
	-B $name-nomendel \
	-T 27 \
	-M $MAP \
	-W 0.5 \
	-O results/$name/wduohmm
	#--noped \

#duohmm --haps results/$name/duohmm -M $MAP -G results/$name/genotyping-errors.txt --output-hap results/$name/duohmm.corrected
shapeit -convert --aligned --input-haps results/$name/wduohmm --output-vcf results/$name/wduohmm.corrected.vcf

<<DUOHMM
seq 1 20 | xargs -P 21 -I{} bash -c  "shapeit -convert --input-graph results/$name/duohmm.graph --output-sample results/$name/sim$name.{} --seed {} && duohmm -H results/$name/sim$name.{} -M $MAP -R results/$name/sim$name.{}.rec"

mapavg.py results/$name/sim${name}.*.rec > duohmm-u-${name}-recombinations.txt
DUOHMM
