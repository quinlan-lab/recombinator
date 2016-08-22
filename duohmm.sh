VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
prefix=simons
PED=~u6000771/Data/ssc_519.ordered.ped
export VCF
export prefix
mkdir -p data/

#seq 1 22 | xargs -P 22 -I {} sh -c "bcftools view -m2 -M2 $VCF -c 1:nref -O b -o data/simons.chr{}.bcf {}"
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

vcf=data/$prefix.chr${chrom}.vcf.gz

vcf=illumina-sites.chr${chrom}.vcf.gz

vcf=data/$prefix.chr${chrom}.filter.vcf.gz
bcf=data/$prefix.chr${chrom}.bcf
bcf=illumina-sites.bcf

vcf=data/$prefix.sites.chrom.vcf.gz
bcftools view  \
	-e  'AVG(FMT/DP)<20 | MIN(FMT/DP)<6 | AVG(FMT/DP)>400 | AVG(FMT/GQ)<40 | MAF[*] < 0.03 | QUAL < 60'  \
	-O z -o $vcf $bcf $chrom


export name=chr${chrom}.${prefix}
D=results/2016_08_12-shapeit-LCR/
D=results/2016_08_21-shapeit-sites/
mkdir -p $D

<<DONE
bcftools view -f .,PASS \
	-e  'AVG(FMT/DP)<20 | MIN(FMT/DP)<6 | AVG(FMT/DP)>400 | AVG(FMT/GQ)<40 | MAF[*] < 0.03 | QUAL < 60'  \
	$bcf | bedtools intersect -sorted -header -a - -b <(zgrep -w ^$chrom data/LCR-hs37d5.bed.gz) -v \
	| bgzip -c > $vcf
exit;

# NOTE: add bcftools view -c 100 to force a higher AF
#(bcftools view -h $VCF;
#zgrep chr22 illumina-1M-duo.bed.gz \
#	| sort -k2,2n \
#	| bedtools merge -i stdin \
#	| sed -s 's/^chr//' | awk '{ print $1":"$2-1"-"$2+1 }' \
#	| gargs -p 20 "bcftools view -H $VCF {}" \
#	| awk 'length($5) == 1 && length($4) == 1' \
#	| awk 'BEGIN{x=0;} $0 ~/^#/{ if(x==0) {print;} next}{x=1; print $0 | "LC_ALL=C sort -u --compress-program gzip --buffer-size 1G -k1,1 -k2,2n"}'
#	) | bcftools view -m2 -M2 -c1 | bgzip -c > $vcf
DONE

plink --keep-allele-order --chr $chrom --geno 0.05 --mind 0.05 --vcf-half-call m --biallelic-only --vcf $vcf --make-bed --out $D/$name

cut -f 1-6 $PED | grep -v ^# > $D/${name}.fam

shapeit \
	--aligned \
	--duohmm \
	-B $D/$name \
	-T 33 \
	-M $MAP \
	-W 1 \
	-O $D/duohmm-$name

shapeit -convert --input-haps $D/duohmm-$name --output-vcf $D/$name.phased.vcf

bgzip -f $D/$name.phased.vcf
tabix $D/$name.phased.vcf.gz

#bash run.sh 

python recombinator.py \
    --min-gq 20 \
    --min-depth 20 \
    --region 22 \
    --vcf $D/$name.phased.vcf.gz \
    --ped $PED \
    --prefix $D/recomb
