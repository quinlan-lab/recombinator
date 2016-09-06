chrom=$1
VCF=~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz
sites=illumina-1M-duo.bed.gz
sites=data/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
(bcftools view -h $VCF;
bedtools merge -i $sites \
   | grep -w ^$chrom \
   | grep -v random \
   | grep -v ^GL \
   | grep -v ^# \
   | sort -k1,1 -k2,2n \
   | sed -s 's/^chr//' | awk '{ print $1":"$2-1"-"$2+1 }' \
   | gargsx -n 100 -p 20 "bcftools view -H $VCF {}" \
   | awk 'length($5) <= 2 1 && length($4) <= 2' \
   | awk 'BEGIN{x=0;} $0 ~/^#/{ if(x==0) {print;} next}{x=1; print $0 | "LC_ALL=C sort -u --compress-program pigz --buffer-size 3G -k1,1 -k2,2n"}'
   ) \
   | bcftools view -m2 -M2 -c1 -O z -o illumina-omni2.5M-sites.$chrom.vcf.gz
