python ~/Projects/src/freebayes/scripts/fasta_generate_regions.py \
	/scratch/ucgd/lustre/u0413537/UGP_Pipeline_Data/references/human_g1k_v37_decoy.fasta 1000000 \
	> regions.txt

PATH=/scratch/ucgd/lustre/u6000771/gem/tools/bin:/scratch/ucgd/lustre/u6000771/gem/data/anaconda/bin/:$PATH:~u6000771/bin

mkdir -p scripts/
mkdir -p logs/
mkdir -p results/

while read region; do
	f=${region//:/-}
	script=scripts/$f.sh
cat << EOF > $script
#!/bin/bash
#SBATCH --account=ucgd-kp
#SBATCH --partition=ucgd-kp

#SBATCH --time=84:00:00
#SBATCH --ntasks=1
#SBATCH -J $f

#SBATCH -o load-$f.out
#SBATCH -e load-$f.err

set -eo pipefail -o nounset

python recombinator.py --min-gq 20 --min-depth 20 --region $region \\
	--vcf ~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz \\
	--ped ~u6000771/Data/ssc_519.ped \
	| bgzip -c > results/$f.bed.gz &

EOF

for i in $(seq 1 23); do
	read region
	f=${region//:/-}
echo "\
python recombinator.py --min-gq 20 --min-depth 20 --region $region \\
	--vcf ~u6000771/Data/519FamiliesUnrecal_snp-recal_indel-recal.vcf.gz \\
	--ped ~u6000771/Data/ssc_519.ped \
	| bgzip -c > results/$f.bed.bz &
" >> $script
done

echo "wait" >> $script

sbatch $script

done <regions.txt
