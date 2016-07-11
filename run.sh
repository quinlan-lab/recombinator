if [[ ! -e regions.txt ]]; then
python ~/Projects/src/freebayes/scripts/fasta_generate_regions.py \
	/scratch/ucgd/lustre/u0413537/UGP_Pipeline_Data/references/human_g1k_v37_decoy.fasta 4000000 \
	> regions.txt
fi

PATH=/scratch/ucgd/lustre/u6000771/gem/tools/bin:/scratch/ucgd/lustre/u6000771/gem/data/anaconda/bin/:$PATH:~u6000771/bin
VCF=results/chr22.simons/chr22.simons.phased.vcf.gz
PED=~u6000771/Data/ssc_519.ped
DATE=2016_07_09

mkdir -p scripts/
mkdir -p logs/
mkdir -p results/$DATE/

while read region; do
	f=${region//:/-}
	script=scripts/recombinator-$f.sh
cat << EOF > $script
#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp

#SBATCH --time=84:00:00
#SBATCH --ntasks=1
#SBATCH -J $f

#SBATCH -o logs/recombinator-$DATE-$f.out
#SBATCH -e logs/recombinator-$DATE-$f.err

set -eo pipefail -o nounset

python recombinator.py --min-gq 20 --min-depth 20 --region $region \\
	--vcf $VCF \\
	--ped $PED \
	--prefix results/$DATE/recomb &

EOF

for i in $(seq 1 23); do
	read region
	f=${region//:/-}
echo "\
python recombinator.py --min-gq 20 --min-depth 20 --region $region \\
	--vcf $VCF \\
	--ped $PED \
	--prefix  results/$DATE/recomb &
" >> $script
done

echo "wait" >> $script

sbatch $script

done <regions.txt
