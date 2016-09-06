for chrom in $(seq 1 22); do

cat <<EOF > scripts/chr$chrom.runner.sh
#!/bin/bash

#SBATCH -e logs/runner.$chrom.%J.err
#SBATCH -o logs/runner.$chrom.%J.out
#SBATCH -J $chrom.run
#SBATCH -t 96:00:00

#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp

set -eo pipefail
bash sites.sh $chrom
bash duohmm.sh $chrom
python run-hmm.py $chrom

EOF

sbatch scripts/chr$chrom.runner.sh
done
