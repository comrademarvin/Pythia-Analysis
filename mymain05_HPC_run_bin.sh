#!/bin/sh
sbatch <<EOT
#!/bin/sh

#SBATCH --account=physics
#SBATCH --partition=ada
#BATCH --time=48:00:00
#SBATCH --nodes=1 --ntasks=15
#SBATCH --job-name="mymain05_HPC_50M_536_$1"
#SBATCH --mail-user=ptgjak001@myuct.ac.za
#SBATCH --mail-type=ALL

module load compilers/gcc820

source /scratch/ptgjak001/root/bin/thisroot.sh

./mymain05_HPC $1 > mymain05_HPC_root_50M_536/myout05_$1
EOT