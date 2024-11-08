#!/bin/sh
sbatch <<EOT
#!/bin/sh

#SBATCH --account=physics
#SBATCH --partition=ada
#BATCH --time=20:00:00
#SBATCH --nodes=1 --ntasks=15
#SBATCH --job-name="pythia_main09_20M_502_tune5_$1"
#SBATCH --mail-user=ptgjak001@myuct.ac.za
#SBATCH --mail-type=ALL

module load compilers/gcc820

source /scratch/ptgjak001/root/bin/thisroot.sh

./mymain09_HPC $1 > mymain09_HPC_root_20M_502_tune5/myout09_$1
EOT