#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-72:00:00
#SBATCH --output=myFBC_H7_150.stdout
#SBATCH --job-name="FBC_H7_150"
#SBATCH -p batch
#SBATCH -A alberlab
#SBATCH -p batch

echo "This assigns 12 threads to this job, creates all output folders, and runs a 2X P3D simulation".

export OMP_NUM_THREADS=12
mkdir yourpath/Animation_Cytoplasm
mkdir yourpath/Node_Data
mkdir yourpath/Nematic_Data
mkdir yourpath/Animation_No_Cytoplasm
mkdir yourpath/Cell_Data
mkdir yourpath/Tissue_Data
./program yourpath/Animation_Cytoplasm yourpath/Node_Data yourpath/Nematic_Data yourpath/Animation_No_Cytoplasm yourpath/Cell_Data yourpath/Tissue_Data -OOP_flag 0 -BFM 556.9484
