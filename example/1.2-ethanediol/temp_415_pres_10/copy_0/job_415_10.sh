#!/bin/bash
#PBS -q long
#PBS -l nodes=1:ppn=28
#PBS -l walltime=10:00:00
#PBS -N 1.2-ethanediol_415_10
#PBS -o /home/st/st_st/st_ac137577/workspace/software/TAMie-force-field/example/1.2-ethanediol/temp_415_pres_10/copy_0/LOG.o 
#PBS -e /home/st/st_st/st_ac137577/workspace/software/TAMie-force-field/example/1.2-ethanediol/temp_415_pres_10/copy_0/LOG.e 
#PBS -l mem=3000mb

module purge
module load mpi/openmpi/3.1-gnu-9.2

# Define the main working path 
WORKING_PATH=/home/st/st_st/st_ac137577/workspace/software/TAMie-force-field/example/1.2-ethanediol/temp_415_pres_10/copy_0
cd $WORKING_PATH

echo "This is the working path: $WORKING_PATH"


# Define the names of each simulation step taken. The folder as well as the output files will be named like this


################################# 
#       00_em       #
#################################
echo ""
echo "Starting ensemble: 00_em"
echo ""

mkdir -p 00_em
cd 00_em

mpirun --bind-to core --map-by core -report-bindings lmp -i em.input

echo "Completed ensemble: 00_em"

cd ../
sleep 10



################################# 
#       01_npt       #
#################################
echo ""
echo "Starting ensemble: 01_npt"
echo ""

mkdir -p 01_npt
cd 01_npt

mpirun --bind-to core --map-by core -report-bindings lmp -i npt.input

echo "Completed ensemble: 01_npt"

cd ../
sleep 10



################################# 
#       02_npt       #
#################################
echo ""
echo "Starting ensemble: 02_npt"
echo ""

mkdir -p 02_npt
cd 02_npt

mpirun --bind-to core --map-by core -report-bindings lmp -i npt.input

echo "Completed ensemble: 02_npt"

cd ../
sleep 10




# End
echo "Ending. Job completed."