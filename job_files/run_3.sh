#!/bin/bash
#PBS -q short
#PBS -l nodes=1:ppn=28
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -N ethanediol_3
#PBS -o /beegfs/work/st_163811/TAMie-force-field/example/pure_ethanediol/sim_3/LOG_ethanediol_3
#PBS -l mem=3000mb

# Load standard enviroment

module purge
module load mpi/openmpi/3.1-gnu-9.2

# Specify job directory and input file

v_dir=/beegfs/work/st_163811/TAMie-force-field/example/pure_ethanediol/sim_3
v_input=lammps.input

cd $v_dir
echo Submitting LAMMPS file: $v_input
mpirun --bind-to core --map-by core -report-bindings lmp -i $v_input -var seed $RANDOM