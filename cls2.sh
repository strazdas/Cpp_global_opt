#!/bin/bash
#
#$ -j y
#$ -S /bin/bash
#$ -cwd

mpirun $HOME/C++_global_opt/mpi_main.out 2
