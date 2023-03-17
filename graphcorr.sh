#!/bin/bash

#------- DEscription of the job -------

#SBATCH --job-name='pugno_gcorr'
#SBATCH --comment='Graphlet Correlation Project'

#------- Parameters -------

#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=bqhn01
#SBATCH --time=10-0:0:0
#SBATCH --mem=32000            # memory per node in MB 


#------- input/output -------

#SBATCH --output="/home/roberto_olayo/Graphlet_CORR_project/output.out"
#SBATCH --error="/home/roberto_olayo/Graphlet_CORR_project/output.err"

#------- Command -------
Rscript Script_simulation_for_cluster_HUB_Random.R