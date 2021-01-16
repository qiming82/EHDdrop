#!/bin/sh
######################################################
#
# EXAMPLE MPICH SCRIPT FOR SGE
# Modified by Basement Supercomputing 1/2/2006 DJE
# To use, change "MPICH_JOB", "NUMBER_OF_CPUS"
# and "MPICH_PROGRAM_NAME" to real values.
#
# Modified for NJIT - AM 18Apr06
#
######################################################
#######################
# Job Name
#######################
#$ -N edrop0
#
###############################
# Use current working directory
###############################
#$ -cwd
#
###############################
# Join stdout and stderr
###############################
#$ -j y
#
###########################################################
# pe request for MPICH. Set your number of processors here.
# Make sure you use the "mpich" parallel environemnt.
###########################################################
# $ -pe mpich 1
#
###########################################################
# Run job through bash shell
###########################################################
# $ -S /bin/bash
#

##############################################################
# Full command path of your program
# Use full pathname to make sure we are using the right mpirun
##############################################################
/home/q/qw6/fwork/Edrop/edropo/e_vdrop