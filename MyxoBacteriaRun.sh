#!/bin/csh

#$ -M aramezan@nd.edu	 # Email address for job notification
#$ -m  abe		 # Send mail when job begins, ends and aborts
#$ -q  long 	 # Specify queue 
#$ -N  MyxoBActeria	# Specify job name

module load gcc
./a.out 
