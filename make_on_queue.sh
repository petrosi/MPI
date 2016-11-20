#!/bin/bash

## Give the Job a descriptive name
#PBS -N makempi

## Output and error files
#PBS -o makempi.out
#PBS -e makempi.err

## How many machines should we get?
#PBS -l nodes=1

## Start 
## Run make in the src folder (modify properly)
cd ~/lab03/mpi/
make

