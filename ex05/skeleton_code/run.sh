#!/bin/bash

N=100

lscpu

rm -f out
mkdir -p out

for i in 1 2 4 8; do

  echo "$i threads run with $N"

  mpirun -n $i ./diffusion 1.0 2.0 $N # add --oversubscribe if you have more than less than 48 cores
done

