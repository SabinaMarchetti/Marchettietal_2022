#!/bin/bash
gcc main.c eq_.c -lgsl -lgslcblas -lm

Rscript ./RCode/0SetParams.R $2 $3
for ((i=0;i<=$1;i++))
    do
    echo ${i}
    ./a.out input_${i} > output_${i}
    Rscript ./RCode/1UpdateModel.R $i $1 $2 $3
done
