#!/bin/bash

for seq in /home/en282/RNA_Seq/JL9_R/FASTQ/021_*.fastq;

do

samp=`basename ${seq}`

echo "${samp}" | awk -F'[_.]' '{print $2}'
    

    #for i in 'seq 1 2';
    #do
    


        #echo "Processing sample ${samp}"
   # done

done
