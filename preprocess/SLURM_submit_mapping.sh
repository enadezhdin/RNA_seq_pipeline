#!/bin/bash

#SBATCH --job-name=B.subtilis_mapping_all       # the name of the job
#SBATCH -o B.subtilis_mapping_all.stdout        # the standard output file
#SBATCH -e B.subtilis_mapping_all.stderr        # the standard error file
#SBATCH -n 32                                      # the number of CPUs for the job
#SBATCH --mem=128G                                 # Allocated memory

for seq in /home/en282/RNA_Seq/JL9_R/FASTQ/*.fastq;

do

samp=`basename ${seq}`


strain=$(echo "${samp}" | awk -F'[_.]' '{printf $1}')
ec_type=$(echo "${samp}" | awk -F'[_.]' '{printf $2}')
rep=$(echo "${samp}" | awk -F'[_.]' '{printf $3}')
R=$(echo "${samp}" | awk -F'[_.]' '{printf substr($4, 2)}')
 
    for str in ${strain};
    do
        for ec in ${ec_type};
        do
            for repl in ${rep};
            do
                for r in ${R};
                do
                      if [ ${r} -eq 1 ]; then
                          echo "Processing sample ${str} ${ec} ${repl} ${R}"
                          salmon quant -i /home/en282/salmon_tutorial/Bacillus_subtilis_index -l A -1 /home/en282/RNA_Seq/JL9_R/FASTQ/${str}_${ec}_${repl}_R1.fastq -2 /home/en282/RNA_Seq/JL9_R/FASTQ/${str}_${ec}_${repl}_R2.fastq -o /home/en282/salmon_mapping_results/B_subtilis_${str}_${ec}_${repl} --seqBias --gcBias

                      fi
                done


            done
        done
    done
done
