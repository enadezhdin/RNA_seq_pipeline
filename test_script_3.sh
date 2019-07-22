#!/bin/bash


#for seq in data/DRR0161{25..30};

for i in 'seq 1 2';
do
    for seq in /home/en282/RNA_Seq/JL9_R/FASTQ/*_R${i}.fastq;
    do
        samp=`basename ${seq}`

        echo "Processing sample ${samp}"
    done
#salmon quant -i athal_index -l A \
#         -1 ${fn}/${samp}_1.fastq.gz \
#         -2 ${fn}/${samp}_2.fastq.gz \
#         -p 8 -o quants/${samp}_quant
done
