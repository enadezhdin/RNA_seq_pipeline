#!/bin/bash
#SBATCH --job-name=B.subtilis_mapping_3_test       # the name of the job
#SBATCH -o B.subtilis_mapping_3_test.stdout        # the standard output file
#SBATCH -e B.subtilis_mapping_3_test.stderr        # the standard error file
#SBATCH -n 16                                      # the number of CPUs for the job
#SBATCH --mem=128G                                 # Allocated memory


#indexin_run
#salmon index -t /home/en282/salmon_tutorial/Bacillus_subtilis.cdna.all.fa.gz -i /home/en282/salmon_tutorial/Bacillus_subtilis_index --type quasi -k 31

#mapping_run

salmon quant -i /home/en282/salmon_tutorial/Bacillus_subtilis_index -l A -1 /home/en282/RNA_Seq/JL9_R/FASTQ/021_B_rep1_R1.fastq -2 /home/en282/RNA_Seq/JL9_R/FASTQ/021_B_rep1_R2.fastq -o /home/en282/salmon_mapping_results/B_subtilis_transcripts_quant_test_3 --seqBias --gcBias
