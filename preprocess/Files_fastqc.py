#File manipulations notebook

from glob import glob
import os
import os.path
import subprocess


path = '/home/en282/RNA_Seq/FASTQ_Generation'

files = glob(path+"/*.fastq")
print(len(files))


for file in files:
   subprocess.call(['fastqc', file])
    #print(file)
