#File manipulations notebook

from glob import glob
import os
import os.path
import subprocess

import sys

path = '/home/en282/RNA_Seq/JL9_R/FASTQ'
files = glob(path+"/021_*/*.fastq")
print (len(files))
#files
for file in files:
    strain = os.path.basename(file).split("_")[1]
    print(strain)




