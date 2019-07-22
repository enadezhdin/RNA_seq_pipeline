#File manipulations notebook

from glob import glob
import os
import os.path
import subprocess
import sys

strain='098'
form='LC'
read='R2'


path = '/home/en282/RNA_Seq/FASTQ_Generation'

files = glob(path+"/{0}_{1}_rep3_L00*/*_{2}_*.fastq".format(strain,form,read))

#strain = os.path.basename(file).split("_")[1]
print(len(files))
print(files)
txt_output=open(path+'/{0}_{1}_rep3_{2}.fastq'.format(strain,form,read), 'a')

for file in files:
   rc=subprocess.call(['cat', file], stdout=txt_output)
txt_output.close()
