#!/bin/bash

for line in /home/en282/salmon_mapping_results/*;
do

    #echo "${line}"
    cd /${line}
    pwd
    pwd >> /home/en282/output.txt
    grep -w compatible_fragment_ratio lib_format_counts.json >> /home/en282/output.txt

    cd ..

done
