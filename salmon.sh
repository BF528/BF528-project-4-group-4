#!/usr/bin/bash

#$ -pe mpi_16_tasks_per_node 16

path="/projectnb/bf528/project_4_scrnaseq/fastq"
cwd="/projectnb/bf528/users/group4/project4/code"

salmon alevin -l ISR -1 $path/SRR3879606/SRR3879606_1_bc.fastq.gz -2 $path/SRR3879606/SRR3879606_2.fastq.gz --tgMap $cwd/transcript.tsv --whitelist $cwd/SRR3879606_whitelist.txt --end 5 --barcodeLength 19 --umiLength 6 --dumpMtx -i $cwd/index -p 16 -o salmon_06_out

