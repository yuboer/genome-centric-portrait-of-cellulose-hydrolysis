#!/bin/bash
#SBATCH --job-name=“dbCAN2_7899_genome”
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=yuboer@northwestern.edu
#SBATCH --output=o.dbCAN
#SBATCH --error=e.dbCAN


module load hmmer/3.1b2


cd /projects/b1052/Wells_b1042/yuboer/ongoing/sea_strains_checkM_results/bins/sea_faa/


#find . -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm ../dbCAN-fam-HMMs.txt $line > ${line}.out; done

find . -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm /home/ywx1845/software/dbCAN_db/dbCAN-HMMdb-V8.txt $line > ${line}.out; done

find . -name "*.out.dm"|while read line ; do sh /home/ywx1845/software/dbCAN/hmmscan-parser.sh $line > ${line}.ps; done

find . -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list

mkdir ../dbCAN_annotation_results

mv *.out.dm.ps > ../dbCAN_annotation_results