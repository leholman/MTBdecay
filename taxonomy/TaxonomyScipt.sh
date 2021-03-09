
TURN THIS INTO A TXT FILE AND UPLOAD IT
###HEADER FILE#####

#!/bin/bash

module load biobuilds
cd ~/scratch/nt_database

##################

MAKE A RESULTS FOLDER

#split fasta

awk -v size=100 -v pre=split -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' dada2OTUs.fasta

#write script for each file 

for i in split.* ;do echo -e "$(cat header.txt) \nblastn -query ../YiMTB/$i -db nt -out ../YiMTB/results/$i.Results.txt -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200 -num_threads 1" > $i.script.sh;done

#submit jobs

for f in *.sh; do echo qsub -l nodes=1:ppn=1 $f | bash; done

## Join files toghetehr - go to results folder

cat * > everyhting.txt


