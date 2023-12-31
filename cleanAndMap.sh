#! /bin/bash 

#clean reads with process_radtags 
cat lists/names.txt |parallel -j 20 "process_radtags -1 raw/{}_R1.fastq.gz -2 raw/{}_R2.fastq.gz -o clean --renz_1 sphI --renz_2 mluCI -c -q"

#small cleaning step
cd clean
rm *rem* process*
cd ../

#map reads against chrysemys reference 
cat lists/names.txt |parallel -j 20 "bwa mem reference/PaintedTurtle.dna.toplevel.fa clean/{}.1.fastq.gz clean/{}.2.fastq.gz |samtools view -b -F 4 |samtools sort > bams/{}.bam"

