#!/bin/bash

# rep 1
fq_dir=/Volumes/BiochemXsan/staff_groups/brownlab/NZGL01919M_no-index-mismatch/NZGL01919M/
r1_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  
r2_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

mkdir -p skewer_rep1_out/
for f in $(cat fastq.ls); do \
    echo "Trimming adapter from ${f}...";
    skewer-0.2.2-linux-x86_64 -N -n -x ${r1_adapter} -y ${r2_adapter} \
    ${fq_dir}${f}_R1_001.fastq.gz ${fq_dir}${f}_R2_001.fastq.gz -o skewer_out/${f};
    pigz skewer_rep1_out/${f}-trimmed-pair*.fastq;
done

cd skewer_rep1_out/
rename 's/C9T78ANXX-1919M-05-2-1_S17_L005-trimmed-pair/A2_rep1_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-06-2-1_S18_L005-trimmed-pair/A2_rep2_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-07-2-1_S19_L005-trimmed-pair/B2_rep1_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-08-2-1_S20_L005-trimmed-pair/B2_rep2_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-09-2-1_S21_L005-trimmed-pair/C2_rep1_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-10-2-1_S22_L005-trimmed-pair/C2_rep2_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-11-2-1_S23_L005-trimmed-pair/pUC57_rep1_/' *.fastq.gz
rename 's/C9T78ANXX-1919M-12-2-1_S24_L005-trimmed-pair/pUC57_rep2_/' *.fastq.gz


# rep 2
fq_dir=/Volumes/userdata/student_users/chunshenlim/virus/doc/hbv/rnaseq/rep2/
r1_adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  
r2_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

mkdir -p skewer_out/
for f in $(cat fastq_rep2.ls); do \
    echo "Trimming adapter from ${f}...";
    skewer-0.2.2-linux-x86_64 -N -n -x ${r1_adapter} -y ${r2_adapter} \
    ${fq_dir}${f}_R1_001.fastq.gz ${fq_dir}${f}_R2_001.fastq.gz -o skewer_out/${f};
    pigz skewer_out/${f}-trimmed-pair*.fastq;
done

cd skewer_out/
rename 's/CC9ARANXX-3621-01-02-01_S13_L007/A2_rep2_/' *.fastq.gz
rename 's/CC9ARANXX-3621-02-02-01_S14_L007/B2_rep2_/' *.fastq.gz
rename 's/CC9ARANXX-3621-04-02-01_S15_L007/D2_rep2_/' *.fastq.gz
rename 's/CC9ARANXX-3621-05-02-01_S16_L007/pUC57_rep2_/' *.fastq.gz
rename 's/CC9ARANXX-3621-13-02-01_S19_L007/C2_rep2_/' *.fastq.gz
rename 's/CC9ARANXX-3621-14-02-01_S20_L007/D3_rep1.1_/' *.fastq.gz
rename 's/CC9ARANXX-3621-15-02-01_S21_L007/D3_rep1.2_/' *.fastq.gz
