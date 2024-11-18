#!/bin/sh

fastqc ./inputs/*.fq.gz

#trim ir adapt remove
#The -q (or --quality-cutoff) parameter can be used to trim low-quality ends from reads
# -m   minimum length 
#--nextseq-trim=20 This feature removes the polyG tails that arise from lack of signal in NextSeq/NovaSeq technologies
# -j Demultiplexing unique dual indices
# -O Minimum adapter overlap

for i in ./inputs/*1.fq.gz
do
   R1=$i
   R2="./inputs/"$(basename $R1 _1.fq.gz)"_2.fq.gz"
   base=$(basename $i _1.fq.gz)
   cutadapt --nextseq-trim=20 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -j 8 \
	-o ./trim/${base}_1.fq.gz \
	-p ./trim/${base}_2.fq.gz \
	${R1} ${R2} &> errors_${base}.log
done

fastqc ./trim/*.fq.gz

#a lot of ggggg at *_2.fq.gz, when removed with --nextseq-trim=20, Per base sequence content at the end: G drops (logical, removed all G repeant from seq ends)

#merge batches

for i in ./trim/*_A_1.fq.gz
do
	base=$(basename $i _A_1.fq.gz)
	zcat ./trim/${base}_A_1.fq.gz ./trim/${base}_B_1.fq.gz ./trim/${base}_C_1.fq.gz | gzip -c > ./trim_batch/${base}_1.fq.gz
	zcat ./trim/${base}_A_2.fq.gz ./trim/${base}_B_2.fq.gz ./trim/${base}_C_2.fq.gz | gzip -c > ./trim_batch/${base}_2.fq.gz
done

fastqc ./trim_batch/*.fq.gz

#Map, view map stats, sort
ref="ref/HV_M3_genome.fa"
for i in ./trim_batch/*_1.fq.gz
do
	base=$(basename $i _1.fq.gz)
	bwa mem -t 10 ${ref} ./trim_batch/${base}_1.fq.gz ./trim_batch/${base}_2.fq.gz 2> ./map/${base}_bwa_log.txt | samtools fixmate -m -@ 4 - tmp.bam
	samtools stats -@ 6 -in tmp.bam > ./map/${base}_bam.txt
	plot-bamstats ./map/${base}_bam.txt -p ./map/plots/${base}_bam
	samtools sort -@ 12 -o ./map/${base}_sorted.bam tmp.bam
done

#index bam using .csi index (.bai dosn't work on big genomes)
for i in ./map/*_sorted.bam
do
	samtools index -c $i
done

#call joined snp (WT_individ + mutant_pool) for mutmap analysis
for j in ./map/*_sorted.bam
do
	i=1
	while [ $i -le 7 ]
	do 
	base=$(basename $j _sorted.bam)
	bcftools mpileup --threads 8 -r ${i}H -a AD,ADF,ADR -Ou -f ./ref/HV_M3_genome.fa ./WT/AII_individ_sorted.bam $j | \
	bcftools call --threads 8 -vm -f GQ,GP -Oz -o ./Mutmap_snp/${base}_SNP_chr$i.vcf.gz
	i=$(($i+1))
	done
done

#plot Mutmap
for i in ./Mutmap_snp/*.vcf.gz
do 
	base=$(basename $i .vcf.gz)
	echo "mutplot -v $i -o ./Mutmap_results/${base} -n 50"
done
