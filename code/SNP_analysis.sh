#!/bin/sh

wget -O ./ref/HV_M3_genome.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/plants/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna_sm.toplevel.fa.gz

gunzip ./ref/HV_M3_genome.fa.gz

bwa index ./ref/HV_M3_genome.fa

#komandos leistos iš AII_pool direktorijos
#QC
fastqc ./inputs/*.fq.gz

#trim ir adapt remove
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -u 10 -U 10 -j 8 \
	-o ./trim/AII_1.fq.gz \
	-p ./trim/AII_2.fq.gz \
	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_1.fq.gz \
	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_2.fq.gz &> errors_AII.log

#QC po trim
fastqc ./trim/*.fq.gz

#Map
bwa mem -t 6 ../ref/HV_genome.fa ./trim/AII_1.fq.gz ./trim/AII_2.fq.gz 2> ./map/AII_bwa_log.txt| samtools fixmate -m -@2 - ./map/AII.bam
samtools stats -in ./map/AII.bam > ./map/AII_bam.txt
plot-bamstats ./map/AII_bam.txt -p ./map/plots/AII_bam

#sort
samtools sort -@6 -o ./map/AII_pool_sorted.bam ./map/AII.bam
samtools stats -@ 6 -in ./map/AII_pool_sorted.bam > ./map/AII_bam_sorted.txt
plot-bamstats ./map/AII_bam_sorted.txt -p ./map/plots_sort/AII_bam

#index and extract non-rek region bam
samtools index -c ./map/AII_pool_sorted.bam
samtools view -@ 6 -b ./map/AII_pool_sorted.bam "5H:490000000-580000000" > ./map/AII_pool_5H.bam

i=1
while [ $i -le 7 ]
do 
	base=$(basename $i)
	bcftools mpileup --threads 4 -r ${base}H -a AD,ADF,ADR -Ou -f ./ref/HV_M3_genome.fa ./map/AII_pool_sorted.bam | \
	bcftools call --threads 4 -vm -f GQ,GP -Oz -o ./results/AII_pool_SNP_chr$i.vcf.gz
	i=$(($i+1))
done

#tw pool mut map
scp -o 'ProxyJump studremdesk@172.16.170.44' bioinformatikai@172.17.168.98:~/mieziai/tw_pool_2023/results/mutmap* .

i=1
while [ $i -le 7 ]
do 
	base=$(basename $i)
	mutplot -v mutmap_${base}H_tw_pool.vcf.gz -o ./mutmap_tw_pool_${base}H -n 50
	i=$(($i+1))
done

#tw2 mutmap su AII individ
scp -o 'ProxyJump studremdesk@172.16.170.44' bioinformatikai@172.17.168.99:~/mieziai/tw_pool/results/mutmap* .

i=1
while [ $i -le 7 ]
do 
	base=$(basename $i)
	mutplot -v mutmap_${base}H_tw2_pool.vcf.gz -o ./mutmap_tw2_pool_${base}H -n 50
	i=$(($i+1))
done

#tw2 pool mutmap
scp -o 'ProxyJump studremdesk@172.16.170.44' bioinformatikai@172.17.168.99:~/mieziai/tw_pool/results/mutmap* .

i=1
while [ $i -le 7 ]
do 
	base=$(basename $i)
	mutplot -v mutmap_${base}H_tw2_pool.vcf.gz -o ./mutmap_tw2_pool_${base}H -n 50
	i=$(($i+1))
done

