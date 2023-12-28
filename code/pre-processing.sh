#!/bin/sh

#mark duplicates

for i in ./BAM/*.bam
   do
      base=$(basename $i _sorted.bam)
   java -jar $PICARD SortSam \
      I=$i \
      O=tmp.bam \
      SORT_ORDER=coordinate
   java -jar $PICARD MarkDuplicates \
      I=tmp.bam \
      O=no_dup_bam/${base}_maked_dub.bam \
      M=no_dup_bam/${base}_marked_dup_metrics.txt
done

#reik specialaus ref indexo .dict
$GATK CreateSequenceDictionary -R ref/HV_M3_genome.fa

#Index known SNP file and fix errors od dublitated ref seq in alt column

cat ref/HV_SNP.vcf | sed -e 's/IBSC_0;//g' | sed -e 's/JHI_0;//g'> ref/HV_SNP_fixed.vcf
bgzip -c ref/HV_SNP_fixed.vcf > ref/HV_SNP_fixed.vcf.gz
tabix -p vcf -C ref/HV_SNP_fixed.vcf.gz
bcftools norm -a ref/HV_SNP_fixed.vcf.gz > ref/HV_SNP_fixed.vcf
$GATK IndexFeatureFile -I ref/HV_SNP_fixed.vcf

# fix the read groups

for i in ./no_dub_bam/*.bam
    do
    base=$(basename $i _maked_dup.bam)
    java -jar $PICARD AddOrReplaceReadGroups \
    I=$i \
    O=./fixed_group_bam/${base}_fixed_group.bam \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGSM=Sample1 \
    RGPU=unit1 \
    CREATE_INDEX=FALSE
done

##nex tike use this

# for i in ./recal_bam/*.bam
#     do
#     base=$(basename $i _recal.bam)
#     java -jar $PICARD AddOrReplaceReadGroups \
#     I=$i \
#     O=./fixed_group_bam/${base}_final.bam \
#     SORT_ORDER=coordinate \
#     RGID=foo \
#     RGLB=bar \
#     RGPL=illumina \
#     RGSM=Sample_${base} \
#     RGPU=unit_${base} \
#     CREATE_INDEX=FALSE
# done

#base recalibration
   
for i in ./fixed_group_bam/*.bam
   do
   base=$(basename $i _fixed_group.bam)
   $GATK BaseRecalibrator \
   -I $i \
   -R ref/HV_M3_genome.fa \
   --known-sites ref/HV_SNP_fixed.vcf\
   -O recal/${base}_recal.table
done

for i in ./fixed_group_bam/*.bam
   do
   base=$(basename $i _fixed_group.bam)
   $GATK ApplyBQSR \
   -R ref/HV_M3_genome.fa\
   -I $i \
   --bqsr-recal-file recal/${base}_recal.table \
   -O recal_bam/${base}_recal.bam \
   -OBI FALSE
done

#sutaisyt dar AII_individ ir AII_pool ir indeksuot .csi
 $GATK ApplyBQSR \
   -R ref/HV_M3_genome.fa\
   -I fixed_group_bam/AII_individ_fixed_group.bam \
   --bqsr-recal-file recal/AII_individ_recal.table \
   -O recal_bam/AII_individ_recal.bam \
   -OBI FALSE

for i in ./fixed_group_bam/*.bam
   do
   base=$(basename $i _fixed_group.bam)
   $GATK ApplyBQSR \
   -R ref/HV_M3_genome.fa\
   -I $i \
   --bqsr-recal-file recal/${base}_recal.table \
   -O recal_bam/${base}_recal.bam \
   -OBI FALSE
done

#fix sample names in bam, mistake was made erlier.

for i in ./recal_bam/*.bam
    do
    base=$(basename $i _recal.bam)
    java -jar $PICARD AddOrReplaceReadGroups \
    I=$i \
    O=./fixed_group_bam/${base}_final.bam \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGSM=Sample_${base} \
    RGPU=unit_${base} \
    CREATE_INDEX=FALSE
done

#komanda su=inoti bam sample name
$GATK GetSampleName \
     -I fixed_group_bam/AII_individ_final.bam \
     -O sample_name.txt

#index to use on igv and mutect2

for i in fixed_group_bam/*.bam
      do
      samtools index -c $i
      done
#use mutect2

# for i in fixed_group_bam/*.bam
#       do
#       base=$(basename $i _final.bam)
#       $GATK Mutect2 \
#      -R ref/HV_M3_genome.fa \
#      -I tumor.bam \
#      -I normal.bam \
#      -normal normal_sample_name \
#      --germline-resource af-only-gnomad.vcf.gz \
#      --panel-of-normals pon.vcf.gz \
#      -O somatic.vcf.gz


for i in fixed_group_bam/*.bam
       do
       base=$(basename $i _final.bam)
       $GATK Mutect2 \
       -R ref/HV_M3_genome.fa  \
       -I $i \
       -O Mutect2/${base}_mutect2.vcf.gz
       --create-output-variant-index FALSE
done



# #test with individs
# $GATK Mutect2 \
#      -R ref/HV_M3_genome.fa \
#      -I fixed_group_bam/tw2_individ_final.bam \
#      -I fixed_group_bam/AII_individ_final.bam \
#      -normal Sample_AII_individ \
#      -O tw2_in_AII_in_somatic.vcf.gz

# java -jar $PICARD ValidateSamFile \
#       I=tw_pool_fixed_group.bam\
#       MODE=SUMMARY



#    $GATK BaseRecalibrator \
#    -I tw_pool_fixed_group.bam \
#    -R ref/HV_M3_genome.fa \
#    --known-sites ref/HV_SNP_fixed.vcf\
#    -O recal/tw_pool_recal2.table

#  $GATK ApplyBQSR \
#    -R ref/HV_M3_genome.fa\
#    -I tw_pool_fixed_group.bam \
#    --bqsr-recal-file recal/tw_pool_recal2.table \
#    -O recal_bam/tw_pool_recal.bam \
#    -OBI FALSE


#     java -jar $PICARD AddOrReplaceReadGroups \
#     I=recal_bam/tw_pool_recal2.bam \
#     O=./fixed_group_bam/tw_pool_final.bam \
#     SORT_ORDER=coordinate \
#     RGID=foo \
#     RGLB=bar \
#     RGPL=illumina \
#     RGSM=Sample_tw_pool \ 
#     RGPU=unit_tw_pool \
#     CREATE_INDEX=FALSE

test