#!/usr/bin/bash

bed2vcf='./support/bedgovcf';
bgzip='./support/bgzip';
tabix='./support/tabix';
kanpig='./support/kanpig';
bedtools=''; # PLUG
samtools=''; # PLUG
bcftools=''; # PLUG
REF_GENOME='./support/GRCh38.p14.primary.fa';
WORK_DIR='tmp_dir'
RES_DIR='fc_kanpig_results'

mkdir -p ${WORK_DIR};
mkdir -p ${RES_DIR};
#rm -rf ${WORK_DIR}/*;
chmod a+w ${WORK_DIR};


# $1 - .bed file from FlyCatcher
# $2 - .bam file for said .bed file
# $3 - .yaml config for bedgovcf (subject to change)
pref=`basename "$1" | sed 's/.bed//g'`;


# Fetch reference sequence for future kanpig use the reference
paste $1 <(bedtools getfasta -fi $REF_GENOME -bed $1 -tab | cut -f 2 | tr '[:lower:]' '[:upper:]' | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, substr($1,1,1)}') | bedtools sort > ${WORK_DIR}/${pref}.fetched.bed;

# Convert .bed + sequences into a VCF. Subject to change if complications arise
$bed2vcf --bed ${WORK_DIR}/${pref}.fetched.bed --config $3 --fai ${REF_GENOME}.fai > ${WORK_DIR}/${pref}.vcf;
$bgzip ${WORK_DIR}/${pref}.vcf -o ${WORK_DIR}/${pref}.vcf.gz;
$tabix -f -p vcf ${WORK_DIR}/${pref}.vcf.gz;

# (Optional): do a kanpig plup for indexing if necessary
$kanpig plup --bam $2 | bedtools sort -header | $bgzip > ${WORK_DIR}/${pref}.plup.gz
$tabix -p ${WORK_DIR}/${pref}.plup.gz;


# Main: kanpig run and filter out bad variants
$kanpig gt --input ${WORK_DIR}/${pref}.vcf.gz \
           --reads $2 --reference ${REF_GENOME} --out ${WORK_DIR}/${pref}.kanpig.out.vcf;
bcftools view -e 'FMT/GT=="./." && FMT/FT==0 && FMT/SQ==0 && FMT/GQ==0 && FMT/PS=="." && FMT/NE==0 && FMT/DP==0 && FMT/AD=="." && FMT/KS=="."' ${WORK_DIR}/${pref}.kanpig.out.vcf -Ov -o ${RES_DIR}/${pref}.filtered.vcf

# clean-up
rm -f ${WORK_DIR}/${pref}.fetched.bed;
