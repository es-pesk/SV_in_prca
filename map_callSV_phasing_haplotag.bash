#!/bin/bash
source ~/.bashrc

#LN289253_chrX_asm.fasta

#singularity pull docker://mkolmogo/hapdup:0.12

file=`pwd`
#t2t="/home/es-pesk/prostate_cancer/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
hg38="/home/es-pesk/prostate_cancer/hg38.fa"
#file="/home/es-pesk/prostate_cancer/BAMs/full_genome"
asm="LN287457_asm_fg.fasta"
eg="LN287457"
#for eg in LN287460 LN289253; do
#for eg in `ls -1 *.bam | sed 's/\.bam//'`; do
    echo ${eg}
    # осознано просрали теги метилирования, вернемся к ним на этапе модкита
#samtools index -@ 10 ${eg}.bam
#samtools fasta -@ 10 ${file}/${eg}.bam > ${file}/${eg}.fa

    ### убираем дупликаты и сборка flye, если это нужно
#    if [[ "$eg" != "LN287458" && "$eg" != "LN287457" ]]; then
#samtools index ${eg}.bam
#samtools fasta ${file}/${eg}.bam > ${file}/${eg}.fa

#conda activate seqkit
#seqkit rmdup -n -o ${file}/x_${eg}.fa ${file}/${eg}.fa
#conda deactivate
#    fi
#conda activate flye
#flye --nano-raw ${file}/x_${eg}.fa \
# --keep-haplotypes \
# --threads 4 \
# --genome-size 3100m \
# --out-dir ${eg}_ \
# --asm-coverage 25
#conda deactivate
#mv ${file}/${eg}_/assembly.fasta ${file}/${eg}_asm.fasta


## до этого шага нужно получить сборку flye.

#conda activate minimap2_env
## картирование фазированных сборок #ref и сборка
#echo "minimap and samtools"
#set -e && samtools sort -@6 <(minimap2 -ax map-ont -t 10 ${hg38} ${file}/${eg}.fa) > ${file}/${eg}_minimap_sort_hg38.bam
#samtools index -@ 10 $file/${eg}_minimap_sort_hg38.bam
#echo "minimap and samtools done"
#conda deactivate
#done

## sniffles for sv -  меняем на 15
#echo "sniffles"
#conda activate sniffles2_env
#for eg in `ls -1 *_minimap_sort_hg38.bam| sed 's/\_minimap_sort_hg38.bam//'`;
#for eg in LN287457 LN287458 LN287459 LN287460 LN287833 LN289253;
#for eg in LN287460 LN289253;
#do #minsvlen cколько ставим? дефолт 50
#sniffles --input ${file}/${eg}_minimap_sort_hg38.bam --threads 10 --reference ${hg38} -v ${file}/${eg}_hg38_4svlen.vcf --minsvlen 4 --allow-overwrite
#done
#conda deactivate
#echo 'sniffles vcf done'


#eg='LN289253'
## склейка vcf - работает
echo 'vcf processing'
conda activate vacmap_env
#for eg in `ls -1 *_minimap_sort_hg38.bam| sed 's/\_minimap_sort_hg38.bam//'`;
#for eg in LN287458 LN287459 LN287833;
#do
#1.убираем аннотацию в снп и инелях
for eg in LN287457 LN287458 LN287459 LN287460 LN287833 LN289253; do
echo ${eg}
#for eg in LN287459 LN287460 LN289253;
#do
#bcftools annotate \
# --force \
# -x INFO,FORMAT/AD,FORMAT/DP,FORMAT/PL \
# -Oz \
# -o ${eg}_snp_noann.vcf.gz ${eg}.analysisReady.annotated_37_38_lift.vcf.gz
#2.убираем аннотацию в sv
#bgzip -f ${file}/${eg}_hg38_4svlen.vcf
#tabix -f ${file}/${eg}_hg38_4svlen.vcf.gz
#echo '#2.убираем аннотацию в sv' ##  этот шаг не используется
#bcftools annotate --force \
#-x INFO,FORMAT/AD,FORMAT/DP,FORMAT/PL \
#-Oz -o ${eg}_sv_noann_4svlen.vcf.gz ${eg}_hg38_4svlen.vcf.gz
##3. меняем название образца
#echo ${eg} > spisok.txt
#bcftools reheader --samples spisok.txt \
#-o ${eg}_hg38_renamed_4svlen.vcf.gz ${eg}_hg38_4svlen.vcf.gz
#tabix -f ${eg}_hg38_renamed_4svlen.vcf.gz

#bcftools reheader --samples spisok.txt \
#-o ${eg}_snp_noann_rename.vcf.gz ${eg}_snp_noann.vcf.gz
#tabix -f ${eg}_snp_noann_rename.vcf.gz

#5 конкатенация
bcftools concat -a -d both -Oz \
 -o ${eg}_sv_snp_4svlen.vcf.gz \
${eg}_snp_noann_rename.vcf.gz \
${eg}_hg38_renamed_4svlen.vcf.gz
tabix -f ${eg}_sv_snp_4svlen.vcf.gz

# сортировка после конкатенирования
bcftools sort -Oz -o ${eg}_sv_snp_sort_4svlen.vcf.gz ${eg}_sv_snp_4svlen.vcf.gz
tabix -f ${eg}_sv_snp_sort_4svlen.vcf.gz
echo 'done'
## фазирование вотсхэп - работает
echo 'whatshapp'
conda activate whatshap_env
whatshap phase --reference ${hg38} \
 --ignore-read-groups \
 --output ${file}/${eg}_hg38_phased_4svlen.vcf  \
${eg}_sv_snp_sort_4svlen.vcf.gz \
${file}/${eg}_minimap_sort_hg38.bam 2> ${eg}_whatshapp_4svlen.log
conda deactivate

#conda activate vacmap_env
bgzip -f ${eg}_hg38_phased_4svlen.vcf
tabix -f ${eg}_hg38_phased_4svlen.vcf.gz
echo 'done'

#оставляем только ps теги и удаляем hp. HP пришли из vcf  генотипов
echo 'rm HP-tag'
bcftools annotate -x FORMAT/HP -Oz -o ${eg}_hg38_phased_corr_4svlen.vcf.gz ${eg}_hg38_phased_4svlen.vcf.gz
tabix -f ${eg}_hg38_phased_corr_4svlen.vcf.gz
echo 'done'

#гаплотег
echo 'haplotag'
conda activate whatshap_env
whatshap haplotag -o ${eg}_haplotagged_4svlen.bam \
--reference ${hg38} \
--output-threads=6 \
--ignore-read-groups \
${file}/${eg}_hg38_phased_corr_4svlen.vcf.gz ${file}/${eg}_minimap_sort_hg38.bam 2> ${eg}_haplotag_4svlen.log

samtools index -@ 6 ${eg}_haplotagged_4svlen.bam
conda deactivate
done
echo 'done'

#modkit  - пример modkit pileup -r assembly.fasta -b phased_reads.bam -o methylation_phased.bed - так мы расфазируем метилирование как раз

#for eg in `ls -1 *_haplotagged.bam|sed 's/\_minimap_sort_hg38.bam//'`;
#do
#modkit pileup -r ${eg}_haplotagged.bam -b phased_reads.bam -o methylation_phased.bed #-t
#done
