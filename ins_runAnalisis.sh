#for i in $(ls -d [2-3]*/*.gz); do echo $i;minimap2 -ax map-ont ./h37Rv/MtubH37Rv.fasta $i > $i.vsH37rv.sam;done;
for i in $(ls -d [2-3]*/*fastq.gz); do echo $i;minimap2 -ax map-ont ./embl_db/Mtuberculosis_ASM19595v2.fa $i > $i.vsH37rv.sam;done;
for i in $(ls -d [2-3]*/*.sam); do echo $i;samtools view -bS $i > $i.bam ;done;
for i in $(ls -d [2-3]*/*sam.bam); do echo $i;samtools sort -o $i.sort.bam $i;done
for i in $(ls -d [2-3]*/*sort.bam); do echo $i;samtools index $i ;done

eval "$(conda shell.bash hook)"
conda activate nanocaller_env
for i in $(ls -d [2-3]*/*sort.bam); do echo $i;python ~/software/NanoCaller_old/scripts/NanoCaller.py -bam $i -ref ./embl_db/Mtuberculosis_ASM19595v2.fa -chrom Chromosome -prefix varcall.$i; done;
conda deactivate

for i in $(ls -d [2-3]*/*varcal.snps.vcf); do echo $i; java -Xmx8g -jar ~/software/snpEff/snpEff.jar Mycobacterium_tuberculosis_h37rv $i > $i.ann.vcf
#java -Xmx8g -jar ~/software/snpEff/snpEff.jar Mycobacterium_tuberculosis_h37rv 3861_varcal.snps.vcf > 3861_varcal.snps.ann.vcf
#java -jar snpEff.jar dump Mycobacterium_tuberculosis_h37rv > check_test
