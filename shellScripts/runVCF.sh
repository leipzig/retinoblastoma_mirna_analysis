#echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hsa_hairpin.dna.fa /storage/Ganguly/gangulyRBhi/RB494N.hairpin.bam /storage/Ganguly/gangulyRBhi/RB495N.hairpin.bam /storage/Ganguly/gangulyRBhi/RB498N.hairpin.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/normal.hairpin.bcf" | qsub
echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hg19.fa /storage/Ganguly/gangulyRBhi/RB494N.hg19.all.bam /storage/Ganguly/gangulyRBhi/RB495N.hg19.all.bam /storage/Ganguly/gangulyRBhi/RB498N.hg19.all.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/normal.hg19.all.bcf" | qsub
#echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hg19.fa /storage/Ganguly/gangulyRBhi/RB494N.hg19.bam /storage/Ganguly/gangulyRBhi/RB495N.hg19.bam /storage/Ganguly/gangulyRBhi/RB498N.hg19.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/normal.hg19.bcf" | qsub
#echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hsa_hairpin.dna.fa /storage/Ganguly/gangulyRBhi/RB494T.hairpin.bam /storage/Ganguly/gangulyRBhi/RB495T.hairpin.bam /storage/Ganguly/gangulyRBhi/RB498T.hairpin.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/tumor.hairpin.bcf" | qsub
echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hg19.fa /storage/Ganguly/gangulyRBhi/RB494T.hg19.all.bam /storage/Ganguly/gangulyRBhi/RB495T.hg19.all.bam /storage/Ganguly/gangulyRBhi/RB498T.hg19.all.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/tumor.hg19.all.bcf" | qsub
#echo "samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hg19.fa /storage/Ganguly/gangulyRBhi/RB494T.hg19.bam /storage/Ganguly/gangulyRBhi/RB495T.hg19.bam /storage/Ganguly/gangulyRBhi/RB498T.hg19.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/tumor.hg19.bcf" | qsub
