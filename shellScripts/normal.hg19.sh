samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hg19.fa /storage/Ganguly/gangulyRBhi/RB494N.hg19.bam /storage/Ganguly/gangulyRBhi/RB495N.hg19.bam /storage/Ganguly/gangulyRBhi/RB498N.hg19.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/normal.hg19.bcf