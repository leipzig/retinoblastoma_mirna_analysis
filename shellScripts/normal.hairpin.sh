samtools mpileup -uDf /storage/Ganguly/gangulyRBhi/refs/hsa_hairpin.dna.fa /storage/Ganguly/gangulyRBhi/RB494N.hairpin.bam /storage/Ganguly/gangulyRBhi/RB495N.hairpin.bam /storage/Ganguly/gangulyRBhi/RB498N.hairpin.bam | bcftools view -bvc - > /storage/Ganguly/gangulyRBhi/normal.hairpin.bcf