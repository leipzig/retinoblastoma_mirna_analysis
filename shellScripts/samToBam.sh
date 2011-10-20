for f in */*/*hg19*.sam; do echo "samtools view -b -S $PWD/$f -T $PWD/refs/hg19.fai > $PWD/$f.bam; samtools sort $PWD/$f.bam $PWD/$f.sorted; samtools index $PWD/$f.sorted.bam" | qsub; done
for f in */*/*hairpin.sam; do echo "samtools view -b -S $PWD/$f -t $PWD/refs/hairpin_dna.fa > $PWD/$f.bam; samtools sort $PWD/$f.bam $PWD/$f.sorted; samtools index $PWD/$f.sorted.bam" | qsub ; done



