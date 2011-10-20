parallel bwa aln -n $mismatches -I $PWD/refs/$ref cutadapt/{}"$settings" ">" sai/{}"$settings.$ref.$miss.sai" ::: $samples \
&& parallel bwa samse -n 1 $PWD/refs/$ref sai/{}"$settings.$ref.$miss.sai" cutadapt/{}"$settings" ">" bwaAlignments/{}"$settings.$ref.$miss.sam" ::: $samples \
&& parallel samtools view -b -S bwaAlignments/{}"$settings.$ref.$miss.sam" ">" bwaAlignments/{}"$settings.$ref.$miss.bam" ::: $samples \
&& parallel --jobs 3 samtools sort -m 137438953472 bwaAlignments/{}"$settings.$ref.$miss.bam" bwaAlignments/{}"$settings.$ref.$miss.sorted" ::: $samples \                                                                                                             
&& parallel samtools index bwaAlignments/{}"$settings.$ref.$miss.sorted.bam" ::: $samples