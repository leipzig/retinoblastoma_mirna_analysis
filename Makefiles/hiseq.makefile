#Directories
FASTQDIR := fastq
BAMDIR := bam
SCRIPTDIR := src
RCSDIR:= rcs

#Parameters
REFS:= refs
GENOME:= hairpin.ndx
TEMPLATE:= hairpin.dna.fa
#assume this is sanger-graded
NOVOCMD := novoalign -m -l 17 -h 60 -t 65 -o sam -o FullNW 

SAMPLES := RB525T WERI01 WERI02

#Target filenames
FASTQ_FILES :=          $(addsuffix .fq,$(addprefix $(FASTQDIR)/,$(SAMPLES)))
SAM_FILES :=            $(addsuffix .sam,$(addprefix $(BAMDIR)/,$(SAMPLES)))
BAM_FILES:=             $(SAM_FILES:.sam=.bam)
BAI_FILES:=             $(BAM_FILES:.bam=.bam.bai)
RCS_FILES:=             $(addsuffix .rcs,$(addprefix $(RCSDIR)/,$(SAMPLES)))



DOWNSTREAM_TARGETS:=  $(SAM_FILES) $(BAM_FILES) $(BAI_FILES)

#targets
default: bai
bai: $(BAI_FILES)
bam: $(BAM_FILES)
sam: $(SAM_FILES)
fastq: $(FASTQ_FILES)
rcs: $(RCS_FILES)

clean:
	rm -f $(DOWNSTREAM_TARGETS)

.PHONY : clean solexa all bai sam bam fastq

#align
$(BAMDIR)/%.sam: $(FASTQDIR)/%.fq
	$(NOVOCMD) -f $< -d $(REFS)/$(GENOME)  > $@

#sam2bam
$(BAMDIR)/%.bam: $(BAMDIR)/%.sam
	samtools view -b -S $< -t $(REFS)/$(TEMPLATE) > $@_tmp
	samtools sort $@_tmp $(BAMDIR)/$*
	rm $@_tmp

#index
$(BAMDIR)/%.bam.bai: $(BAMDIR)/%.bam
	samtools index $<

$(RCSDIR)/%.rcs: $(FASTQDIR)/%.fq
	cat $< | fastx_trimmer -l 26 -Q 33 | fastx_collapser -Q 33 | fasta_formatter -t | perl -ne 'm/\d+\-(\d+)\t(\S+)/;print $$2."\t".$$1."\n";' > $@
