#Parameters
ALIGNERS:= novo
PARAMSETS:= loose tight
REFS:= refs
REFGENOMES:= hairpin hg19.ambig tRNAs
STRATEGIES:= all random



#Target filenames
FASTQ_FILES :=          $(addsuffix .fq,$(addprefix $(FASTQDIR)/,$(SAMPLES)))
COUNT_FILES :=		$(FASTQ_FILES:.fq=.cnt)
SAMS :=                 $(addsuffix .sam,$(SAMPLES))
tRNA_SAM_FILES:=	$(addprefix $(BAMDIR)/$(ALIGNERS)/tight/tRNAs/all/,$(SAMS))
SAM_FILES    :=		$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/tight/hg19.ambig/$(strat)/,$(SAMS)))))) $(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/loose/hairpin/$(strat)/,$(SAMS)))))) $(tRNA_SAM_FILES)


BAM_FILES:=             $(SAM_FILES:.sam=.bam)
BAI_FILES:=             $(BAM_FILES:.bam=.bam.bai)
RCS_FILES:=             $(addsuffix .rcs,$(addprefix $(RCSDIR)/,$(SAMPLES)))

#convenience
tRNA_FILES:=	       $(tRNA_SAM_FILES:.sam=.bam.bai)
DOWNSTREAM_TARGETS:=  $(SAM_FILES) $(BAM_FILES) $(BAI_FILES)

#targets
default: bai count notmirna
bai: $(BAI_FILES)
bam: $(BAM_FILES)
sam: $(SAM_FILES)
fastq: linkfastq $(FASTQ_FILES)
count: $(COUNT_FILES)
rcs: $(RCS_FILES)
tRNA: $(tRNA_FILES)
clean:
	rm -f $(DOWNSTREAM_TARGETS)

.PHONY : linkfastq clean solexa all bai sam bam fastq

define align
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam: $(FASTQDIR)/%.fq
	mkdir -p $(BAMDIR)/$(1)/$(2)/$(3)/$(4)
	$($(1)_$(2)) -r $(4) -f $$< -d $(REFS)/$(3).ndx > $$@    
endef

$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call align,$(aligner),$(paramSet),$(ref),$(strat)))))))

#we could limit this to just the ref but this is an easier copy-paste job from align
#-@ INT    number of sorting threads [1]
#         -m INT    max memory per thread; suffix K/M/G recognized [768M]
define sam2bam
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam_tmp:  $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam
	$(SAMTOOLS) view -b -S $$< -t $(REFS)/$(3).fa > $$@
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam:  $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam_tmp
	$(SAMTOOLS) sort -@ 12 -m 4G $$< $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/$$*
endef

$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call sam2bam,$(aligner),$(paramSet),$(ref),$(strat)))))))

#index
%.bam.bai: %.bam
	$(SAMTOOLS) index $<

$(RCSDIR)/%.rcs: $(FASTQDIR)/%.fq
	cat $< | fastx_trimmer -l 26 -Q 33 | fastx_collapser -Q 33 | fasta_formatter -t | perl -ne 'm/\d+\-(\d+)\t(\S+)/;print $$2."\t".$$1."\n";' > $@

%.cnt:%.fq
	../exe/fastq-grep -c '.*' $< > $@

#fasta, in case some other program wants that
$(FASTADIR)/%.fq:$(DECODEDIR)/%.fq
	fastq_to_fasta < $< > $@


#unique
define unique
 $(BAMDIR)/$(1)/$(2)/$(3)/unique/%.bam:  $(BAMDIR)/$(1)/$(2)/$(3)/all/%.bam
	mkdir -p $(BAMDIR)/$(1)/$(2)/$(3)/unique
	$(SAMTOOLS) view -b -q 1 $$<  > $$@
endef



notmirnafiles:=$(addsuffix .notmirna.sorted.txt,$(addprefix $(BAMDIR)/novo/tight/hg19.ambig/all/,$(SAMPLES)))
notmirna:$(notmirnafiles) $(REFS)/hsa.chr.gff
%.notmirna.txt:%.bam
	$(BEDTOOLS) intersect -v -abam $< -b $(REFS)/hsa.chr.gff | $(SAMTOOLS) view - | cut -f 10 | uniq > $@
%.notmirna.sorted.txt:%.notmirna.txt
	cat $< | sort -u -S 150G -T /nas/is1/leipzig/ > $@
