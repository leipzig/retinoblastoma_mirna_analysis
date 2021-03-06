#Directories
TOP := $(shell pwd)
FASTQDIR := $(TOP)/fastq
BAMDIR := $(TOP)/bam
SCRIPTDIR := $(TOP)/src
RCSDIR:= $(TOP)/rcs

#Programs
SAM_POSS:= /share/apps/bin/samtools /usr/bin/samtools
SAMTOOLS:= $(wildcard $(SAM_POSS))

NOVO_POSS:= /share/apps/bin/novoalign /usr/bin/novoalign
NOVOALIGN:= $(wildcard $(NOVO_POSS))

BEDTOOLS:= $(TOP)/../exe/bedtools

#Parameters
ALIGNERS:= novo
PARAMSETS:= loose tight
REFS:= $(TOP)/refs
REFGENOMES:= hairpin hg19 hg19.ambig
STRATEGIES:= all none random


#let's remove -m for variome crashes
#assume this is sanger-graded
novo_loose := $(NOVOALIGN)  -l 17 -h 60 -t 65 -o sam -o FullNW 
novo_tight := $(NOVOALIGN)  -l 17 -h 0 -t 0 -o sam -o FullNW

SAMPLES := RB525T01 RB525T02 RB525T03 WERI01 WERI02 WERI03

#Target filenames
FASTQ_FILES :=          $(addsuffix .fq,$(addprefix $(FASTQDIR)/,$(SAMPLES)))
COUNT_FILES :=		$(FASTQ_FILES:.fq=.cnt)
SAMS :=                 $(addsuffix .sam,$(SAMPLES))
BAMS :=	                $(addsuffix .bam,$(SAMPLES))
SAM_FILES    :=		$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/$(paramSet)/$(ref)/$(strat)/,$(SAMS))))))
BAM_FILES:=             $(SAM_FILES:.sam=.bam)
UNIQUE_FILES:=		$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/$(paramSet)/$(ref)/unique/,$(BAMS)))))
BAI_FILES:=             $(BAM_FILES:.bam=.bam.bai) $(UNIQUE_FILES:.bam=.bam.bai)

RCS_FILES:=             $(addsuffix .rcs,$(addprefix $(RCSDIR)/,$(SAMPLES)))



DOWNSTREAM_TARGETS:=  $(SAM_FILES) $(BAM_FILES) $(BAI_FILES)

#targets
default: bai rcs count notmirna
bai: $(BAI_FILES)
bam: $(BAM_FILES)
sam: $(SAM_FILES)
fastq: $(FASTQ_FILES)
count: $(COUNT_FILES)
rcs: $(RCS_FILES)
unique:$(UNIQUE_FILES)
clean:
	rm -f $(DOWNSTREAM_TARGETS)

.PHONY : clean solexa all bai sam bam fastq

define align
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam: $(FASTQDIR)/%.fq
	mkdir -p $(BAMDIR)/$(1)/$(2)/$(3)/$(4)
	$($(1)_$(2)) -r $(4) -f $$< -d $(REFS)/$(3).ndx > $$@    
endef

$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call align,$(aligner),$(paramSet),$(ref),$(strat)))))))

#we could limit this to just the ref but this is an easier copy-paste job from align
define sam2bam
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam_tmp:  $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam
	$(SAMTOOLS) view -b -S $$< -t $(REFS)/$(3).fa > $$@
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam:  $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam_tmp
	$(SAMTOOLS) sort $$< $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/$$*
endef

$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call sam2bam,$(aligner),$(paramSet),$(ref),$(strat)))))))

#unique
define unique
 $(BAMDIR)/$(1)/$(2)/$(3)/unique/%.bam:  $(BAMDIR)/$(1)/$(2)/$(3)/all/%.bam
	mkdir -p $(BAMDIR)/$(1)/$(2)/$(3)/unique
	$(SAMTOOLS) view -b -q 1 $$<  > $$@
endef

$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call unique,$(aligner),$(paramSet),$(ref))))))

#index
%.bam.bai: %.bam
	$(SAMTOOLS) index $<

$(RCSDIR)/%.rcs: $(FASTQDIR)/%.fq
	cat $< | fastx_trimmer -l 26 -Q 33 | fastx_collapser -Q 33 | fasta_formatter -t | perl -ne 'm/\d+\-(\d+)\t(\S+)/;print $$2."\t".$$1."\n";' > $@


%.cnt:%.fq
	../exe/fastq-grep -c '.*' $< > $@

notmirnafiles:=$(addsuffix .notmirna.txt,$(addprefix $(BAMDIR)/novo/tight/hg19.ambig/all/,$(SAMPLES)))
notmirna:$(notmirnafiles)
%.notmirna.txt:%.bam
	$(BEDTOOLS) intersect -v -abam $< -b $(REFS)/hsa.chr.gff | samtools view - | cut -f 10 | uniq > $@
