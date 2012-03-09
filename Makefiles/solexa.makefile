#Directories
TOP := $(shell pwd)
FASTQDIR := $(TOP)/fastq
BAMDIR := $(TOP)/bam
SCRIPTDIR := $(TOP)/src
RCSDIR:= $(TOP)/rcs

SOURCEDIR := $(TOP)/sources
DECODEDDIR := $(TOP)/decoded
FASTADIR := $(TOP)/fasta
TRIMMEDDIR := $(TOP)/trimmed

#Programs
SAMTOOLS:= /share/apps/bin/samtools
NOVOALIGN:= /share/apps/bin/novoalign

#Parameters
ALIGNERS:= novo
PARAMSETS:= loose tight
REFS:= $(TOP)/refs
REFGENOMES:= hairpin hg19 hg19.ambig
STRATEGIES:= all none random

MIN_LENGTH := 15

#assume this is sanger-graded
novo_loose := $(NOVOALIGN)  -l 17 -h 60 -t 65 -o sam -o FullNW 
novo_tight := $(NOVOALIGN)  -l 17 -h 0 -t 0 -o sam -o FullNW

#Sample info, some emerge from debarcoding, some are ready to go
#samples with barcodes
BC_SAMPLES := FGC0036_s_2.WERI FGC0031_s_8.WERI
#the un-barcoded sequences in a lane shared with barcoded sequences
UNBC_SAMPLES := FGC0036_s_2.Y79
#samples in lanes where no barcode was used
NONBC_SAMPLES := FGC0036_s_1.Y79 FGC0042_s_1.normal
SAMPLES := $(BC_SAMPLES) $(UNBC_SAMPLES) $(NONBC_SAMPLES)

BARCODED_LANES:= FGC0036_s_2 FGC0031_s_8

NONBARCODED_LANES:= FGC0036_s_1 FGC0042_s_1

#Target filenames
SOLEXA_FILES :=         $(wildcard $(SOURCEDIR)/*.txt.gz)
UNCOMPRESSED_FILES :=   $(SOLEXA_FILES:.txt.gz=.txt)
FASTQ_FILES :=          $(addprefix $(FASTQDIR)/,$(notdir $(UNCOMPRESSED_FILES:.txt=.fq)))
TRIMMED_FILES :=        $(addprefix $(TRIMMEDDIR)/,$(notdir $(FASTQ_FILES:.fq=.$(MIN_LENGTH).fq)))
BARCODED_SEQUENCES:=    $(addsuffix .$(MIN_LENGTH).fq,$(addprefix $(TRIMMEDDIR)/,$(BARCODED_LANES)))
NONBARCODED_SEQUENCES:= $(addsuffix .$(MIN_LENGTH).fq,$(addprefix $(TRIMMEDDIR)/,$(NONBARCODED_LANES)))
DEBARCODED_FILES:=      $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(BC_SAMPLES)))
READY_FILES:=           $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(NONBC_SAMPLES)))
DECODED_FILES :=        $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(SAMPLES)))
COUNT_FILES :=		$(DECODED_FILES:.fq=.cnt)
SAMS :=                 $(addsuffix .sam,$(SAMPLES))
BAMS :=	                $(addsuffix .bam,$(SAMPLES))
SAM_FILES    :=		$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/$(paramSet)/$(ref)/$(strat)/,$(SAMS))))))
BAM_FILES:=             $(SAM_FILES:.sam=.bam)
UNIQUE_FILES:=		$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(addprefix $(BAMDIR)/$(aligner)/$(paramSet)/$(ref)/unique/,$(BAMS)))))

BAI_FILES:=             $(BAM_FILES:.bam=.bam.bai) $(UNIQUE_FILES:.bam=.bam.bai)

RCS_FILES:=             $(addsuffix .rcs,$(addprefix $(RCSDIR)/,$(SAMPLES)))

DOWNSTREAM_TARGETS:= $(FASTQ_FILES) $(TRIMMED_FILES) $(DECODED_FILES) $(SAM_FILES)

#targets
default: bai rcs count
bai: $(BAI_FILES)
bam: $(BAM_FILES)
sam: $(SAM_FILES)
decoded: $(DECODED_FILES)
trimmed: $(TRIMMED_FILES)
fastq: $(FASTQ_FILES)
uncompressed: $(UNCOMPRESSED_FILES)
count: $(COUNT_FILES)
rcs: $(RCS_FILES)
unique:$(UNIQUE_FILES)

solexa:
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0036_s_1_sequence.txt.gz $(SOURCEDIR)/FGC0036_s_1.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0031_s_8_sequence.txt.gz $(SOURCEDIR)/FGC0031_s_8.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0036_s_2_sequence.txt.gz $(SOURCEDIR)/FGC0036_s_2.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/589-Ganguly-Chao-Solexa_smRNA/basic/Solexa/FGC0042_s_1_sequence.txt.gz $(SOURCEDIR)/FGC0042_s_1.txt.gz

clean:
	rm -f $(DOWNSTREAM_TARGETS)

.PHONY : clean solexa all bai sam bam fastq trimmed decoded


#Rules
$(SOURCEDIR)/%.txt:$(SOURCEDIR)/%.txt.gz
	gunzip -c $< > $@

#convert to fastq
#convert to sanger phred scores
$(FASTQDIR)/%.fq:$(SOURCEDIR)/%.txt
	python $(SCRIPTDIR)/python/FGC2fastq.py < $< > $@_tmp
	python $(SCRIPTDIR)/python/Solexa2Sanger.py $@_tmp $@
	rm $@_tmp

#trim adapters
$(TRIMMEDDIR)/%.$(MIN_LENGTH).fq:$(FASTQDIR)/%.fq
	cutadapt -a TCGTATGCCGTCTTCTGCTTG -m $(MIN_LENGTH) $< > $@

#FGC0036_s_2 is WERI end-barcoded with GTCT and Y79 not barcoded
#FGC0031_s_8 is WERI end-barcoded with GTCT
$(DEBARCODED_FILES): $(DECODEDDIR)/%.WERI.fq : $(TRIMMEDDIR)/%.$(MIN_LENGTH).fq
	cat $< | fastx_barcode_splitter.pl --bcfile barcodes.txt --eol --prefix $*. --suffix .decoded.tmp
	fastx_trimmer -Q33 -t 4 < $*.WERI.decoded.tmp > $@
	rm $*.WERI.decoded.tmp

#unmatched seqs are actually Y79
$(DECODEDDIR)/FGC0036_s_2.Y79.fq : $(DECODEDDIR)/FGC0036_s_2.WERI.fq
	cp -u FGC0036_s_2.unmatched.decoded.tmp $@

#FGC0036_s_1 is Y79 not barcoded
#FGC0042_s_1 is normal not barcoded
$(READY_FILES): $(NONBARCODED_SEQUENCES)
	cp $< $@

#fasta, in case some other program wants that
$(FASTADIR)/%.fq:$(DECODEDIR)/%.fq
	fastq_to_fasta < $< > $@

define align
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam: $(DECODEDDIR)/%.fq
	mkdir -p $(BAMDIR)/$(1)/$(2)/$(3)/$(4)
	$($(1)_$(2)) -r $(4) -f $$< -d $(REFS)/$(3).ndx > $$@    
endef

$(foreach strat,$(STRATEGIES),$(foreach ref,$(REFGENOMES),$(foreach paramSet,$(PARAMSETS),$(foreach aligner,$(ALIGNERS),$(eval $(call align,$(aligner),$(paramSet),$(ref),$(strat)))))))

#we could limit this to just the ref but this is an easier copy-paste job from align
define sam2bam
 $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.bam:  $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/%.sam
	$(SAMTOOLS) view -b -S $$< -t $(REFS)/$(3).fa > $$@_tmp
	$(SAMTOOLS) sort $$@_tmp $(BAMDIR)/$(1)/$(2)/$(3)/$(4)/$$*
	rm $$@_tmp
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


