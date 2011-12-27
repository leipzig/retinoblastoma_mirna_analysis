#Directories
SOURCEDIR := sources
FASTQDIR := fastq
DECODEDDIR := decoded
BAMDIR := bam
SCRIPTDIR := scripts
TRIMMEDDIR := trimmed

#Parameters
MIN_LENGTH := 15
REFS:= /home/leipzig/ganguly/gangulyRBhi/data/refs
GENOME:= hairpin.ndx
TEMPLATE:= hairpin.dna.fa
#assume this is sanger-graded
NOVOCMD := novoalign -m -l 17 -h 60 -t 65 -o sam -o FullNW 

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
SAM_FILES :=            $(addsuffix .sam,$(addprefix $(BAMDIR)/,$(SAMPLES)))
BAM_FILES:=             $(SAM_FILES:.sam=.bam)
BAI_FILES:=             $(BAM_FILES:.bam=.bam.bai)

DOWNSTREAM_TARGETS:= $(FASTQ_FILES) $(TRIMMED_FILES) $(DECODED_FILES) $(SAM_FILES)

#targets
default: bai
bai: $(BAI_FILES)
bam: $(BAM_FILES)
sam: $(SAM_FILES)
decoded: $(DECODED_FILES)
trimmed: $(TRIMMED_FILES)
fastq: $(FASTQ_FILES)
uncompressed: $(UNCOMPRESSED_FILES)


clean:
	rm -f $(DOWNSTREAM_TARGETS)

.PHONY : clean all bai sam bam fastq trimmed decoded

#Explicit rules prevent intermediates from being deleted
#$(SAM_FILES):$(DECODED_FILES)
#$(BAM_FILES):$(SAM_FILES)
#$(BAI_FILES):$(BAM_FILES)

#Rules
$(SOURCEDIR)/%.txt:$(SOURCEDIR)/%.txt.gz
	gunzip -c $< > $@

#convert to fastq
#convert to sanger phred scores
$(FASTQDIR)/%.fq:$(SOURCEDIR)/%_sequence.txt
	python $(SCRIPTDIR)/FGC2fastq.py $< > $@_tmp
	python $(SCRIPTDIR)/Solexa2Sanger.py $@_tmp $@
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
	mv FGC0036_s_2.unmatched.decoded.fq $@

#FGC0036_s_1 is Y79 not barcoded
#FGC0042_s_1 is normal not barcoded
$(READY_FILES): $(NONBARCODED_SEQUENCES)
	cp $< $@

#align
$(BAMDIR)/%.sam: $(DECODEDDIR)/%.fq
	$(NOVOCMD) -f $< -d $(REFS)/$(GENOME)  > $@

#sam2bam
$(BAMDIR)/%.bam: $(BAMDIR)/%.sam
	samtools view -b -S $< -t $(REFS)/$(TEMPLATE) > $@_tmp
	samtools sort $@_tmp $(BAMDIR)/$*
	rm $@_tmp

#index
$(BAMDIR)/%.bam.bai: $(BAMDIR)/%.bam
	samtools index $<