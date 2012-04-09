include ../progdir.mk


MIN_LENGTH := 15

#assume this is sanger-graded
novo_loose := $(NOVOALIGN)  -l 17 -h 60 -t 60 -o sam -o FullNW 
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
PREFASTQ_FILES :=          $(addprefix $(PREFASTQDIR)/,$(notdir $(UNCOMPRESSED_FILES:.txt=.fq)))
TRIMMED_FILES :=        $(addprefix $(TRIMMEDDIR)/,$(notdir $(PREFASTQ_FILES:.fq=.$(MIN_LENGTH).fq)))
BARCODED_SEQUENCES:=    $(addsuffix .$(MIN_LENGTH).fq,$(addprefix $(TRIMMEDDIR)/,$(BARCODED_LANES)))
NONBARCODED_SEQUENCES:= $(addsuffix .$(MIN_LENGTH).fq,$(addprefix $(TRIMMEDDIR)/,$(NONBARCODED_LANES)))
DEBARCODED_FILES:=      $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(BC_SAMPLES)))
READY_FILES:=           $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(NONBC_SAMPLES)))
DECODED_FILES :=        $(addsuffix .fq,$(addprefix $(DECODEDDIR)/,$(SAMPLES)))
#decodeddir is the same as fastqdir which is used in shared.mk

decoded: $(DECODED_FILES)
trimmed: $(TRIMMED_FILES)
prefastq: $(PREFASTQ_FILES)
uncompressed: $(UNCOMPRESSED_FILES)

unique:$(UNIQUE_FILES)

solexa:
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0036_s_1_sequence.txt.gz $(SOURCEDIR)/FGC0036_s_1.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0031_s_8_sequence.txt.gz $(SOURCEDIR)/FGC0031_s_8.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/559-Ganguly-Chao-Solexa/basic/Solexa/FGC0036_s_2_sequence.txt.gz $(SOURCEDIR)/FGC0036_s_2.txt.gz
	ln -s -f /nas/is1/leipzig/Ganguly/RB_miRNA_raw_data/589-Ganguly-Chao-Solexa_smRNA/basic/Solexa/FGC0042_s_1_sequence.txt.gz $(SOURCEDIR)/FGC0042_s_1.txt.gz


#Rules
$(SOURCEDIR)/%.txt:$(SOURCEDIR)/%.txt.gz
	gunzip -c $< > $@

#convert to fastq
#convert to sanger phred scores
$(PREFASTQDIR)/%.fq:$(SOURCEDIR)/%.txt
	python $(SCRIPTDIR)/python/FGC2fastq.py < $< > $@_tmp
	python $(SCRIPTDIR)/python/Solexa2Sanger.py $@_tmp $@
	rm $@_tmp

#trim adapters
$(TRIMMEDDIR)/%.$(MIN_LENGTH).fq:$(PREFASTQDIR)/%.fq
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


include ../shared.mk
