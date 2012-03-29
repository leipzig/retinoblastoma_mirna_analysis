#Directories
TOP := $(shell pwd)
FASTQDIR := $(TOP)/fastq
BAMDIR := $(TOP)/bam
SCRIPTDIR := $(TOP)/src
RCSDIR:= $(TOP)/rcs

#these are Solexa-specific
SOURCEDIR := $(TOP)/sources
DECODEDDIR := $(TOP)/decoded
FASTADIR := $(TOP)/fasta
TRIMMEDDIR := $(TOP)/trimmed

#Programs
SAM_POSS:= /share/apps/bin/samtools /usr/bin/samtools
SAMTOOLS:= $(wildcard $(SAM_POSS))
NOVO_POSS:= /share/apps/bin/novoalign /usr/bin/novoalign
NOVOALIGN:= $(wildcard $(NOVO_POSS))
BEDTOOLS:= $(TOP)/../exe/bedtools/bedtools
