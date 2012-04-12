#Directories
TOP := $(shell pwd)
FASTQDIR := $(TOP)/fastq
BAMDIR := $(TOP)/bam
SCRIPTDIR := $(TOP)/src
RCSDIR:= $(TOP)/rcs

#these are Solexa-specific
SOURCEDIR := $(TOP)/sources
DECODEDDIR := $(TOP)/fastq
FASTADIR := $(TOP)/fasta
TRIMMEDDIR := $(TOP)/trimmed
PREFASTQDIR := $(TOP)/prefastq

#Programs
SAM_POSS:= /nas/is1/leipzig/exe/samtools
# /share/apps/bin/samtools /usr/bin/samtools
SAMTOOLS:= $(wildcard $(SAM_POSS))
NOVO_POSS:= /share/apps/bin/novoalign /usr/bin/novoalign
NOVOALIGN:= $(wildcard $(NOVO_POSS))
BEDTOOLS:= $(TOP)/../exe/bedtools/bedtools
