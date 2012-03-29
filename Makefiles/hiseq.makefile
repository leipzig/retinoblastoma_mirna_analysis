fastq:
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_7/HWI-ST431_52_7_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_7_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_6/HWI-ST431_52_6_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_6_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_5/HWI-ST431_52_5_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_5_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_4/HWI-ST431_52_4_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_4_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_3/HWI-ST431_52_3_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_3_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_2/HWI-ST431_52_2_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_2_1.export.txt.fq
	ln -s /nas/is1/leipzig/Ganguly/gangulyRBhi/flowcell17/lane_1/HWI-ST431_52_1_1.export.txt.fq  $$PWD/fastq/HWI-ST431_52_1_1.export.txt.fq
	ln -s $$PWD/fastq/HWI-ST431_53_8_1.export.txt.fq  $$PWD/fastq/RB525T.fq
	ln -s $$PWD/fastq/HWI-ST431_52_7_1.export.txt.fq  $$PWD/fastq/RB517T.fq
	ln -s $$PWD/fastq/HWI-ST431_52_6_1.export.txt.fq  $$PWD/fastq/RB498T.fq
	ln -s $$PWD/fastq/HWI-ST431_52_5_1.export.txt.fq  $$PWD/fastq/RB498N.fq
	ln -s $$PWD/fastq/HWI-ST431_52_4_1.export.txt.fq  $$PWD/fastq/RB495T.fq
	ln -s $$PWD/fastq/HWI-ST431_52_3_1.export.txt.fq  $$PWD/fastq/RB495N.fq
	ln -s $$PWD/fastq/HWI-ST431_52_2_1.export.txt.fq  $$PWD/fastq/RB494T.fq
	ln -s $$PWD/fastq/HWI-ST431_52_1_1.export.txt.fq  $$PWD/fastq/RB494N.fq

include ../progdir.mk

#assume this is illumina-graded and not adapter-trim
novo_loose := $(NOVOALIGN)  -l 17 -h 60 -t 60 -o sam -o FullNW -a ATCTCGTATGCCGTCTTCTGCTTG  -F ILMFQ
novo_tight := $(NOVOALIGN)  -l 17 -h 0 -t 0 -o sam -o FullNW -a ATCTCGTATGCCGTCTTCTGCTTG  -F ILMFQ
SAMPLES := RB494N RB494T RB495N RB495T RB498N RB498T RB517T RB525T 

include ../shared.mk