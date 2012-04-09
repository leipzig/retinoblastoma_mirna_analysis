SUBDIRS = doSolexa doHiseq doPGM

#.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)    

$(SUBDIRS):
	$(MAKE) -C $@

doSolexa: doPGM


STRINGENCIES := loose tight
SHARING_MODES := maximum gene-est-only
COUNTS:= 50 100 200

OUTPUTDIR:=results/reports
DESCRIPT:=src/R/diffExp.Rnw
R_CMD:=/usr/bin/R --vanilla --slave
CONFIG:=src/R/dirConfig.R
diffexp:=$(foreach strng,$(STRINGENCIES),$(foreach sm,$(SHARING_MODES),$(foreach cnt,$(COUNTS),$(OUTPUTDIR)/de$(strng)$(sm)$(cnt).pdf)))

RCSDIR:=rcs

diffexp: $(diffexp)

NOVO_POSS:= /share/apps/bin/novoindex /usr/bin/novoindex
NOVOINDEX:= $(wildcard $(NOVO_POSS))

define report
 results/reports/de$(1)$(2)$(3).tex: src/R/diffExp.Rnw
	echo "configFile<-'$(CONFIG)';alignmentStringency='$(1)';sharingMode='$(2)';minCount=$(3);Sweave('src/R/diffExp.Rnw',output='results/reports/de$(1)$(2)$(3).tex')" | $(R_CMD)
endef

%.pdf:%.tex
	pdflatex -interaction nonstopmode -output-directory $(<D) $<
	pdflatex -interaction nonstopmode -output-directory $(<D) $<
	rm -fv $*.tex $*.aux $*.log $*.toc $*.out
	rm -fv $(*F)-*.pdf

$(foreach strng,$(STRINGENCIES),$(foreach sm,$(SHARING_MODES),$(foreach cnt,$(COUNTS),$(eval $(call report,$(strng),$(sm),$(cnt))))))

#11 RCS files:
#6 paired
#517 525
#WERI
#Y-79
#42
RB494N:=doHiseq/fastq/RB494N.fq
RB494T:=doHiseq/fastq/RB494T.fq
RB495N:=doHiseq/fastq/RB495N.fq
RB495T:=doHiseq/fastq/RB495T.fq
RB498N:=doHiseq/fastq/RB498N.fq 
RB498T:=doHiseq/fastq/RB498T.fq
RB517T:=doHiseq/fastq/RB517T.fq
WERI  :=doSolexa/decoded/FGC0031_s_8.WERI.fq doSolexa/decoded/FGC0036_s_2.WERI.fq doPGM/fastq/WERI01.fq doPGM/fastq/WERI02.fq doPGM/fastq/WERI03.fq
Y79   :=doSolexa/decoded/FGC0036_s_1.Y79.fq doSolexa/decoded/FGC0036_s_2.Y79.fq
NRML  :=doSolexa/decoded/FGC0042_s_1.normal.fq
RB525T :=doHiseq/fastq/RB525T.fq doPGM/fastq/RB525T01.fq  doPGM/fastq/RB525T02.fq  doPGM/fastq/RB525T03.fq 

CONCATGROUPS:=RB494N RB494T RB495N RB495T RB498N RB498T RB517T WERI Y79 NRML RB525T
CONCATDIR:=concats

CONCATTARGETS:=$(foreach concatGroup,$(CONCATGROUPS),$(CONCATDIR)/$(concatGroup).fq)
RCSTARGETS:=$(foreach concatGroup,$(CONCATGROUPS),$(RCSDIR)/$(concatGroup).rcs)

concats:$(CONCATTARGETS)
rcs:$(RCSTARGETS)

define makeConcats
 $(CONCATDIR)/$(1).fq:$($(1))
	cat $$+ > $$@
endef

$(foreach concatGroup,$(CONCATGROUPS),$(eval $(call makeConcats,$(concatGroup))))

$(RCSDIR)/%.rcs: $(CONCATDIR)/%.fq
	cat $< | fastx_trimmer -l 26 -Q 33 | fastx_collapser -Q 33 | fasta_formatter -t | perl -ne 'm/\d+\-(\d+)\t(\S+)/;print $$2."\t".$$1."\n";' > $@

refs/hairpin.fa.gz:
	wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz refs/hairpin.fa.gz

refs/hairpin.fa:refs/hairpin.fa.gz
	gunzip -c refs/hairpin.fa.gz > refs/newhairpin.fa
	perl -ne 'unless(/^>/){s/U/T/g;}print' < refs/newhairpin.fa > refs/hairpin.dna.fa
	perl -ne 'BEGIN{$$/=">"}s/>//g;if(/hsa/){print ">".$$_}' refs/hairpin.dna.fa > refs/hsa_hairpin.dna.fa
	mv refs/hsa_hairpin.dna.fa refs/hairpin.fa
	$(NOVOINDEX) refs/hairpin.ndx refs/hairpin.fa

tRNAs.txt:tRNAs.unsorted.txt
	cat tRNAs.unsorted.txt | sort -S 10G --batch-size 64 -u -T /nas/is1/leipzig/doHiseq/ > tRNAs.txt

tRNAs.unsorted.txt:
	for f in do*/bam/novo/tight/tRNAs/all/*bam; do samtools view -F 4 $$f | cut -f 10 >> tRNAs.unsorted.txt; done

notmirna.txt:
	sort -m -u do*/bam/novo/tight/hg19.ambig/all/*.notmirna.sorted.txt > notmirna.txt

refs/hsa.chr.gff:
	wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff refs/hsa.gff
	perl -ne 'unless(/^#/){print "chr".$_}' < refs/hsa.gff > refs/hsa.chr.gff
