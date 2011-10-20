#!/bin/sh
bwa aln -t 8 /share/apps/genome/human/bwa/GRCh37/human_g1k_v37.fasta /storage/Ganguly/gangulyRBhi/RB494N.fq > /storage/Ganguly/gangulyRBhi/RB494N.fq.sai
bwa samse /share/apps/genome/human/bwa/GRCh37/human_g1k_v37.fasta /storage/Ganguly/gangulyRBhi/RB494N.fq.sai /storage/Ganguly/gangulyRBhi/RB494N.fq > /storage/Ganguly/gangulyRBhi/RB494N.fq.bwa.sam
samtools view -bt /share/apps/genome/human/bwa/GRCh37/human_g1k_v37.fai -o /storage/Ganguly/gangulyRBhi/RB494N.fq.bwa.bam /storage/Ganguly/gangulyRBhi/RB494N.fq.bwa.sam
samtools sort /storage/Ganguly/gangulyRBhi/RB494N.fq.bwa.bam /storage/Ganguly/gangulyRBhi/RB494N.fq.sorted.bwa
samtools index /storage/Ganguly/gangulyRBhi/file.sorted.bwa.bam

