for f in *fa; do cat $f | fastx_trimmer -l 26 | fastx_collapser | fasta_formatter -t | perl -ne 'm/\d+\-(\d+)\t(\S+)/;print $2."\t".$1."\n";' > $f.col & done
