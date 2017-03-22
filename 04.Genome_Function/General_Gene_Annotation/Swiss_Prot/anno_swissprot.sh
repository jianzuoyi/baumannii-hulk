blastp -query ~/baumannii/03.Genome_Component/Gene/Ab.gmhmmp.pep -db ~/data/db/swissprot/swissprot -evalue 1e-5 -outfmt 11 -num_threads 12 > Ab2swissprot.ASN
blast_formatter -archive Ab2swissprot.ASN -outfmt 5  -max_target_seqs 1 > blast.xml
python ~/tools/blastxml_to_tabular.py -c "qseqid,pident,evalue,sseqid,salltitles" blast.xml > blast.tab

