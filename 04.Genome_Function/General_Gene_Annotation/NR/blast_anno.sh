#blastp -query ~/baumannii/03.Genome_Component/Gene/Ab.gmhmmp.pep -db ~/data/db/nr/nr -evalue 1e-5 -outfmt 11 -num_threads 24 > blast.ASN
blast_formatter -archive blast.ASN -outfmt 5  -max_target_seqs 1 > blast.xml
python blastxml_to_tabular.py -c "qseqid,pident,evalue,sseqid,salltitles" blast.xml > blast.tab
cat blast.tab | awk -F "]" '{print $1"]"}' | awk '!arr[$1]++' > Ab.gmhmmp.pep.blast.tab

