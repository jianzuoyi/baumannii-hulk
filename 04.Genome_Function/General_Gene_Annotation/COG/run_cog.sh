wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt /home/zuoyi/db/COG/
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog /home/zuoyi/db/COG/
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG/myva /home/zuoyi/db/COG/
makeblastdb -in ~/db/COG/myva -dbtype prot -titile cog -parse_seqids -out ~/db/COG/cog
blastp -query CP019114-5.fa -db ~/db/COG/cog -evalue 1e-5 -out CP019114-5.blast.tab -outfmt 6 -max_target_seqs 1 -num_threads 12
