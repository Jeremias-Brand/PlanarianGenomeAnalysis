FILE=$1
FASTA=$2
INPUT="${FILE/.gtf}"
KEEP=$3

~/bin/TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl $INPUT.gtf $FASTA > $INPUT.fasta
~/bin/TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl $INPUT.gtf > $INPUT.ali.gff3

~/bin/TransDecoder-v5.5.0/TransDecoder.LongOrfs -S -m 120 -t $INPUT.fasta

~/bin/TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl \
     $INPUT.fasta.transdecoder_dir/longest_orfs.gff3 \
     $INPUT.ali.gff3 \
     $INPUT.fasta > $INPUT.rawCDS.gff3

mv $INPUT.fasta.transdecoder_dir/longest_orfs.pep $INPUT.CDS.pep
mv $INPUT.fasta.transdecoder_dir/longest_orfs.cds $INPUT.CDS.cds

grep ">" $INPUT.fasta > $INPUT.names

# This step below is only required with Transdecoder versions up to v5.5.0
# (the bug has been fixed in v5.7.0)
Rscript ~/scripts/fixCDS_raw.R $INPUT.rawCDS.gff3 $INPUT.CDS.cds

rm -rf $INPUT.fasta $INPUT.ali.gff3 $INPUT.fasta.transdecoder_dir* pipeliner.*
if [[ "$KEEP" != "keeptmp" ]]; then rm $INPUT.rawCDS.gff3; fi

