FILE=$1
OUTPUT=$2
ROOTNAME=$3

INPUT="${FILE/.gff3}"

awk '$7 == "+" {print $0}' $INPUT.gff3 > $INPUT.pos.gff3
awk '$7 == "-" {print $0}' $INPUT.gff3 > $INPUT.neg.gff3

gffread -o $INPUT.pos.gtf -T --sort-alpha --cluster-only $INPUT.pos.gff3
gffread -o $INPUT.neg.gtf -T --sort-alpha --cluster-only $INPUT.neg.gff3

Rscript ~/scripts/tidyUpRename.R $INPUT.pos.gtf $INPUT.neg.gtf $OUTPUT $ROOTNAME

rm $INPUT.pos.* $INPUT.neg.*
