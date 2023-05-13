#!/usr/bin/env bash
liftover="/home/hisakatha/softwares/ucsc/liftOver"
chain="/home/hisakatha/softwares/ucsc/mm9ToMm10.over.chain.gz"
domains="GSE63525_mouse_lymphoblasts_Arrowhead_domainlist.txt"
loops="GSE63525_mouse_lymphoblasts_HiCCUPS_looplist.txt"

# Convert domains into bed and then lift it over
cat "$domains" | tail -n+2 | cut -f 1-3 -d$'\t' > "$domains.bed"
$liftover "$domains.bed" $chain "$domains.mm10.bed" "$domains.mm10_unmapped.bed"
cat "$domains.mm10.bed" | sed -e '1 i\chr1\tx1\tx2' > "$domains.mm10.bed.w_header"

# Convert loops into two beds, lift it over, and then merge
cat "$loops" | tail -n+2 | awk -v OFS="\t" '{print $1,$2,$3,"region"NR}' > "$loops.first.bed"
cat "$loops" | tail -n+2 | awk -v OFS="\t" '{print $4,$5,$6,"region"NR}' > "$loops.second.bed"
for ord in first second; do
    $liftover "$loops.$ord.bed" $chain "$loops.$ord.mm10.bed" "$loops.$ord.mm10_unmapped.bed"
done

# Join two beds based on the fourth column (name)
# xsv: A CSV toolkit. https://github.com/BurntSushi/xsv
xsv join 4 "$loops.first.mm10.bed" 4 "$loops.second.mm10.bed" -d "\t" | xsv select 1-3,5-8 | xsv fmt -t "\t" > "$loops.mm10.bedpe"
cat "$loops.mm10.bedpe" | sed -e '1 i\chr1\tx1\tx2\tchr2\ty1\ty2\tname' > "$loops.mm10.bedpe.w_header"

