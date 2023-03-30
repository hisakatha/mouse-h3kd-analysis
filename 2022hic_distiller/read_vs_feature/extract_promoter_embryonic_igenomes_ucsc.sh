#!/usr/bin/env bash
set -u
# rg: ripgrep, a kind of grep
export PATH="/home/g57008/.cargo/bin:/work/ga17/share/tools/bedtools2/bin:$PATH"
export LC_ALL=C
# TPM list for embryonic genes
target_list="Embryonicgene_5stage_average_tpm.tsv"
# Gene annotation in Illumina iGenomes mouse UCSC mm10 reference
annot="/work/ga17/share/references/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf.gz"
# Chromosome lengths in UCSC format
fai="/work/ga17/share/references/ucsc_mm10/mm10.chrom.sizes"
# Length of output promoter regions
length=500

zcat $annot | awk '$3=="exon"' |
    rg -F -f <(cat $target_list | tail -n+2 | cut -f 1 | sed -E 's/(.+)/gene_name "\1";/') |
    awk -v FS="\t" -v OFS="\t" '{gene_name = gensub(/.*gene_name "([^";]+)";.*/, "\\1", 1, $9); print $1,$4-1,$5,gene_name,$6,$7}' |
    sort -k4,4 -k1,1 | bedtools groupby -g 1,4 -c 2,3 -o min,max -full |
    awk -v OFS="\t" '{print $1,$7,$8,$4,$5,$6}' | bedtools flank -i - -g $fai -l $length -r 0 -s |
    sort -k1,1 -k2,2n > promoter.embryonic_igenomes_ucsc.bed

cat $target_list | tail -n+2 | cut -f 1 | rg -v -F -f <(zcat $annot | awk '$3=="exon"' | sed -E 's/.*gene_name "([^;]+)".*/\1/') > genes_not_in_igenomes_ucsc.txt

