#!/usr/bin/env bash
set -e

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate star279a
star=STAR
ulimit -n 10000
#tmpdir=$(mktemp -d /tmp/star.XXXXXX)
# make a temporary directory using process ID (accessible as $$)
tmpdir=/tmp/star.$$

if [[ $(uname) == "Darwin" ]]; then
    # MacOS
    zcat=gzcat
else
    # Linux
    zcat=zcat
fi
echo "[INFO] .gz extractor: $zcat"

genome_dir="$HOME/Desktop/reference/Mus_musculus/Ensembl_direct/GRCm38/star_index_150/genome"
fastq_dir="/Users/shigenseigyo/Desktop/analysis/personal/funaya/Dr_Funaya_210210_remove_rib_Fastq/"
fastq_prefix="$(basename $(pwd))"
suffix1=".fastq.1.gz"
suffix2=".fastq.2.gz"
$star --runThreadN 8 --genomeDir $genome_dir --outSAMtype BAM SortedByCoordinate --outTmpDir=$tmpdir \
 --readFilesIn "$fastq_dir$fastq_prefix$suffix1" "$fastq_dir$fastq_prefix$suffix2" --readFilesCommand $zcat \
 --outFilterMultimapNmax 10 --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.2 --outSAMstrandField intronMotif --sjdbScore 2
