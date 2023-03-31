#!/usr/bin/env bash
ln ../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gff3.gz
ln ../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gff3.gz.tbi

ln ../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene_and_nc.gff3.gz
ln ../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene_and_nc.gff3.gz.tbi

ln ../../mouse_ref/Ensembl_direct/GRCm38/Sequence/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
ln ../../mouse_ref/Ensembl_direct/GRCm38/Sequence/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai

for d in ../NoRibs_*; do
    ln $d/Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed read_template.$(basename $d).bed
    ln $d/Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed.tdf read_template.$(basename $d).bed.tdf
    ln $d/Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed.cov.bigWig read_template.$(basename $d).bed.cov.bigWig
    ln $d/Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed.cpm.bigWig read_template.$(basename $d).bed.cpm.bigWig
done
