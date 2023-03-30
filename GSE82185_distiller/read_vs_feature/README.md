# Reproduction steps
## Get interactions between promoter and enhancer regions
1. `make` to get extended Hi-C read regions
1. `./count_interaction.*.sh` for each combination of samples and annotations. Do not run scripts in parallel for the same sample.
1. `./filter_intrachr_1k.sh *.feature_interaction_detail.tsv.gz`
1. `./get_enhancer_ratio*.py`

