#/bin/bash

VEP=/opt/ensembl-vep/vep
DATA=/home/projects/archive/cisplatin/VEP

in_rsID=$DATA/kalra_hg38_for_VEP.aa
out_rsID=$DATA/kalra_hg38_for_VEP.aa_results
CA=/home/projects/archive/cisplatin/VEP

perl $VEP \
--input_file $in_rsID \
--format "vcf" \
--output_file $out_rsID \
--dir_cache $CA \
--offline \
--cache
