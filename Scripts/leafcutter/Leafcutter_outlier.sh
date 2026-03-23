#!/usr/bin/bash
#SBATCH --job-name=leaf
#SBATCH --output=leaf.%j.out
#SBATCH --error=leaf.%j.err

module load leafcutter/0.2.9
module load R

###########
# Depleted
###########

# leafcutter outlier detection
Rscript /scg/apps/software/leafcutter/0.2.9/scripts/leafcutterMD.R --num_threads 8 -o D D_perind_numers.counts.gz

# annotate 
Rscript annotate_counts_with_verdict.R -a annotation_codes/gencode_v48 -o D_pVals_bh_sig_hits.annotated.txt D_pVals_bh_sig_hits.txt.gz

###########
# Standard
###########

# leafcutter outlier detection
Rscript /scg/apps/software/leafcutter/0.2.9/scripts/leafcutterMD.R --num_threads 8 -o WM WM_perind_numers.counts.gz

# annotate 
Rscript annotate_counts_with_verdict.R -a annotation_codes/gencode_v48 -o WM_pVals_bh_sig_hits.annotated.txt  WM_pVals_bh_sig_hits.txt.gz