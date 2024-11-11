#!/bin/bash

# Run experiments testing model-X knockoffs for inference of the BEELINE toy networks.
# conda activate BEELINE
 #./initialize.sh

python BLRunner.py --config config-files/knockoff_experiments/dyn-BFC.yaml

# for toy in BFC BF CY LI LL TF
# do
#     python BLRunner.py --config config-files/knockoff_experiments/dyn-${toy}.yaml
#     # python BLEvaluator.py --config config-files/knockoff_experiments/dyn-${toy}.yaml --auc
#     # python BLEvaluator.py --config config-files/knockoff_experiments/dyn-${toy}.yaml --ufdr
#     # python BLEvaluator.py --config config-files/knockoff_experiments/dyn-${toy}.yaml --fdr
# done

# Rscript plotKnockoffEvaluation.R
