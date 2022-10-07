## BEELINE benchmarking of Model-X knockoffs

We made minimal modifications to BEELINE in order to benchmark various model-X knockoffs.

Our changes:

- evaluate multiple parameter settings of the same model (issue [here](https://github.com/Murali-group/Beeline/issues/59))
- add the GeneNet network inference method
- add a knockoff-based network inference method
- add datasets with protein concentration and RNA production rate revealed (small modification of BoolODE)
- benchmark FDR and undirected FDR in addition to AUPR et cetera
- Add a little R script for the final plots

To run our experiments, use `test_knockoffs.sh`. For more info about the project, see out [project main repo](https://github.com/ekernf01/knockoffs_paper). For more info about BEELINE, see the [original](https://github.com/Murali-group/Beeline/).
To run our experiments:

- Install Docker and Conda
- Set up the BEELINE conda environment via `setupAnacondaVENV.sh` 
- Run `test_knockoffs.sh`. 

For more info about this project, see [the central repo](https://github.com/ekernf01/knockoffs_paper) . For more info about BEELINE, see the [original](https://github.com/Murali-group/Beeline/).
