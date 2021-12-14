## BEELINE benchmarking of Model-X knockoffs

We mmade mminimal mmodifications to BEELINE in order to benchmark various model-X knockoffs.

The main changes we made were:

- evaluate multiple parameter settings of the same model (https://github.com/Murali-group/Beeline/issues/59)
- add a knockoff-based network inference method
- add datasets with protein concentration and RNA production rate revealed
- benchmark FDR and undirected FDR in addition to AUPR et cetera

To run our experiments, use `test_knockoffs.sh`. For more info about the project, see https://github.com/ekernf01/knockoffs_paper . 
For more info about BEELINE, see the [original](https://github.com/Murali-group/Beeline/).
