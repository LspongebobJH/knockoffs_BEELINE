
FROM r-base:4.0.2

LABEL Maintainer="Eric Kernfeld <ekernfe1@jh.edu>"

USER root

WORKDIR /

# slow but usually cached 
RUN R -e "install.packages(c('knockoff', 'data.table'))"
RUN R -e "install.packages(c('magrittr'))"
RUN mkdir data/
RUN apt-get update && apt-get install time
RUN R -e "install.packages(c('FNN'))"
RUN R -e "install.packages(c('optparse'))"
RUN R -e "install.packages(c('ggplot2'))"
RUN R -e "install.packages(c('quantreg'))"
RUN R -e "install.packages(c('tidyverse'))"
RUN R -e "install.packages(c('mclust'))"
RUN R -e "install.packages(c('GeneNet'))"

# frequently altered
COPY rlookc_0.1.0.tar.gz /
RUN R -e "install.packages('rlookc_0.1.0.tar.gz', type = 'source')"
COPY runGeneNet.R /

