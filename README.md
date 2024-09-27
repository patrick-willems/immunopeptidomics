# Maximizing immunopeptidomics-based bacterial epitope discovery by multiple search engines and rescoring
This repository provides a computational workflow to search, rescore, quantify and perform downstream immuno-informatics (HLA binding prediction, Gibbs cluster) and plotting.

## External tools
This processing is currently developed using tools listed belows. Not all tools are required to run part of the pipeline workflow. For instance, only searching and re-scoring can be performed using the search engine of choice

### Search engines
- MSFragger-4.1
- Sage
- Comet
- PEAKS 12

### Immunoinformatics
- GibbsCluster 2.0
- NetMHCpan-4.1
- NetMHCIIpan-4.3

### Python packages
- ms2rescore
- pandas
- pyteomics
- matplotlib
- seaborn
- venny4py
