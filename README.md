# Maximizing immunopeptidomics-based bacterial epitope discovery by multiple search engines and rescoring
This repository provides a computational workflow to search, rescore, quantify and perform downstream immuno-informatics (HLA binding prediction, Gibbs cluster) and plotting.

## External tools
This processing is currently developed using tools listed belows. Not all tools are required to run part of the pipeline workflow. For instance, only searching and re-scoring can be performed using the search engine of choice

### Search engines
- MSFragger-4.1, https://github.com/Nesvilab/MSFragger/
- Sage, https://github.com/lazear/sage
- Comet, https://uwpr.github.io/Comet/
- PEAKS 12, https://www.bioinfor.com/

### Immunoinformatics (stand-alone versions available)
- GibbsCluster 2.0, https://services.healthtech.dtu.dk/services/GibbsCluster-2.0/
- NetMHCpan-4.1, https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
- NetMHCIIpan-4.3, https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/

### Quantification
- FlashLFQ, https://github.com/smith-chem-wisc/FlashLFQ
  
### Python packages
- ms2rescore
- pandas
- pyteomics
- matplotlib
- seaborn
- Logomaker
- venny4py
