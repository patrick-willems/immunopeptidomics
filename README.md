# Maximizing immunopeptidomics-based bacterial epitope discovery by multiple search engines and rescoring
This repository provides a computational workflow to search, rescore, quantify and perform downstream immuno-informatics (HLA binding prediction, Gibbs cluster) and plotting.

## Workflows
There are four consecutive steps in the workflow that can be executed:
1. Search (-search)
Currently supports MSFragger, Comet and Sage searches. PSMs and search engine scoring features are used for mokapot processing.
 
2. Rescore (-rescore)
TSV inputs are generated from search results per engine and MS2Rescore/TIMS2Rescore is ran per engine.

3. Quantify (-quant)
FlashLFQ is ran for the identified PSMs. Currently only supports Orbitrap (.RAW) data!

4. Report (-report)
Results for all search engines are integrated at a mokapot 1% peptide-level FDR (after rescoring). A variety of immuno-informatic tools and quality plots are generated. 
This includes NetMHCpan-based MHC binding predictions and GibbsCluster2.0, as well as peptide length histograms, identifications per runs, logo's etc. plotted using matplotlib and seaborn.

## Input Parameter Reference

### Arguments for Input and Workflow Determination:
- `-i`, `--in_dir` **(required)**:  
  Directory with the `.d` or `.raw` files, or preprocessed `.mzML` files.
  
- `-search`:  
  Perform database searches. Specify engines (e.g. `-comet`).
  
- `-rescore`:  
  Run MS2Rescore.
  
- `-report`:  
  Perform post-hoc reports (plots, tables, etc.).
  
- `-quant`:  
  Enable FlashLFQ quantification. Only for Orbitrap data currently!

### Search Engines and Database:
**For everyone workflow (search/rescore/quant/report) the search engine results that you want to consider have to be specified.**

- `-comet`:  
  Run Comet and/or parse existing Comet search results.
  
- `-frag`:  
  Run MSFragger and/or parse existing MSFragger search results.
  
- `-sage`:  
  Run Sage and/or parse existing Sage search results.
  
- `-peaks`:  
  Read instructions [here](https://github.com/patrick-willems/immunopeptidomics/blob/main/README.md#peaks-studio-12).
  
- `-fa`, `--fasta` **(required)**:  
  Fasta file required for search. Decoys need a 'rev_' prefix, they are generated automatically if not present.

### Modification Options:
- `-mod` **(default: 'mod')**:  
  Choose modification type from:
  - `'nomod'`: Only Methionine oxidation.
  - `'mod'`: Methionine oxidation, Cysteinylation, pyroGlus, N-terminal acetylation.
  - `'TMT10'`: Fixed TMT10-plex N-termini and Lys, variable methionine oxidation and cysteinylation.
  - `'TMT16'`: Fixed TMT16-plex N-termini and Lys, variable methionine oxidation and cysteinylation.

### Reporting Options:
- `-len` **(default: '9')**:  
  Immunopeptide length (default: 9 for human MHCI). Set an interval using `'8-11'`.

- `-HLA` **(default: 'human')**:  
  For instance `'JY'`, `'A549'`, `'HeLa'`, or `'mouse'`. Specify your own in dictionary in core.py.

- `-gibbs` **(default: 'auto')**:  
  Number of maximal clusters to run with GibbsCluster2. Options:
  - `'skip'`: Turn off Gibbs clustering.
  - `'2'`, `'3'`, `'4'`, `'5'`, `'auto'`: Auto will use the cluster (from 2 to 5) with the highest KLD.

- `-contaminant` **(default: 'CON__')**:  
  Contaminant match pattern, default is `'CON__'` (MaxQuant).

- `-plot_chimera`:  
  Plot annotated spectra for spectra assigned to multiple peptides by different engines.

**Only relevant to bacterial infection set-ups:**

- `-host` **(default: '0')**:  
  Input the species tag (e.g., `HUMAN`) for background check (BLASTP-like). Plots all spectra of peptides not matching this species (unless `-no_spectra_plots` is used).

- `-no_spectra_plots`:  
  Do not plot annotated MS2 spectra.

## Compatibility
This pipeline is currently only tested and used within a Linux Ubuntu environment.
The python script can be easily transformed in a GUI using for instance streamlit.
For MSFragger and Sage searches take up a lot of memory (we use a 500 GB RAM server), though a splitted database search can be performed. 

### PEAKS Studio 12
It is advised to first run the search with MSFragger-4.1 (`-search -frag`) and use the preprocessed mzML also used for MSFragger/Sage/Comet searches, this will facilitate later PSM result integration. When loading data within PEAKS 12, select the box for preprocessed data (no refinement or chimera options). For rescoring, PEAKS Studio 12 configuration (peaks.conf) has to be adapted to export decoys by updating `"show-debug-panel" : true,` and adding below the line `"export-decoy" : true,`. Importantly, **all PSMs have to be exported for rescoring**, this can be achieved by setting the required -lg10P threshold to 0, thus not applying any scoring/FDR thresholds. The exported PSM report (db.psms.csv) has to be placed in the correct output folders (../report/search_res/PEAKS/). Then, PEAKS rescoring can be done by specifying `-search -rescore -peaks`. In a next step, rescored results can be considered together with other search engines by only running the reporting using `-report -sage -comet -frag -peaks`.

## External tools
This processing is currently developed using tools listed below that have to be downloaded/installed. Not all tools are required to run part of the pipeline workflow. For instance, only searching and re-scoring can be performed using the search engine of choice

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
- ms2rescore, https://github.com/compomics/ms2rescore
- mokapot, https://github.com/wfondrie/mokapot
- pandas, https://github.com/pandas-dev/pandas
- pyteomics, https://github.com/levitsky/pyteomics
- matplotlib, https://github.com/matplotlib/matplotlib
- seaborn, https://github.com/mwaskom/seaborn
- Logomaker, https://github.com/jbkinney/logomaker
- venny4py, https://github.com/timyerg/venny4py
