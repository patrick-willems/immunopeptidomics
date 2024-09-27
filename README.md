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

## Input parameter reference (also see -h)

# Arguments for input and workflow determination
"-i", "--in_dir", required=True, type=Path, help="Directory with the .d/.raw/.mzML files.")
    parser.add_argument('-search', action='store_true', help="Perform database searches. Specify engines! E.g. -comet).")
    parser.add_argument('-rescore', action='store_true', help="Run MS2Rescore.")
    parser.add_argument('-report', action='store_true', help="Perform post-hoc reports (plots, tables, ..).")
    parser.add_argument('-quant', action='store_true', help="Enable FlashLFQ quantification.")

    # Search engines and database
    parser.add_argument('-comet', action='store_true', help="Run Comet and/or parse existing Comet search results.")
    parser.add_argument('-frag', action='store_true', help="Run MSFragger and/or parse existing MSFragger search results.")
    parser.add_argument('-sage', action='store_true', help="Run Sage and/or parse existing Sage search results.")
    parser.add_argument('-peaks', action='store_true', help="Also read-in PEAKS results (put 'db.psms.csv' in /report/search_res/PEAKS/!")
    parser.add_argument("-fa", "--fasta", required=True, type=Path, help="Fasta file, required for search.'")
    mods = ['mod','NL','nomod','rank2','TMT16','TMT10','semi','mhcii','lowres','restricted','tryptic','tryptic_lowres','glygly']
    parser.add_argument("-mod", default='mod', required=False, type=str, choices=mods, help="Choose from: 'nomod' (only M oxidation), mod (M oxidation, Cysteinylation, pyroGlu's, Nt Acetylation), TMT10 or TMT16.'")
    # Reporting options
    parser.add_argument("-len", default='9', required=False, type=str, help="Immunopeptide length. Default 9 for human MHCI. Can set an interval as '8-11'.")
    parser.add_argument("-HLA", default="human", required=False, type=str, help="Choose either 'JY', 'HeLa', or 'mouse'...")
    parser.add_argument("-gibbs", default='auto', required=False, type=str, choices=['skip','2','3','4','5','auto'], help="Number of maximal clusters to run with GibbsCluster2. Default is 'auto', which will use cluster (from 2 to 5) with highest KLD. Set to 'skip' to turn off Gibbs clustering.'")
    parser.add_argument("-host", default='0', required=False, type=str, help="Input the species tag (e.g. HUMAN) to do a background check (BLASTP-like) and plot all spectra of peptides not matching to this species (unless -no_spectra_plots).'")
    parser.add_argument("-contaminant", default='CON__', required=False, type=str, help="Contaminant match pattern, default is 'CON__' (MaxQuant).'")
    parser.add_argument("-no_spectra_plots", action='store_true', help="Do not plot annotated MS2 spectra.'")
    parser.add_argument("-plot_chimera", action='store_true', help="Plot annotated spectra for spectra assigned to multiple peptides by different engines.'")


## Compatibility
This pipeline is currently only tested and used within a Linux Ubuntu environment.
The python script can be easily transformed in a GUI using for instance streamlit.
For MSFragger and Sage searches take up a lot of memory (we use a 500 GB RAM server), though a splitted database search can be performed.  

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
