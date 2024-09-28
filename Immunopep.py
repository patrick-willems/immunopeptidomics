"""Impens immunopeptidomics multi-search engine pipeline."""

import os
import sys
import logging
import argparse
import pickle
import shutil
from pathlib import Path
from core import *      # Import all pipeline functions.

# Run the Impens immunopeptidomics pipeline.
def main():
    """
    Run the immunopeptidomics pipeline.
    For more information see https://github.com/patrick-willems/immunopeptidomics/.
    """

    # Parse the CLI arguments.
    parser = argument_parser()
    args = parser.parse_args()
    indir = args.in_dir

    # Check CLI arguments and necessary input files.
    CLI_check(args)

    # Set-up the folder structure and log file.
    prepare_folders(args)

    # Check the input MS files, determine instrument and load annotation.
    instrument, runs = determine_runs(args)

    # Check the FASTA file and append decoys/contaminants if necessary
    fasta_file, prot_info = prepare_fasta(args)
    
    # Run search engine and mokapot.
    if args.search:
        make_folder(f'{indir}/report/search_res/') # Prepare input folder.
        
        # Run engines: MSFragger is used to generate mzML for Comet, Sage and PEAKS Studio!
        if args.frag: run_MSFragger(args, fasta_file, instrument, runs) # MSFragger-4.1
        if args.comet: run_Comet(args, fasta_file, instrument, runs)    # Comet
        if args.sage: run_Sage(args, fasta_file, instrument, runs)      # Sage v0.14.7

    # Gather scan information.
    scan_folder = f'{indir}/report/scans'
    if os.path.exists(f'{scan_folder}/scans.pkl') and os.path.exists(f'{scan_folder}/index2scan.pkl'):
        logging.info(f'Gathering previously saved scan information.')
        with open(f'{scan_folder}/scans.pkl', 'rb') as scan_file: scans = pickle.load(scan_file)
        with open(f'{scan_folder}/index2scan.pkl', 'rb') as index_file: index2scan = pickle.load(index_file)
    else:
        scans, index2scan = read_mzML(args, runs,)   # Parse and save scan info from mzMLs.
    
    # Run mokapot on PEAKS results file - we need the scans information for this.
    if args.search and args.peaks:
        run_PEAKS(args, index2scan)

    # MS2Rescore.
    if args.rescore:
        make_folder(f'{indir}/report/ms2rescore/') # Prepare input folder.
        
        # Re-score per search engine: slight differences for input conversions.
        if args.frag: rescore_MSFragger(args, instrument, scans)         # MSFragger-4.1
        if args.comet: rescore_Comet(args, instrument, scans)                   # Comet
        if args.sage: rescore_Sage(args, instrument, scans)                     # Sage v0.14.7
        if args.peaks: rescore_PEAKS(args, instrument, scans)                   # PEAKS Studio 12

    # FlashLFQ quantification.
    if args.quant and instrument == 'orbitrap': run_FlashLFQ(args, runs)

    # Quit if no reporting required
    if not args.report:
        logging.info('No report specified, quitting here! Otherwise specify -report!')
        quit()

    # Read now all the mokapot output before and after rescoring.
    peptides, psms, scans = read_results(args, scans, runs, prot_info, instrument)

    # Background check to host species if specified.
    if args.host != '0':
        peptides = run_BLASTP(args, peptides)

    # Make plot folder and remove previous plots if any.
    if os.path.exists(f'{indir}/report/plots/'):
        shutil.rmtree(f'{indir}/report/plots/')    # Remove prior plots!
    make_folder(f'{indir}/report/plots/')   # Many subfolders with different plots will be stored here.
    
    # Polygon plot for timsTOF data.
    if instrument == 'timsTOF':
        polygon_plot(args, scans)
    
    # Peptide length histogram plot.
    histogram_plot(args, peptides)
    
    # Run netMHCpan-4.1 and/or netMHCIIpan-4.3.
    if args.HLA != 'skip':
        peptides, scans = run_netMHCpan(args, peptides, scans)
        netMHCpan_plot(peptides, args)

    # Identification - MS2Rescore plots.
    rescore_plot(args)                      # Ridgeplot target-decoy per search engine.
    rescoring_features(args, peptides, instrument)      # Plot rescoring features performance.
    peptide_barplot(args, peptides)         # Barplot peptide identifications per search engine before and after rescoring.
    numbers_per_run(args, peptides)         # Number of PSMs/peptides per run.
    VennOverlap(peptides, args)             # Peptide overlap per search engine.

    # Run GibbsCluster2.0.
    if args.gibbs != 'skip':
        peptides = run_Gibbs(args, peptides)
        gibbs_plot(args, peptides)

    # Following are only relevant for infection experiments.
    if args.host != '0':
        # Plot spectra for non-host/contaminant proteins.
        if not args.no_spectra_plots:
            plot_spectra(args, scans, instrument)

        # Plot number of epitopes per antigen.
        antigen_plot(args, peptides)

    # Print out scans, PSMs and peptides overview
    scans.to_csv(f'{indir}/report/scans/scans.txt',sep='\t',index=False)
    psms.to_csv(f'{indir}/report/psms.txt',sep='\t')
    peptides.to_csv(f'{indir}/report/peptides.txt',sep='\t',index=False)

    # DONE!
    logging.info('DONE!')

# Parse CLI arguments.
def argument_parser() -> argparse.ArgumentParser :
    """
    Parses command-line arguments for the Impens immunoproteomics pipeline.

    The parser is designed to configure input files, specify search engines, 
    enable rescore, generate reports, and handle quantification settings.
    Main workflow steps are -search/-rescore/-quant/-report.
    For each workflow the search engine(s) of interest have to be specified.
    Returns:
        argparse.ArgumentParser
            Configured argument parser with all options.
    """
    
    parser = argparse.ArgumentParser(
            description="Impens immunoproteomics pipeline.",
            epilog="Contact: Patrick Willems (patwille.willems [at] ugent.be)"
    )

    # Input and workflow determination
    parser.add_argument(
        "-i", "--in_dir", required=True, type=Path, 
        help="Directory containing the .d/.raw/.mzML files."
    )
    parser.add_argument(
        "-search", action="store_true", 
        help="Perform database searches. Specify search engines (e.g., -comet)."
    )
    parser.add_argument(
        "-rescore", action="store_true", 
        help="Run MS2Rescore to refine peptide-spectrum match (PSM) confidence."
    )
    parser.add_argument(
        "-report", action="store_true", 
        help="Generate post-hoc reports including plots and tables."
    )
    parser.add_argument(
        "-quant", action="store_true", 
        help="Enable FlashLFQ quantification for protein or peptide-level data."
    )

    # SEARCH ENGINES
    parser.add_argument(
        "-comet", action="store_true", 
        help="Run Comet and/or parse existing Comet search results."
    )
    parser.add_argument(
        "-frag", action="store_true", 
        help="Run MSFragger and/or parse existing MSFragger search results."
    )
    parser.add_argument(
        "-sage", action="store_true", 
        help="Run Sage and/or parse existing Sage search results."
    )
    parser.add_argument(
        "-peaks", action="store_true", 
        help="Read PEAKS results (place 'db.psms.csv' in /report/search_res/PEAKS/)."
    )
    parser.add_argument(
        "-fa", "--fasta", required=True, type=Path, 
        help="FASTA file required for the database search."
    )
    mods = ['mod', 'NL', 'nomod', 'TMT16', 'TMT10', 'mhcii', 'lowres']
    parser.add_argument(
        "-mod", default="mod", required=False, type=str, choices=mods, 
        help=(
            "Select modification types: 'nomod' (M oxidation), 'mod' (M oxidation, "
            "Cysteinylation, pyroGlu's, Nt acetylation), 'TMT10', 'TMT16', 'mhcii' "
            "(longer peptides), or 'lowres' (low-resolution ion trap data)."
        )
    )

    # REPORTING
    parser.add_argument(
        "-len", default="9", required=False, type=str, 
        help="Specify immunopeptide length (default: 9 for human MHCI). Can be set as an interval (e.g., '8-11')."
    )
    parser.add_argument(
        "-HLA", default="human", required=False, type=str, 
        help="Select HLA type: 'JY', 'HeLa', or 'mouse'."
    )
    parser.add_argument(
        "-gibbs", default="auto", required=False, type=str, choices=["skip", "2", "3", "4", "5", "auto"], 
        help=(
            "Specify the number of clusters for GibbsCluster2 (default: 'auto', "
            "which uses the cluster with the highest KLD). Use 'skip' to disable Gibbs clustering."
        )
    )
    parser.add_argument(
        "-plot_chimera", action="store_true",
        help="Plot annotated spectra for spectra assigned to multiple peptides by different engines."
    )
    parser.add_argument(
        "-contaminant", default="CON__", required=False, type=str, 
        help="Pattern to identify contaminants (default: 'CON__', as used in MaxQuant)."
    )
   
    # Only relevant for infection set-ups:
    parser.add_argument(
        "-no_spectra_plots", action="store_true", 
        help="Do not plot annotated MS2 spectra."
    )
    parser.add_argument(
        "-host", default="0", required=False, type=str,
        help=(
            "Input the species tag (e.g., HUMAN) to perform a background check and plot spectra "
            "of peptides not matching the species (unless -no_spectra_plots is specified)."
        )
    )

    return parser

def CLI_check(args):
    """
    Validates command-line arguments for the Impens immunoproteomics pipeline.

    This function performs several checks to ensure the required inputs are specified and valid.
    
    Parameters:
        args : argparse.Namespace
            Parsed command-line arguments.

    Raises:
        SystemExit: If any of the required conditions are not met.
    """
    
    # Required input folder or FASTA does not exist.
    if not os.path.exists(args.in_dir):
        logging.critical(f'Specified input directory \'{args.in_dir}\' (with MS files) does not exist! Exiting.')
        raise SystemExit
    if not os.path.exists(args.fasta):
        logging.critical(f'Specified FASTA \'{args.fasta}\' does not exist! Exiting.')
        raise SystemExit
    
    # No search engine specified.
    if not (args.sage or args.peaks or args.comet or args.frag):
        logging.critical("At least a single search engine should be specified: -frag, -comet, -sage or -peaks.")
        raise SystemExit

    # No mode specified.
    if not (args.search or args.rescore or args.quant or args.report):
        logging.critical("At least select one process to be executed: -search, -rescore, -quant or -report.")
        raise SystemExit

    # If a search has to be conducted a search engine should be enabled!
    if args.search and not (args.comet or args.peaks or args.sage or args.frag):
        logging.critical("Specify a search engine to be used (e.g. -comet).")
        raise SystemExit

# Make folders
def make_folder(folder: str) -> None:
    """
    Creates a new directory if it does not already exist and sets its permissions.

    Parameters:
        folder : str
            The path of the folder to be created.
    """
    
    if not os.path.exists(folder):
        os.makedirs(folder)
        os.chmod(folder, 0o777)

def prepare_folders(args) -> None:
    """
    Prepares the necessary directories and initializes logging for the immunopeptidomics pipeline.

    Parameters:
        args: argparse.Namespace
            Command-line arguments.
    Returns:
        None.
    """

    #Generate the report folder
    make_folder(f'{args.in_dir}/report/')

    #Generate the new log file
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        handlers=[
            logging.StreamHandler(sys.stdout),  # Output to STDOUT and log file.
            logging.FileHandler(f'{args.in_dir}/report/immunopeptidomics.log', mode='w')
        ]
    )
    logger = logging.getLogger(__name__)

    #Print the full argument list to the log file.
    invoked_cmd = ' '.join(sys.argv)
    arguments = '\n'.join(f'{arg:20}: {getattr(args, arg)}' for arg in vars(args))
    logging.info(f'Command used for immunopeptidomics pipeline:\n'
                f'{invoked_cmd}'
                f'Command line arguments:\n'
                f'{arguments}')

if __name__ == "__main__":
    main()
