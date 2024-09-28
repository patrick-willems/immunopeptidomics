''' Core script containing all the code for searches, rescoring, and quantification '''

import os
import re
import sys
import logging
import glob
import shutil
import pandas as pd
import pymzml
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import logomaker as lm
import subprocess
from typing import Tuple
from pyteomics import mgf
from pyteomics import proforma
from collections import Counter
from venny4py.venny4py import venny4py
from matplotlib.gridspec import GridSpec
from io import StringIO
logging.getLogger('matplotlib').setLevel(logging.WARNING)   # Otherwise font type spam.
logger = logging.getLogger(__name__)

# HLA dictionary for netMHCpan-4.1 and netMHCIIpan-4.3 alleles per cell type.
hla_dict = {
    "netMHCpan_human": "HLA-A01:01,HLA-B15:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01",
    "netMHCpan_mouse": "H-2-Db,H-2-Dd,H-2-Dq,H-2-Kb,H-2-Kd,H-2-Kk,H-2-Kq,H-2-Ld,H-2-Lq",
    "netMHCpan_C57Bl6": "H-2-Kb,H-2-Db",
    "netMHCpan_JY": "HLA-A02:01,HLA-B07:02,HLA-C07:02",
    "netMHCpan_HeLa": "HLA-A03:19,HLA-A68:02,HLA-B15:03,HLA-C12:03",
    "netMHCpan_U937": "HLA-A03:01,HLA-A31:01,HLA-B18:01,HLA-B51:01,HLA-C01:02,HLA-C07:01",
    "netMHCpan_H1650": "HLA-A02:01,HLA-B15:18,HLA-C07:04",
    "netMHCpan_HCT116": "HLA-A01:01,HLA-A02:01,HLA-B45:01,HLA-B18:01,HLA-C07:01,HLA-C05:01",
    "netMHCpan_THP1": "HLA-A02:01,HLA-B15:11,HLA-C03:03",
    "netMHCpan_HEK293": "HLA-A02:01,HLA-A03:01,HLA-B07:02,HLA-C07:02",
    "netMHCpan_A549": "HLA-A25:01,HLA-A30:01,HLA-B18:01,HLA-B44:03,HLA-C12:03,HLA-C16:01",
    "netMHCpan_H9": "HLA-A02:01,HLA-A03:01,HLA-B35:03,HLA-B44:02,HLA-C04:01,HLA-C07:04",
    "netMHCIIpan_JY": "HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10402,HLA-DQA10103-DQB10302,HLA-DQA10103-DQB10603,HLA-DQA10301-DQB10302,HLA-DQA10301-DQB10603,DRB1_0404,DRB1_1301,DRB4_0104,DRB5_0105",
    "netMHCIIpan_U937": "DRB1_1601,DRB1_0101,DRB1_1501,DRB5_0101,HLA-DQA10101-DQB10501,HLA-DQA10102-DQB10501,HLA-DQA10101-DQB10602,HLA-DPA10202-DPB10402,HLA-DPA10202-DPB10201,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10401",
    "netMHCIIpan_THP1": "DRB1_0101,DRB1_1501,DRB5_0101,HLA-DQA10101-DQB10501,HLA-DQA10102-DQB10501,HLA-DQA10102-DQB10602,HLA-DQA10101-DQB10602,HLA-DPA10202-DPB10402,HLA-DPA10202-DPB10201,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10402",
    "netMHCIIpan_human": "DRB1_0101,DRB3_0101,DRB4_0101,DRB5_0101",
    "netMHCIIpan_HeLa": "DRB1_0101,DRB3_0101,DRB4_0101,DRB5_0101",
    "netMHCIIpan_mouse": "H-2-IAb",
    "netMHCIIpan_C57Bl6": "H-2-IAb",
    "netMHCIIpan_H9": "DRB1_1501,DRB1_1601",
    "netMHCIIpan_HCT116": "DRB1_0305,DRB1_1102,HLA-DQA10502-DQB10202",
    "netMHCIIpan_HEK293": "DRB1_1501,HLA-DQA10102-DQB10602",
    "netMHCIIpan_A549": "DRB1_0701,DRB1_1104,HLA-DQA10201-DQB10202,HLA-DQA10201-DQB10301,HLA-DQA10505-DQB10202,HLA-DQA10505-DQB10301"
}

def detect_instrument_from_mzml(mzml_file: str) -> str:
    """
    Detect the instrument type from an mzML file by reading its contents.

    Parameters:
    ----------
    mzml_file : str
        Path of mzML file to read in.

    Returns:
    -------
        A string indicating the instrument type, either 'orbitrap' or 'timsTOF'.
    """
    with open(mzml_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= 20000:
                return 'orbitrap'
            if "ion mobility" in line:
                return 'timsTOF'
    return 'orbitrap'

# Determine the input file types.
def determine_runs(args):
    """
    Determine MS input files, their filetype and provided annotation.

    This function checks the input files (.RAW, .d, .mzML) and identifies the instrument type 
    (Orbitrap or timsTOF) based on the provided MS files. If provided by the user, an annotation
    file will be parsed to add sample names and condition. If the condition is 'uninfected',
    this is used later on for bacterial peptide filtering.

    Parameters:
    ----------
    args : argparse.Namespace
        The arguments parsed by the argparse module.

    Returns:
    -------
    instrument : str
        A string indicating the instrument type, either 'Orbitrap' or 'timsTOF'.
    
    runs : pd.DataFrame
        A pandas DataFrame containing information about the MS input files/runs.
    """

    logger.info(f'Determining input file types and possible annotation.')
    indir = args.in_dir
    ms_data = {}   # Store MS file info here.

    # Check directory for input files.
    raw_files = list(glob.iglob(os.path.join(indir,'*raw')))         # Orbitrap files
    mzml_files = list(glob.iglob(os.path.join(indir,'*mzML')))       # Preprocessed mzML files
    mzml_files = [mzml for mzml in mzml_files if 'calibrated.mzML' not in mzml]                 # Ignore intermediary files
    d_folders = [ f.path for f in os.scandir(indir) if f.is_dir() and f.path.endswith('.d') ]    # timsTOF file directories
    
    # Check for valid input files.
    if not any([d_folders, mzml_files, raw_files]):
        logging.critical("No orbitrap (.raw), timsTOF (.d) or preprocessed data (.mzML) found in input directory.")
        raise FileNotFoundError("No MS files found in the input directory.")

    # Check for mixed instrument types.
    if raw_files and d_folders:
        logging.critical("Orbitrap and timsTOF input data found in input directory! Should be single instrument data as input.")
        raise ValueError("Mixed instrument data detected.")

    # Determine instrument type.
    if d_folders:
        ms_data.update({'MS_files': sorted(d_folders)})
        instrument = 'timsTOF'
    elif raw_files:
        ms_data.update({'MS_files': sorted(raw_files)})
        instrument = 'orbitrap'
    else: instrument = None

    # Preprocessed (.mzML) data found, double check if conform to raw MS data.
    if mzml_files:
        #There is no raw MS data: determine file type from mzML.
        if not instrument:
            instrument = detect_instrument_from_mzml(mzml_files[0])
        
        # If also raw MS data found, compare to raw data (if any).
        else:
            ms_basenames = [os.path.splitext(f)[0] for f in (raw_files or d_folders)]
            mzml_basenames = [os.path.splitext(f)[0] for f in mzml_files]
            for filename in mzml_basenames:
                if filename not in ms_basenames:
                    logging.critical(f'mzML file {filename} does not matching a raw MS file in the input directory.')
                    raise ValueError(f"mzML file mismatch: {filename}.")

        #Add mzML files as raw data.
        if 'MS_files' not in ms_data:
            ms_data.update({'MS_files': sorted(mzml_files)})

    # Initiate pandas dataframe with runs and add the basename.
    runs = pd.DataFrame(ms_data)
    runs['Filenames'] = runs.iloc[:,0].apply(lambda x: os.path.splitext(os.path.basename(x))[0])    # Filenames.
    runs['mzML_files'] = runs['MS_files'].apply(lambda x: os.path.splitext(x)[0] + '.mzML')         # Preprocessed mzML.
    logging.info(f"Found {len(runs)} {instrument} file(s).")

    # Append annotation.
    annotation_f = os.path.join(indir, 'annotation.txt')
    if os.path.exists(annotation_f):
        logging.info('Annotation file found, will append to runs.')
        annotation = pd.read_csv(annotation_f, sep='\t')
        annotation.columns = ['MS_file','name','condition']
        annotation['Filenames'] = annotation['MS_file'].apply(lambda x: os.path.splitext(x)[0]) #Remove extension to match by filename

        # If our MS files are specified within the annotation file its good to go.
        if runs['Filenames'].isin(annotation['Filenames']).all():
            runs = pd.merge(runs,annotation, on='Filenames', how='inner')
            if len(annotation) > len(runs):
                logging.warning('There were more MS files specified in the annotation file than there were MS files in the input directory.')
        else:
            logging.warning('Some of the MS files were not specified in annotation file, no annotation will be added.')

    return instrument, runs

def parse_protein_header(header, args):
    """
    Parse FASTA protein header line to extract relevant fields such as protein ID, description, gene name, 
    accession number, and respective species. Recommended to use UniProtKB!

    Parameters:
        header : str
            FASTA protein header.
        args : argparse.Namespace
            The arguments parsed by the argparse module.

    Returns:
        dict
            A dictionary containing the description, gene, UniProtKB accession and species.
    """

    # Get Protein ID, description, gene and accession from UniProtKB formatted headers.
    prot_id = header.split()[0][1:]
    descr = re.search(r"^\S+ (.+)", header).group(1) if re.search(r"^\S+ (.+)", header) else prot_id
    gene = re.search(r"GN=(\S+)", header).group(1) if re.search(r"GN=(\S+)", header) else prot_id
    accession = prot_id.split('|')[1] if prot_id.count('|') == 2 else prot_id

    # Assign species/contaminant protein. Contaminant tag default 'CON__' (MaxQuant).
    if 'OS=' in header:
        species = prot_id.split('_')[-1]
    elif args.contaminant in header:
        species = 'contaminant'
    else:
        species = 'unknown'
    
    return {'description': descr, 'gene': gene, 'acc': accession, 'species': species}

def make_folder(folder):
    """
    Creates a directory if it does not exist and give full read/write/execute permissions.

    Parameters:
        folder (str): The path to the folder to be created.

    Returns:
        None
    """

    if not os.path.exists(folder):
        os.mkdir(folder)
        os.chmod(folder, 0o777)     # Give user permissions

def check_file(file):
    """
    Checks whether the specified file exists.

    Parameters:
        file (str): The path to the file to check.

    Returns:
        None
    """

    if not os.path.exists(file):
        logging.critical(f"Path {file} not found!")
        raise SystemExit

def prepare_fasta(args):
    """
    Check target/decoy numbers and parse protein information from specified FASTA file.
    If no decoys are found a concatenated target-decoy databse will be made.

    Parameters:
        args : argparse.Namespace
            The arguments parsed by the argparse module.

    Returns:
        tuple: (fasta, prot_info)
            fasta: The path to the FASTA file (modified if decoys are generated).
            prot_info: A dictionary containing protein information (description, gene, accession, species).
    """

    # First check existence specified FASTA file.
    fasta = str(args.fasta)
    check_file(fasta)

    # Initiate dictionaries and counts.
    target, decoy = [0]*2   # Count how many target and decoy sequences present.
    prot_info = {}          # Here all protein information is stored.
    target_seq = {}         # If no decoys, we need to revert these and append.

    # Read FASTA file and gather protein information
    with open(fasta, 'r') as in_f:
        prot_id = None
        for line in in_f:
            line = line.rstrip()
            if line.startswith('>rev_'):
                decoy += 1
            elif line.startswith('>'):
                target += 1
                # Save protein info (description, species, ..) and initiate sequence storage.
                prot_id = line.split()[0][1:]                   # Extract protein ID without '>'.
                prot_id_peaks = line.split()[0][4:]             # For PEAKS, which chops off 'sp|' or 'tr|'
                header_info = parse_protein_header(line, args)
                prot_info[prot_id] = header_info
                prot_info[prot_id_peaks] = header_info
                target_seq[prot_id] = ''
            elif prot_id:
                target_seq[prot_id] += line

    logging.info(f'{os.path.basename(args.fasta)} contains {decoy} decoy and {target} target proteins.')

    # Generate and concatenate decoys if none are found
    if decoy == 0:
        fasta = re.sub(r'\.\w+$', '_tda.fasta', fasta)
        logging.warning('No decoys found, generating reverse target-decoy database.')

        with open(fasta, 'w') as out_fasta:
            for prot_id, seq in target_seq.items():
                rev_seq = seq[::-1] # Reverse
                out_fasta.write(f'>{prot_id} {prot_info[prot_id]["description"]} GN={prot_info[prot_id]["gene"]}\n'
                                f'{seq}\n'   # Target
                                f'>rev_{prot_id} decoy_{prot_id}\n'
                                f'{rev_seq}\n') # Decoy
        logging.info(f'Generated {os.path.basename(fasta)}!')

    return fasta, prot_info

def run_mokapot(pin_folder, args) -> None:
    """
    Concatenates PIN files from a specified folder, performs engine-specific adjustments, 
    and runs mokapot on search engine features.

    Parameters:
        pin_folder : str
            Path to the folder containing engine-specific PIN files.
        args : argparse.Namespace
            Command-line arguments.

    Returns:
        None
    """

    mokapot_folder = os.path.join(args.in_dir, 'report', 'mokapot')
    make_folder(mokapot_folder)
    engine = os.path.split(pin_folder)[-1]
    logging.info(f'Running mokapot for {engine} results - see log at report/mokapot/{engine}_search_log.txt.')
    
    # Concatenate all PIN files and keep the header of the first file.
    pin_file = f'{mokapot_folder}/{engine}_search_combined.pin'
    concat_cmd = f"awk 'NR == 1 || FNR > 1'  {pin_folder}/*.pin > {pin_file}"
    subprocess.run(concat_cmd, shell=True, check=True)

    # Search engine-specific adjustments PIN file (Comet and MSFragger).
    if engine == 'Comet':           # Merge trailing tab-delimited alternative proteins.
        with open(pin_file, 'r') as pin, open(f'{mokapot_folder}/tmp.pin', 'w') as out:
            for i, line in enumerate(pin):
                splitted = line.split('\t')
                if i == 0:
                    columns = len(splitted)-1    # Number of columns in the header.
                out.write('\t'.join(splitted[:columns]) + '\t' + ';'.join(splitted[columns:]))
        shutil.move(f'{mokapot_folder}/tmp.pin', pin_file)
    elif engine == 'MSFragger':       # Remove trailing semicolon.
        os.system(f"awk '{{sub(/;$/, \"\"); print}}' {pin_file} > {mokapot_folder}/tmp.pin")
        shutil.move(f'{mokapot_folder}/tmp.pin', pin_file)

    # Run mokapot. We use a test FDR of 5% as this is sometimes required for samples with low peptide identifications.
    mokapot_log = f'{mokapot_folder}/{engine}_search_log.txt'
    try:
        mokapot_cmd = f"mokapot --keep_decoys --test_fdr 0.05 -d {mokapot_folder} -r {engine}_search {pin_file} > {mokapot_log} 2>&1"
        subprocess.run(mokapot_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Mokapot failed for {engine}: {e}")
        raise

def run_MSFragger(args, fasta_file, instrument, runs) -> None:
    """
    Executes MSFragger search with the specified parameters, processes search results.
    Pre-made parameter files are present in tools/MSFragger-4.1/params/.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        fasta_file : str
            Path to the FASTA file to use for the MSFragger search.
        instrument : str
            Instrument being orbitrap/timsTOF, important for parameter file selection.
        runs : pd.DataFrame
            DataFrame containing information about the runs/files to be searched.

    Returns:
        None
    """

    logging.info('Starting MSFragger-4.1 search!')

    # Prepare folder and check executable!
    indir = args.in_dir
    msfragger_folder = f'{indir}/report/search_res/MSFragger'
    make_folder(msfragger_folder)
    exe = 'external_tools/MSFragger-4.1/MSFragger-4.1.jar'
    check_file(exe)

    # Check parameter file and copy to output folder with correct FASTA specified.
    param_file = f'external_tools/MSFragger-4.1/params/{instrument}_{args.mod}.params'
    check_file(param_file)
    os.system(f"sed 's#your_fasta#{fasta_file}#g' {param_file} > {msfragger_folder}/used.params")

    # Construct command
    cmd = f'java -Dfile.encoding=UTF-8 -Xmx490g -jar {exe} {msfragger_folder}/used.params '
    
    # Always use raw MS files if available - otherwise mzMLs.
    if 'MS_files' in runs: cmd += ' '.join(runs['MS_files'].tolist())
    else: cmd += ' '.join(runs['mzML_files'].tolist())
    
    # Execute
    log_file = os.path.join(msfragger_folder, 'search_log.txt')
    logging.info(f'Executing MSFragger - see full log in {log_file}')
    try:
        subprocess.run(f'{cmd} > {log_file} 2>&1', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"MSFragger failed: {e}")
        raise

    # Clean-up pepindices and pepXML generated.
    map(os.remove, glob.glob(f'{fasta_file}*pepindex'))
    map(os.remove, glob.glob(f'{indir}/*pepXML'))

    # Move percolator input files to search results folder.
    for pin_file in glob.glob(f'{indir}/*.pin'):
        shutil.move(pin_file, f'{msfragger_folder}/{os.path.basename(pin_file)}')

    # Keep calibrated mzMLs or otherwise uncalibrated ones (sometimes not enough data for calibration).
    for mzml in glob.glob(f'{indir}/*_calibrated.mzML'):
        final_mzml = mzml.replace('_calibrated','') 
        shutil.move(mzml, final_mzml)
    for mzml in glob.glob(f'{indir}/*_uncalibrated.mzML'):
        final_mzml = mzml.replace('_uncalibrated','')   
        if os.path.exists(final_mzml):  # We have a calibrated version - remove.
            os.remove(mzml)             
        else:
            shutil.move(mzml, final_mzml)
    
    # Run mokapot on combined PIN file.
    run_mokapot(msfragger_folder, args)

def run_Comet(args, fasta_file, instrument, runs) -> None:
    """
    Executes a Comet search using specified parameters and mzML files.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        fasta_file : str
            Path to the FASTA file to use for the MSFragger search.
        instrument : str
            Instrument being orbitrap/timsTOF, important for parameter file selection.
        runs : pd.DataFrame
            DataFrame containing information about the runs/files to be searched.

    Returns:
        None.
    """

    logging.info('Starting Comet v2024.01.0 search!')
    indir = args.in_dir

    # Check if mzMLs were generated by MSFragger.
    mzML_files = runs['mzML_files']
    for mzML in mzML_files:
        if not os.path.exists(mzML):
            logging.critical(f"{mzML} not found - first run MSFragger.")
            raise SystemExit

    # Prepare folder and check executable!
    comet_folder = f'{indir}/report/search_res/Comet'
    make_folder(comet_folder)
    exe = 'external_tools/Comet/comet.linux.exe'
    check_file(exe)

    # Check parameter file and copy to output folder with correct FASTA specified.
    param_file = f'external_tools/Comet/params/{instrument}_{args.mod}.params'
    check_file(param_file)
    shutil.copy(param_file,f'{comet_folder}/used.params')

    # Construct and execute command.
    log_file = f'{comet_folder}/search_log.txt'
    logging.info(f'Executing Comet - see full log in {log_file}.')
    cmd = f'{exe} -D{fasta_file} -P{param_file} ' + ' '.join(mzML_files)
    try:
        subprocess.run(f'{cmd} > {log_file} 2>&1', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Comet failed: {e}")
        raise

    # Move output PIN files to report folder and run mokapot
    for pin_file in glob.glob(f'{indir}/*.pin'):
        shutil.move(pin_file, f'{comet_folder}/{os.path.basename(pin_file)}')
    
    # Run mokapot.
    run_mokapot(comet_folder, args)

def run_Sage(args, fasta_file, instrument, runs) -> None:
    """
    Executes a Sage search using specified parameters and mzML files.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        fasta_file : str
            Path to the FASTA file to use for the MSFragger search.
        instrument : str
            Instrument being orbitrap/timsTOF, important for parameter file selection.
        runs : pd.DataFrame
            DataFrame containing information about the runs/files to be searched.

    Returns:
        None.
    """

    logging.info('Starting Sage v0.14.7 search!')
    indir = args.in_dir

    # Check if mzMLs were generated by MSFragger.
    mzML_files = runs['mzML_files']
    for mzML in mzML_files:
        if not os.path.exists(mzML):
            logging.critical(f"{mzML} not found - first run MSFragger.")
            raise SystemExit

    # Prepare folder and check software!
    sage_folder = f'{indir}/report/search_res/Sage'
    make_folder(sage_folder)
    exe = 'external_tools/Sage/sage'
    check_file(exe)

    # Retrieve calibrated MS2 tolerance and TopN peaks settings from MSFragger search
    log_path = f'{indir}/report/search_res/MSFragger/search_log.txt'
    if os.path.exists(log_path):
        with open(log_path, 'r') as in_file:
            for line in in_file:
                if 'New fragment_mass_tolerance' in line: ms2_tol = line.rstrip().split()[3]
                elif 'New use_topN_peaks' in line: max_peaks = line.rstrip().split()[3]
    
    # Check parameter file and copy to output folder with correct FASTA specified.
    param_file = f'external_tools/Sage/params/{instrument}_{args.mod}.json'
    check_file(param_file)
    
    # Adapt search configuration with calibrate search settings if found.
    if ms2_tol is not None:
        logging.info(f'Using MSFragger calibrated search settings: fragment ion tolerance of {ms2_tol} ppm and top {max_peaks} peaks in spectrum.')
        with open(param_file, 'r') as in_f, open(f'{sage_folder}/used_params.json', 'w') as out_f:
            for line in in_f:
                if 'fragment_tol' in line: out_f.write(f'  "fragment_tol": {{"ppm": [-{ms2_tol},{ms2_tol}]}},\n')
                elif 'max_peaks' in line: out_f.write(f'  "max_peaks": {max_peaks},\n')
                else: out_f.write(line)
    else:
        shutil.copy(param_file,f'{sage_folder}/used_params.json')

    # Construct and execute command.
    cpus = int(os.cpu_count() / 2)
    cmd = f'{exe} --batch-size {cpus} --write-pin -o {sage_folder} -f {fasta_file} {sage_folder}/used_params.json ' + ' '.join(mzML_files)
    os.environ['RUST_MIN_STACK'] = '4194304'    # Set to stack size 4 MB to avoid overflow error.
    log_file = f'{sage_folder}/search_log.txt'
    logging.info(f'Executing Sage - see full log in {sage_folder}.')
    try:
        subprocess.run(f'{cmd} > {log_file} 2>&1', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Sage failed: {e}")
        raise

    # Sage combines a single PIN output file for all runs, so just run mokapot here.
    run_mokapot(sage_folder, args)

# Convert db.psms.csv to mokapot input file.
def run_PEAKS(args, index2scan) -> None:
    """
    Generate a percolator input file for PEAKS from exported PSM report (db.psms.csv).
    This db.psms.csv has to placed manually in the correct folder!

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        index2scan : dict
            Spectrum index to scan title, parsed before by mzML reading.

    Returns:
        None.
    """
    
    # PEAKS Studio PSM report should be exported (db.psms.csv) - this only tested for PEAKS 12 with decoys exported!
    peaks_folder = f'{args.in_dir}/report/search_res/PEAKS'
    csv_file = f'{peaks_folder}/db.psms.csv' 
    check_file(csv_file)

    # Make additional columns for correct PIN format.
    df = pd.read_csv(csv_file,sep=',')
    df['Run'] = df['Source File'].apply(lambda x: x.replace('.mzML',''))    # PEAKS should be ran on preprocessed mzML!
    df['scan_id'] = df['Run'] + '_' + df['Scan'].astype(str)
            
    # PEAKS works with spectrum indices, we convert here as we saved this during mzML reading with pyteomics.
    df['ScanNr'] = df['scan_id'].apply(lambda x: index2scan.get(x))
    missing_scans = df['ScanNr'].isna().sum()
    converted_scans = df['ScanNr'].notna().sum()
    if (converted_scans / len(df)) < 0.98:   # Likely something wrong if less than 98% of the scans are converted!  
        raise ValueError("PEAKS scan indices were not matched to a stored scan title! Did you search mzML in PEAKS Studio? Try re-run.")
    elif missing_scans > 0:
        logging.warning(f'{missing_scans} / {len(df)} scan results from PEAKS not assigned to saved ScanNr parsed from mzMLs.')
    df = df[df['ScanNr'].notna()]

    # Append additional columns.
    df['SpecId'] = df['Run'] + '_' + df['ScanNr'] 
    df['Label'] = df['decoy'].map({True: -1, False: 1})
    df['Proteins'] = df['Accession'].apply(lambda x: x.replace('#DECOY#','rev_'))   # Change decoy prefix to be consistent.

    # Ion mobility scoring feature will only be present for timsTOF data.
    selected_columns = ['SpecId', 'Label', 'ScanNr', '-10LgP', 'Tag Length', 'Delta RT', 'MS2 Correlation', 'ppm', 'Peptide', 'Proteins']
    if 'Delta 1/k0' in df.columns:
        selected_columns.append('Delta 1/k0')
    mokapot_input = df[selected_columns]

    # Save input and run mokapot.
    mokapot_input.to_csv(f'{peaks_folder}/db.psms.converted.pin',sep='\t',index=False)
    run_mokapot(f'{args.in_dir}/report/search_res/PEAKS', args)

def parse_precursor(spectrum) -> Tuple[float, str]:
    """
    Parses the precursor information from a given spectrum object.

    This function attempts to extract the m/z value and charge state of the selected precursor.
    If the normal method of extraction fails due to missing precursor information,
    it falls back to parsing the string representation of the spectrum using regular expressions.

    Parameters:
        spectrum
            Pyteomics spectrum object containing spectral data.

    Returns:
        tuple: A tuple containing:
            prec_mz : float
                The m/z value of the precursor.
            charge : str
                The charge state of the precursor.
    """
    try:    # Normal syntax.
        precursor = spectrum.selected_precursors[0]
        prec_mz = precursor['mz']
        charge = precursor['charge']
    except: # Sometimes precursor not found - regex fix.
        spectrum_string = spectrum.to_string().decode('utf-8')
        prec_mz_match = re.search(r'name="selected ion m/z" value="([\d.]+)"', spectrum_string)
        charge_match = re.search(r'name="charge state" value="(\d+)"', spectrum_string)
        if prec_mz_match and charge_match:
            prec_mz = float(prec_mz_match.group(1))
            charge = charge_match.group(1)
        else:
            logging.warning("Problem with mzML parsing of spectra in run {run} scan {scan}.")
            raise SystemExit

    return prec_mz, charge

def write_MGF(args, mzml, run) -> None:
    """
    Converts an mzML file to MGF format and writes it to disk.

    This function reads a specified mzML file, extracts the precursor information
    and spectral data for each spectrum, and writes it to a MGF file in the scans directory.
    These are handy for MS2Rescore, spectrum plotting and data backup.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        mzml : str
            The path to the mzML file to be converted.
        run: str
            Basename of mzML file.

    Returns:
        None.
    """
    
    with open(f'{args.in_dir}/report/scans/{run}.mgf','w') as mgf_out:
        mzml_run = pymzml.run.Reader(mzml)
        for spectrum in mzml_run:
            # Get precursor info.
            prec_mz, charge = parse_precursor(spectrum)

            #Write the MGF
            mgf_out.write(
                f"BEGIN IONS\n"
                f"TITLE={run}_scan={str(spectrum.id_dict['scan'])}_z={charge}\n"
                f"RTINSECONDS={(spectrum.scan_time_in_minutes()*60):.6f}\n"
                f"PEPMASS={prec_mz:.4f}\n"
                f"CHARGE={charge}+\n"
            )
            for mz, intensity in zip(spectrum.mz,spectrum.i):
                mgf_out.write(f"{mz:.6f}\t{intensity}\n")       # Six digits to reduce file size
            mgf_out.write('END IONS\n')

def read_mzML(args, runs) -> dict:
    """
    Reads scan information from mzML files and save MGF files.

    This function parses each mzML file generated by MSFragger. It extracts spectrum index, titles,
    MS2-level scan data, precursor information, and ion mobility.
    It saves the scan and index2scan dictionaries as pickle files to bypass this lengthy step if the
    script is repeated. The index2scan is used to convert PEAKS results only.

    Parameters:
        runs: dict
            A dictionary containing 'mzML_files', which lists the paths to mzML files.
        args : argparse.Namespace
            Command-line arguments.

    Returns:
        dict: A tuple containing two dictionaries:
            scans
                Dictionary of spectra (run_scanNumber) and associated scan information.
            index2scan
                Dictionary of spectrum index to scan number (for PEAKS results later).
    """

    logging.info('Reading scan info from mzML files and generating MGFs.')
    scans = {}
    index2scan = {}

    # Making MGF files for simplified MS2Rescore spectrum parsing.
    scan_folder = f'{args.in_dir}/report/scans/'
    make_folder(scan_folder)

    # Read through mzML files - anticipating MSFragger generated mzML file format.
    for mzml in runs['mzML_files']:
        run = os.path.splitext(os.path.basename(mzml))[0]
        logging.info(f'Now doing {run}..')

        # Write MGF for simplified parsing MS2Rescore and plots later.
        if not os.path.exists(f'{scan_folder}/{run}.mgf'):
            write_MGF(args, mzml, run)

        # Read mzML for the scans dictionary.
        mzml_run = pymzml.run.Reader(mzml)
        for spectrum in mzml_run:
            if spectrum['ms level'] != 2: continue          # Only MS2-level scans.
            scan = str(spectrum.ID)
            specid = f"{run}_{scan}"
            ion_mobility = spectrum.get('MS:1002815',0)     # Default to zero if missing.
            prec_mz, charge = parse_precursor(spectrum)
            
            # PEAKS scan dict: (index+1) --> Scan
            index2scan.update({f'{run}_{(spectrum.index + 1)}': scan})
            
            # Update scans dict.
            scans.update({specid: {
                'run': run,
                'scan': scan,
                'spectrum_index': spectrum.index, 
                'rt': str(spectrum.scan_time_in_minutes()),
                'mz': str(prec_mz),
                'z': str(charge),
                'im': str(ion_mobility),
                'rescore': f'{run}_scan={scan}_z={charge}_seq='}    # Seq append later for unique PSM identifier.
            })

    # Save this so this can be parsed immediately upon next run of the script if necessary.
    with open(f'{args.in_dir}/report/scans/scans.pkl', 'wb') as scans_file: pickle.dump(scans, scans_file)
    with open(f'{args.in_dir}/report/scans/index2scan.pkl', 'wb') as index_file: pickle.dump(index2scan, index_file)
    
    return scans, index2scan

def get_peptidoform(peptide) -> str:
    """
    Converts modified peptide sequence of search engines to PSI-MOD peptidoform notation.
    These peptidoform notations are compatible with MS2Rescore/TIMS2Rescore.
    A variety of modification tags were encountered dependent per search engine, although
    the mod_mapping dictionary has to be extended if additional PTMs are considered in the search.
    It also handles N-terminal modifications by reformatting them appropriately.

    Parameters:
        peptide : str
            Modified peptide sequence outputted by search engine

    Returns:
        peptide : str
            Peptide adapted to PSI-MOD peptidoform notation.
    """

    # Modification specifications vary according search engine.
    mod_mapping = {
        '[+57.0215]': '[Carbamidomethyl]',
        '+57.021465': 'Carbamidomethyl',
        '[+15.9949]': '[Oxidation]',
        '[+42.010567]': '[Acetyl]',
        '[+119.0041]': '[Cysteinylation]',
        'n[304.2072]': '[TMTpro]-',
        'n[+304.20715]': '[TMTpro]-',
        '[304.2072]': '[TMTpro]',
        '[+304.20715]': '[TMTpro]',
        'n[229.1629]': '[TMT6plex]-',
        '[229.1629]': '[TMT6plex]',
        '[+229.16293]': '[TMT6plex]',
        '(+229.16)': '[TMT6plex]',
        'n[+229.1629]': '[TMT6plex]-',
        '[+229.1629]': '[TMT6plex]',
        '[15.9949]': '[Oxidation]',
        '[57.0215]': '[Carbamidomethyl]',
        '[119.0041]': '[Cysteinylation]',
        '[-17.0265]': '[Gln->pyro-Glu]',
        '[-18.0106]': '[Glu->pyro-Glu]',
        'n[42.0106]': '[Acetyl]-',
        '[-17.026548]': '[Gln->pyro-Glu]',
        '[-18.010565]': '[Glu->pyro-Glu]',
        '(+15.99)': '[Oxidation]',
        '(+42.01)': '[Acetyl]',
        '(-17.03)': '[Gln->pyro-Glu]',
        '(-18.01)': '[Glu->pyro-Glu]',
        '(+119.00)': '[Cysteinylation]',
    }
    
    # Replace the modification names using dict.
    for mod, psimod in mod_mapping.items():
        peptide = peptide.replace(mod, psimod)
    
    # Replace N-terminal modifications as A[TMT6plex]XX to [TMT6plex]-AXX.
    for nt_mod in ['TMT6plex', 'TMTpro', 'Acetyl']:
        pattern = rf'^([A-Z])\[({nt_mod})\]'
        peptide = re.sub(pattern, r'[\2]-\1', peptide)

    return peptide

def prepare_tsv(psms, scans):
    """
    Prepares a DataFrame containing TSV input to be used for MS2Rescore/TIMS2Rescore.
    This function constructs a DataFrame containing various basic PSM metrics (RT, IM, mz, etc.).

    Parameters:
        psms : pd.DataFrame
            DataFrame containing peptide spectrum match data.
        scans : dict
            Dictionary containing scan information indexed by spectrum ID.

    Returns:
        tsv_df : pd.DataFrame
            A DataFrame with TSV input for rescoring.
    """

    # Initiate tsv_df dataframe.
    columns = ['peptidoform','peptide','spectrum_id','run','is_decoy','score','qvalue',
            'pep','precursor_mz','retention_time','protein_list','rank']
    tsv_df = pd.DataFrame(columns=columns)
    
    # Get charge from saved scans dict.
    psms['z'] = psms['specid'].apply(lambda x: scans[x]['z'])

    # Populate the dataframe.
    tsv_df['peptidoform'] = psms['Peptide'].apply(lambda x: get_peptidoform(x)) + '/' + psms['z']   # e.g. AAM[Oxidation]ADSAA/2
    tsv_df['peptide'] = psms['Peptide'].apply(lambda x: re.sub('[^A-Z]+', '', x))                   # e.g. AAMADSAAA
    tsv_df['spectrum_id'] = psms['specid'].apply(lambda x: scans[x]['rescore']) + tsv_df['peptide'] #Filename_scan=xx_z=2_seq=AAMADSAA
    tsv_df['run'] = psms['run']
    tsv_df['is_decoy'] = ~psms['Label']
    tsv_df['score'] = psms['mokapot score']
    tsv_df['qvalue'] = psms['mokapot q-value']
    tsv_df['pep'] = psms['mokapot PEP']
    tsv_df['ion_mobility'] = psms['specid'].apply(lambda x: scans[x]['im'])             # Get ion mobility.
    tsv_df['precursor_mz'] = psms['specid'].apply(lambda x: scans[x]['mz'])             # Get m/z.
    tsv_df['retention_time'] = psms['specid'].apply(lambda x: scans[x]['rt'])           # Get RT.
    tsv_df['protein_list'] = psms['Proteins'].apply(lambda x: x.split(';'))
    tsv_df['rank'] = 1

    return tsv_df

# Read mokapot PSMs from target and decoys
def read_mokapot_search(args, engine: str) -> pd.DataFrame:
    """
    Read and combine mokapot target and decoy mokapot PSM outputs.

    Parameters:
        in_dir : str
            The input directory containing the mokapot report files.
        engine : str
            The name of the search engine used (e.g. 'Comet').

    Returns:
        pd.DataFrame
            A combined DataFrame containing both target and decoy PSMs.
    """

    # File paths.
    mokapot_dir = f'{args.in_dir}/report/mokapot'
    target_file = f'{mokapot_dir}/{engine}_search.mokapot.psms.txt'
    decoy_file = f'{mokapot_dir}/{engine}_search.mokapot.decoy.psms.txt'
    
    # Read target PSMs
    check_file(target_file)
    psms_target = pd.read_csv(target_file,sep='\t')
    
    # Read decoy PSMs
    check_file(decoy_file)
    psms_decoy = pd.read_csv(decoy_file,sep='\t')

    #Concatenate and return.
    psms = pd.concat([psms_target, psms_decoy], ignore_index=True)
    
    return psms

def run_MS2Rescore(tsv_df: pd.DataFrame, features: pd.DataFrame, instrument:str, args, engine: str) -> None:
    """
    Run MS2Rescore, the following steps are followed:
    1. Extend mokapot PSM info with search engine scoring features.
    2. Prepare input and configuration file.
    3. Run MS2Rescore (or TIMS2Rescore for timsTOF data).
    4. Log output and move output files.

    Parameters:
        tsv_df : pd.DataFrame
            DataFrame containing the mokapot results.
        features : pd.DataFrame
            DataFrame containing the search engine scoring features.
        args : argparse.Namespace
            Command-line arguments.
        instrument : str
            The MS instrument, either orbitrap or timsTOF.
        engine : str
            The name of the search engine used.

    Returns:
        None.
    """
    
    ms2rescore_folder = f'{args.in_dir}/report/ms2rescore/'

    # Re-shuffle indices and join the mokapot results with PIN features.
    features.set_index('spectrum_id', inplace=True)
    features = features.add_prefix('rescoring:')
    tsv_df['spectrum_id_index'] = tsv_df['spectrum_id']
    tsv_df.set_index('spectrum_id_index', inplace=True)
    tsv_joint = tsv_df.join(features, how='inner')
    
    # Specify the rank of the search engine if within scoring features (otherwise defaults to 1).
    if 'rescoring:rank' in tsv_joint:
        tsv_join = tsv_joint.pop('rescoring:rank')

    # Print out MS2Rescore input file.
    ms2rescore_input = f'{ms2rescore_folder}/{engine}_rescore_input.tsv'
    tsv_joint.to_csv(ms2rescore_input,sep='\t',index=False)

    # Configuration file.
    config_file = f'external_tools/ms2rescore/{instrument}_{args.mod}.json'
    check_file(config_file)

    # Construct command and run.
    log_file = f'{ms2rescore_folder}/{engine}_pipeline_log.txt'
    cmd = f'ms2rescore -t tsv -s {args.in_dir}/report/scans/ -p {ms2rescore_input} -c {config_file} -o {ms2rescore_folder}/{engine}'
    logging.info(f'Running MS2Rescore for {engine} - see log {log_file}.\nExecuted command: {cmd}')
    try:
        subprocess.run(f'{cmd} > {log_file} 2>&1', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"MS2Rescore failed: {e}")
        raise

    # Print to log
    with open(log_file, 'r') as log_in:
        for line in log_in:
            if 'core // Identified' in line:
                gain = line.rstrip().split('// ')[-1]
                logging.info(f'{gain} PSMs at 1% FDR after rescoring.')

    # Move the mokapot output files: add ms2rescore to avoid confusion with files before rescoring.
    for mokapot_file in glob.glob(f'{ms2rescore_folder}/{engine}.mokapot*txt'):
        destination = f'{args.in_dir}/report/mokapot/' + os.path.basename(mokapot_file.replace(engine,engine + '_ms2rescore'))
        shutil.move(mokapot_file, destination)

def rescore_MSFragger(args, instrument: str, scans: dict) -> None:
    """
    Prepares and executes MS2Rescore for MSFragger search engine results.
    This function processes mokapot PSM output from MSFragger, prepares input files for MS2Rescore,
    appends rescoring features from the corresponding PIN file, and then runs MS2Rescore.
    
    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        instrument : str
            The MS instrument, either orbitrap or timsTOF.
        scans : dict
            Dictionary containing scan information.
    
    Returns:
        None.
    """

    logging.info('Preparing MS2Rescore input for MSFragger results.')
    indir = args.in_dir

    # Check and read mokapot target PSM output.
    psms = read_mokapot_search(args,'MSFragger')
    psms['run'] = psms['SpecId'].apply(lambda x: '.'.join(x.split('.')[:-3]))   # e.g. H083972_uPAC_10_CMB-623_FRIMP_ImPep_LFQ_HCT116-1.12141.12141.3_1
    psms['specid'] = psms['run'] + '_' + psms['ScanNr'].astype(str)             # To retrieve from scans dictionary.
    psms['Peptide'] = psms['Peptide'].apply(lambda x: x[2:-3])                  # e.g. R.SHYEEGPGKNLPFSVENKWS3.L
    
    # Prepare initial columns for TSV input.
    tsv_df = prepare_tsv(psms, scans)

    # Now append search engine features from the percolator input file.
    pin_file = f'{indir}/report/mokapot/MSFragger_search_combined.pin'
    check_file(pin_file)
    features = pd.read_csv(pin_file,sep='\t',dtype=str)
    features['run'] = features['SpecId'].apply(lambda x: '.'.join(x.split('.')[:-3]))   # e.g. H083972_uPAC_10_CMB-623_FRIMP_ImPep_LFQ_HCT116-1.12141.12141.3_1
    features['Peptide'] = features['Peptide'].apply(lambda x: re.sub('[^A-Z]+', '', x[2:-3]))                  # e.g. R.SHYEEGPGKNLPFSVENKWS3.L
    features['specid'] = features['run'] + '_' + features['ScanNr']                     # To retrieve from scans dictionary.
    features['spectrum_id'] = features['specid'].apply(lambda x: scans[x]['rescore']) + features['Peptide']
    features = features.drop(['run','specid','ScanNr','SpecId','Label','Peptide','Proteins'], axis=1)   # Unnecessary columns

    # Run MS2Rescore
    run_MS2Rescore(tsv_df,features,instrument,args,'MSFragger')

def rescore_Comet(args, instrument: str, scans: dict) -> None:
    """
    Prepares and executes MS2Rescore for Comet search engine results.
    This function processes mokapot PSM output from Comet, prepares input files for MS2Rescore,
    appends rescoring features from the corresponding PIN file, and then runs MS2Rescore.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        instrument : str
            The MS instrument, either orbitrap or timsTOF.
        scans : dict
            Dictionary containing scan information.

    Returns:
        None.
    """

    logging.info('Preparing MS2Rescore input for Comet results.')
    indir = args.in_dir

    # Prepare basic TSV input from mokapot PSM reports.
    psms = read_mokapot_search(args,'Comet')
    psms['run'] = psms['SpecId'].apply(lambda x: '_'.join(os.path.basename(x).split('_')[:-3])) # e.g. /folder/H083972_uPAC_10_CMB-623_FRIMP_ImPep_LFQ_HCT116-1_744_2_1
    psms['specid'] = psms['run'] + '_' + psms['ScanNr'].astype(str)                         # To retrieve from scans dictionary.
    psms['Peptide'] = psms['Peptide'].apply(lambda x: x[2:-2])                              # e.g. A.IQRTPKIQVYSRHPAENG.K
    tsv_df = prepare_tsv(psms, scans)

    # Comet does not report fixed modifications, add TMT modifications if necessary.
    if args.mod == 'TMT16':
        tsv_df['peptidoform'] = tsv_df['peptidoform'].apply( lambda x: '[TMTpro]-' + x.replace('K','K[TMTpro]') )
    elif args.mod == 'TMT10':
        tsv_df['peptidoform'] = tsv_df['peptidoform'].apply( lambda x: '[TMT6plex]-' + x.replace('K','K[TMT6plex]') )

    # Now append search engine features from the percolator input file.
    pin_file = f'{indir}/report/mokapot/Comet_search_combined.pin'
    check_file(pin_file)
    features = pd.read_csv(pin_file,sep='\t',dtype=str)
    features['run'] = features['SpecId'].apply(lambda x: '_'.join(os.path.basename(x).split('_')[:-3])) # Extract run.
    features['Peptide'] = features['Peptide'].apply(lambda x: re.sub('[^A-Z]+', '', x[2:-2]))           # Get plain peptide sequence.
    features['specid'] = features['run'] + '_' + features['ScanNr']                                     # To retrieve from scans dictionary.
    features['spectrum_id'] = features['specid'].apply(lambda x: scans[x]['rescore']) + features['Peptide']
    features = features.drop(['run','specid','ScanNr','SpecId','Label','Peptide','Proteins'], axis=1)   # Unnecessary columns

    # Run MS2Rescore
    run_MS2Rescore(tsv_df,features,instrument,args,'Comet')

def rescore_Sage(args, instrument: str, scans: dict) -> None:
    """
    Prepares and executes MS2Rescore for Sage search engine results.
    This function processes mokapot PSM output from Sage, prepares input files for MS2Rescore,
    appends rescoring features from the corresponding PIN file, and then runs MS2Rescore.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        instrument : str
            The MS instrument, either orbitrap or timsTOF.
        scans : dict
            Dictionary containing scan information.

    Returns:
        None.
    """

    logging.info('Preparing MS2Rescore input for Sage results.')
    indir = args.in_dir

    # Prepare basic TSV input from mokapot PSM reports.
    psms = read_mokapot_search(args,'Sage')
    psms['run'] = psms['FileName'].apply(lambda x: x.replace('.mzML', ''))
    psms['specid'] = psms['run'] + '_' + psms['ScanNr'].astype(str)       # To retrieve from scans dictionary.
    tsv_df = prepare_tsv(psms, scans)

    # Now append search engine features from the percolator input file.
    pin_file = f'{indir}/report/mokapot/Sage_search_combined.pin'
    check_file(pin_file)
    features = pd.read_csv(pin_file,sep='\t',dtype=str)
    features['run'] = features['FileName'].apply(lambda x: x.replace('.mzML', ''))
    features['Peptide'] = features['Peptide'].apply(lambda x: re.sub('[^A-Z]+', '', x))
    features['specid'] = features['run'] + '_' + features['ScanNr']       # To retrieve from scans dictionary.
    features['spectrum_id'] = features['specid'].apply(lambda x: scans[x]['rescore']) + features['Peptide']
    features = features.drop(['run','specid','ScanNr','SpecId','Label','Peptide','Proteins','FileName'], axis=1)    # Unnecessary columns

    # Run MS2Rescore
    run_MS2Rescore(tsv_df,features,instrument,args,'Sage')

def rescore_PEAKS(args, instrument: str, scans: dict) -> None:
    """
    Prepares and executes MS2Rescore for PEAKS search engine results.
    This function processes mokapot PSM output from PEAKS, prepares input files for MS2Rescore,
    appends rescoring features from the corresponding PIN file, and then runs MS2Rescore.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        instrument : str
            The MS instrument, either orbitrap or timsTOF.
        scans : dict
            Dictionary containing scan information.

    Returns:
        None.
    """

    logging.info('Preparing MS2Rescore input for PEAKS results.')
    indir = args.in_dir

    # Prepare basic TSV input from mokapot PSM reports.
    psms = read_mokapot_search(args,'PEAKS')
    psms['run'] = psms['SpecId'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    psms['specid'] = psms['SpecId']       # To retrieve from scans dictionary.
    tsv_df = prepare_tsv(psms, scans)

    # Now append search engine features from the percolator input file.
    pin_file = f'{indir}/report/mokapot/PEAKS_search_combined.pin'
    check_file(pin_file)
    features = pd.read_csv(pin_file,sep='\t',dtype=str)
    features['Peptide'] = features['Peptide'].apply(lambda x: re.sub('[^A-Z]+', '', x))
    features['spectrum_id'] = features['SpecId'].apply(lambda x: scans[x]['rescore']) + features['Peptide']
    features = features.drop(['ScanNr','SpecId','Label','Peptide','Proteins','Delta RT'], axis=1) #Delta RT has missing values for double mods on single AA.

    # Run MS2Rescore
    run_MS2Rescore(tsv_df,features,instrument,args,'PEAKS')

def run_FlashLFQ(args,  runs: pd.DataFrame) -> None:
    """
    Prepares and executes FlashLFQ quantification based on Mokapot peptide results.

    This function reads Mokapot peptide output files, processes them to create an input 
    DataFrame for FlashLFQ, and then executes FlashLFQ using the specified command line 
    options. It handles multiple runs, deduplicates peptide entries, and saves the 
    generated input file in a specified directory.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        runs : pd.DataFrame
            A DataFrame containing run information, specifically the MS file paths.

    Returns:
        None
    """

    logging.info('Reading mokapot peptide results for FlashLFQ input.')
    indir = args.in_dir
    quant_dir = f'{indir}/report/quant'
    make_folder(quant_dir)

    # Initiate FlashLFQ input dataframe.
    columns = ['File Name','Scan Retention Time','Precursor Charge','Base Sequence','Full Sequence','Peptide Monoisotopic Mass','Protein Accession'] 
    flashlfq = pd.DataFrame(columns=columns)

    # Add required identified peptide info.
    peptide_files = glob.glob(f'{indir}/report/mokapot/*ms2rescore*mokapot.peptides.txt')
    for mokapot_peptide in peptide_files:
        pepfile = pd.read_csv(mokapot_peptide,sep='\t')
        pepfile['Full Sequence'] = pepfile['peptide']
        pepfile['Base Sequence'] = pepfile['peptide'].apply( lambda x: re.sub(r'\[.*?\]', '', x).replace('-','') )
        pepfile['File Name'] = pepfile['run']
        pepfile['Peptide Monoisotopic Mass'] = pepfile['calcmass']
        pepfile['Protein Accession'] = pepfile['protein_list'].apply( lambda x: x.strip('[]').replace("'","") )
        pepfile['Precursor Charge'] = pepfile['charge']
        pepfile['Scan Retention Time']  = pepfile['retention_time']

        # Reduce the pepfile dataframe to only the specified columns and append.
        pepfile = pepfile[columns]
        flashlfq = pd.concat([flashlfq, pepfile], ignore_index=True)
    
    # Remove duplications introduced by the multiple search engines.
    flashlfq = flashlfq.groupby(['Full Sequence', 'Precursor Charge', 'File Name']).agg({
        'Scan Retention Time': 'first',
        'Base Sequence': 'first', 
        'Peptide Monoisotopic Mass': 'first',
        'Protein Accession': 'first'
    }).reset_index()

    # Save input file.
    flashlfq.to_csv(f'{quant_dir}/flashlfq_input.tsv',sep='\t', index=False)
    
    # Temporarily move the MS files - otherwise error with mzML files.
    for run in runs['MS_files'].tolist(): shutil.move(run,quant_dir)

    #Run Flash LFQ - MBR only when more than single run.
    log_file = f'{quant_dir}/quant_log.txt'
    logging.info(f'Running FlashLFQ - see log at {log_file}.')
    cmd = f'dotnet ./external_tools/FlashLFQ/CMD.dll --idt {quant_dir}/flashlfq_input.tsv --thr {int(os.cpu_count())} --chg --ppm 5 --ath --rep {quant_dir}'
    cmd += ' --mbr true' if len(runs) > 1 else ' --mbr false' 
    try:
        subprocess.run(f'{cmd} > {log_file} 2>&1', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"FlashLFQ failed: {e}")
        raise

    # Place the MS files back!
    for filename in runs['Filenames'].tolist():
        shutil.move(f'{quant_dir}/{filename}.raw', indir)

# Count quantifications per condition
def quant_condition(df: pd.DataFrame, quantification: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """
    Counts quantifications per condition (loaded annotation file).

    Parameters:
        df : pd.DataFrame
            A DataFrame containing peptide quantifications, where rows represent peptides and
            columns (runs) represent intensity values.
        
        quantification : pd.DataFrame
            A DataFrame with peptide information, including the 'Base Sequence' column that
            serves as a reference for peptide identification.
        
        runs : pd.DataFrame
            A DataFrame containing run metadata, including filenames and conditions for each run.
            Must include the 'Filenames' and 'condition' columns.

    Returns:
        pd.DataFrame
            A DataFrame with aggregated peptide quantifications. It includes an additional column 
            for each condition that counts the number of positive intensities per peptide for that 
            condition.
    """
    # Add and aggregate per peptide sequence - index for later merge.
    df['peptide'] = quantification['Base Sequence'].values
    df = df.groupby('peptide').sum(numeric_only=True).reset_index()

    # Now consider per condition if specified in annotation file.
    if 'condition' in runs:
        run2cond = dict(zip('Intensity_'+runs['Filenames'], runs['condition']))         #  Run to condition mapper.
        for cond in runs['condition']:
            cond_cols = [col for col in df.columns if run2cond.get(col) == cond]        # Intensity cols for this condition
            df[f'Intensity_{cond}'] = (df[cond_cols] > 0).sum(axis=1)                   # Count if non-zero intensity

    return df

# Read quantifications.
def read_quant(quant_file: str, runs: pd.DataFrame):
    """
    Read FlashLFQ output, will distinguish quantification with and without matching-
    between-runs (MBR). Also counts the number of quantifications per condition which
    is useful for later filtering of bacterial peptides.
    
    Parameters:
        quant_file : str
            The path to the quantification file, expected to be in tab-separated format.
        
        runs : pd.DataFrame
            A DataFrame containing metadata about the runs, including filenames and conditions.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]
            MBR : peptide quantifications with MBR.
            MSMS : peptide quantifications without considering MBR
    """

    logging.info('Appending FlashLFQ quantifications to peptide dataframe.')
    quantification = pd.read_csv(quant_file,sep='\t')

    # Get detection type and intensity columns.
    detection_cols = [col for col in quantification.columns if col.startswith('Detection Type')]
    intensity_cols = [col.replace('Detection Type', 'Intensity') for col in detection_cols]

    # Mapping dict run to name, if no name use filename.
    run2name = dict(zip('Intensity_' + runs['Filenames'], runs.get('name', runs['Filenames'])))

    # MSMS-only quantifications
    MSMS = quantification.copy()                              # Duplicate.
    for det_col, int_col in zip(detection_cols, intensity_cols):    # Set MBR detections to zero.
        MSMS.loc[MSMS[det_col] == 'MBR', int_col] = 0
    MSMS = MSMS[intensity_cols]                         # Intensities df only.
    MSMS = quant_condition(MSMS, quantification, runs)  # Add number of quantifications per condition.
    MSMS.rename(columns={col: f'MSMS_{run2name[col]}' for col in MSMS.columns if col in run2name}, inplace=True)    # Rename.
    MSMS = MSMS.sort_index(axis=1)

    # MBR quantifications
    MBR = quantification[intensity_cols].copy()
    MBR = quant_condition(MBR, quantification, runs)
    MBR.rename(columns={col: f'MBR_{run2name[col]}' for col in MBR.columns if col in run2name}, inplace=True)   # Rename.
    MBR = MBR.sort_index(axis=1)

    return MSMS, MBR

def pepform2seq(pepform: str) -> str:
    """
    Convert peptidoform string to plain peptide sequence.

    Parameters:
        pepform : str
            Peptidoform string.

    Returns:
        seq : str
            Plain peptide sequence.
    """

    pepform_noMod = re.sub(r'\[.*?\]', '', pepform)
    seq = re.sub(r'[^A-Z]', '', pepform_noMod)
    return seq

def sort_protein_list(protein_list: str) -> str:
    """
    Convert a protein list string to a sorted, semicolon-separated string, excluding decoy proteins.


    Parameters:
        protein_list : str
            A string representing a list of proteins, typically in the format:
            '[protein1, protein2, rev_protein3, ...]'.

    Returns:
        str
            A sorted, semicolon-separated string of proteins excluding decoy proteins.
    """
    
    list_string = protein_list.strip('[]').replace("'","").replace(', ', ';')
    sorted_proteins = sorted(prot for prot in list_string.split(';') if not prot.startswith('rev_'))
    return ';'.join(sorted_proteins)

# Read mokapot output basic - applies for peptide and PSM files.
def read_psms(args, engine: str, psms_all: pd.DataFrame, peptides_all: pd.DataFrame, prot_info: dict, instrument: str) -> pd.DataFrame:
    """
    Reads PSM and peptide identification results from the specified search engine output. Filter
    target peptides at 1% peptide FDR, and aggregate relevant information for downstream 
    analysis.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        engine : str
            Name of the search engine used (e.g., MSFragger, Comet).
        psms_all : pd.DataFrame
            DataFrame containing previously processed PSM-level results for concatenation.
        peptides_all : pd.DataFrame
            DataFrame containing previously processed peptide-level results for concatenation.
        prot_info : dict
            Dictionary containing gene, species, and description data for proteins based on FASTA parsing.
        instrument : str
            Instrument type (e.g., 'orbitrap', 'timsTOF') used for analysis.

    Returns:
        psm_all : pd.DataFrame
            Updated peptide-level data aggregated per search engine.
        peptides_all : pd.DataFrame
            Updated peptide-level data aggregated per search engine.
    """

    # Check for MS2Rescore PSM output.
    tsv_file = f'{args.in_dir}/report/ms2rescore/{engine}.psms.tsv'
    check_file(tsv_file)
    logging.info(f'Reading identifications from {engine}.')
    df = pd.read_csv(tsv_file,sep='\t')

    # Filter target PSMs both at 1% PSM- and peptide-level FDR.
    df = df[ (df['qvalue'] <= 0.01) &
             (df['meta:peptide_qvalue'] <= 0.01) &
             (~df['is_decoy'].astype(bool))]    # Filter target PSM at 1% PSM and Peptide Q-value.
    
    # Add columns
    df['peptide'] = df['peptidoform'].apply(pepform2seq) # Pepform to sequence.
    df['length'] = df['peptide'].apply(lambda x: len(x))
    df['engine'] = engine
    df['scan_id'] = df['run'] + '_' + df['spectrum_id'].astype(str)

    # Sort on lowest PSM q-value and, if identical, on the posterior error probability.
    df = df.sort_values(by=['qvalue', 'pep'], ascending=[True, True])

    # Store the PSM and Peptide level mokapot score, q-value and posterior error probabilities.
    df = df.rename(columns={'qvalue': 'psm_qval', 'score': 'psm_score', 'pep': 'psm_PEP',
        'meta:peptide_qvalue': 'peptide_qval', 'meta:peptide_score': 'peptide_score','meta:peptide_pep': 'peptide_PEP'})

    # Sort protein list and get genes, species and descriptions from info stored at FASTA reading.
    df['protein_list'] = df['protein_list'].apply(sort_protein_list)
    df['gene'] =  df['protein_list'].apply( lambda x: ';'.join([prot_info[prot_id]['gene'] for prot_id in x.split(';')]) )
    df['species'] = df['protein_list'].apply( lambda x: ';'.join(set(sorted([prot_info[prot_id]['species'] for prot_id in x.split(';')]))) )
    df['description leftmost'] = df['protein_list'].apply(lambda x: prot_info[ x.split(';')[0] ]['description'])

    # Specify whether defined as immunopeptide.
    if '-' in args.len:
        start,end = args.len.split('-')
        immunolengths = list(range(int(start),int(end)+1))
    else: immunolengths = [int(args.len)]
    df['immuno'] = df['peptide'].str.len().isin(immunolengths)

    # Filter with the required columns for downstream reporting (PSM and peptide level reports).
    columns = ['scan_id','run','peptide','length','immuno','peptidoform','protein_list','gene','species','description leftmost',
        'engine','psm_qval','psm_score','psm_PEP','peptide_qval','peptide_score','peptide_PEP','precursor_mz','retention_time','ion_mobility',
        'rescoring:spec_pearson','rescoring:rt_diff_best','rescoring:predicted_retention_time','rescoring:observed_retention_time']
    psms = df[columns]
   
    # Peptides per search engine: aggregate PSMs per each unique peptide sequence.
    agg_funcs = {
        'length': 'first',
        'immuno': 'first',
        'peptidoform': lambda x: ';'.join(pd.Series.unique(x)), # Can have multiple peptidoforms!
        'protein_list': 'first',
        'gene': 'first',
        'species': 'first',
        'description leftmost': 'first',
        'engine': 'first',
        'peptide_qval': 'min',
        'peptide_score': 'max',
        'peptide_PEP': 'min',
        'rescoring:spec_pearson': 'max',
        'rescoring:rt_diff_best': 'min'
    }
    peptides = psms.groupby('peptide').agg(agg_funcs).reset_index()
    
    # Concatenate and return the results.
    psms_all = pd.concat([psms_all,psms],ignore_index=True)
    peptides_all = pd.concat([peptides_all,peptides],ignore_index=True)
    
    return psms_all, peptides_all

# Count how many chimeric PSMs for a peptide sequence.
def count_chimera(chimera_psms) -> dict:
    """
    Count the number of ambiguous PSMs for each peptide and track their co-identified peptides.

    Parameters:
        chimera_psms : pd.DataFrame
            A DataFrame containing a 'peptide' column with semicolon-separated peptide sequences for each PSM.

    Returns:
        chimera_count : dict
            A dictionary listing the total number of ambiguous PSMs per peptide.
        chimera_pep : dict
            A dictionary listing the co-identified 'chimeric' peptides and their PSMs.
    """

    chimera_count = {}
    chimera_pep = {}

    for peptides in chimera_psms['peptide']:
        peptide_list = peptides.split(';')          # These peptides were ambiguously assigned to one spectrum.
        for peptide in peptide_list:
            if peptide not in chimera_count:        # Initiate count.
                chimera_count.update({peptide: 1})
                chimera_pep.update({peptide: {}})
            else:
                chimera_count[peptide]  += 1
            
            # Keep track of the chimeric hits.
            for peptide2 in peptide_list:
                if peptide == peptide2: continue
                else:
                    if peptide2 not in chimera_pep[peptide]:
                        chimera_pep[peptide].update({peptide2: 1})
                    else:
                        chimera_pep[peptide][peptide2] += 1
    
    # Restructure the dictionary.
    for peptide in chimera_pep:
        chimera_pep[peptide] = [f"{pep}({count})" for pep, count in chimera_pep[peptide].items()]

    return chimera_count, chimera_pep

# Read all search engine results and add to the scans dictionary.
def read_results(args, scans, runs, prot_info, instrument):
    """
    Reads and processes PSM and peptide identification results from MS2Rescore across different search engines. 
    Aggregates the PSMs and peptides, flags chimeric PSMs, integrates quantification results if available, 
    and appends PSM information to scan data.

    Parameters:
        args: argparse.Namespace
            Command-line arguments.
        scans: dict
            Dictionary storing scan-level information, indexed by scan ID.
        runs: pd.DataFrame
            DataFrame containing run metadata such as filenames, conditions, and experiment names.
        prot_info: dict
            Dictionary containing protein information including gene, species, and descriptions.
        instrument: str
            The MS instrument used, either orbitrap or timsTOF.

    Returns:
        peptides: pd.DataFrame
            DataFrame containing aggregated and de-duplicated peptide-level data, along with PSM counts and quantifications.
        psms: pd.DataFrame
            DataFrame containing aggregated PSM-level data across different search engines.
        scans: pd.DataFrame
            DataFrame of scan-level information with appended PSM identifications.
    """

    logging.info('Reading MS2Rescore PSM results (PSM and peptide Q-value < 1%).')
    indir = args.in_dir

    # Store available mokapot PSM results after rescoring.
    columns = ['scan_id','run','peptide','length','immuno','peptidoform','protein_list','gene','species','description leftmost',
        'engine','psm_qval','psm_score','psm_PEP','peptide_qval','peptide_score','peptide_PEP',
        'precursor_mz','retention_time','ion_mobility',
        'rescoring:spec_pearson','rescoring:rt_diff_best','rescoring:predicted_retention_time','rescoring:observed_retention_time']
    psms_all = pd.DataFrame(columns=columns)
    psms_all['immuno'] = psms_all['immuno'].astype(bool)

    # Store peptides.
    columns = ['peptide','length','immuno','peptidoform','protein_list','gene','species','description leftmost',
        'engine','peptide_qval','peptide_score','peptide_PEP','rescoring:spec_pearson','rescoring:rt_diff_best']
    peptides_all = pd.DataFrame(columns=columns)
    peptides_all['immuno'] = peptides_all['immuno'].astype(bool)

    # Read MS2Rescore PSM output for specified search engines.
    if args.frag: psms_all, peptides_all = read_psms(args, 'MSFragger', psms_all, peptides_all, prot_info, instrument)
    if args.comet: psms_all, peptides_all = read_psms(args, 'Comet', psms_all, peptides_all, prot_info, instrument)
    if args.sage: psms_all, peptides_all = read_psms(args, 'Sage', psms_all, peptides_all, prot_info, instrument)
    if args.peaks: psms_all, peptides_all = read_psms(args, 'PEAKS', psms_all, peptides_all, prot_info, instrument)

    # Aggregate PSMs across search engines.
    logging.info('Integrating PSMs from search engines.')
    agg_funcs = {
        'run': 'first', 
        'peptide': lambda x: ';'.join(pd.Series.unique(x)),     # Possibly a PSM is assigned to multiple peptides.
        'length': 'first',
        'immuno': 'first',
        'peptidoform': list,
        'protein_list': 'first',
        'gene': 'first',
        'species': 'first',
        'description leftmost': 'first',
        'engine': list,
        'psm_qval': list,
        'psm_score': list,
        'psm_PEP': list,
        'peptide_qval': list,
        'peptide_score': list,
        'peptide_PEP': list,
        'precursor_mz': 'first',
        'retention_time': 'first',
        'ion_mobility': 'first',
        'rescoring:spec_pearson': list,
        'rescoring:rt_diff_best': list,
        'rescoring:predicted_retention_time': list,
        'rescoring:observed_retention_time': list
    }
    psms_agg_engine = psms_all.groupby('scan_id').agg(agg_funcs).reset_index()
    
    # Subdivide ambiguous PSMs (scans with > 1 peptide) and unique PSMs (single matched peptide).
    chimera_psms = psms_agg_engine[psms_agg_engine['peptide'].str.contains(';')]
    chimera_count, chimera_pep = count_chimera(chimera_psms)
    chimera_psms.to_csv(f'{indir}/report/search_res/chimeric_PSMs.txt',sep='\t')
    psms = psms_agg_engine[~psms_agg_engine['peptide'].str.contains(';')].copy()
    logging.info(f'There were {len(psms)} PSMs to single peptides and {len(chimera_psms)} chimeric PSMs.')

    # Append run name and condition if specified.
    if 'condition' in runs.columns:
        psms = psms.merge(runs[['Filenames','name','condition']], left_on='run', right_on='Filenames', how='left')
        psms['run'] = psms['name']    # More informative than filename.
        psms = psms.drop(columns=['Filenames','name'])
    elif 'name' in runs.columns:
        psms = psms.merge(runs[['Filenames','name']], left_on='run', right_on='Filenames', how='left')
        psms['run'] = psms['name']    # More informative than filename.
        psms = psms.drop(columns=['Filenames','name'])
    
    # Aggregate now the peptides across search engines.
    logging.info('Integrating peptide identifications from search engines.')
    agg_funcs = {
        'length': 'first',
        'immuno': 'first',
        'protein_list': 'first',
        'gene': 'first',
        'species': 'first',
        'description leftmost': 'first',
        'peptidoform': lambda x: ';'.join(pd.Series.unique(x)),
        'engine': list,
        'peptide_qval': list,
        'peptide_score': list,
        'peptide_PEP': list,
        'rescoring:spec_pearson': 'max',
        'rescoring:rt_diff_best': 'min'
    }
    peptides = peptides_all.groupby('peptide').agg(agg_funcs).reset_index()

    # De-duplicate peptidoforms, count engine/peptidoform, and get lowest q-value, PEP and highest score
    peptides['peptidoform'] = peptides['peptidoform'].apply(lambda x: ';'.join(list(set(x.split(';')))))
    peptides['#peptidoform'] = peptides['peptidoform'].apply(lambda x: x.count(';') + 1)
    peptides['#engine'] = peptides['engine'].apply(len)
    peptides['peptide_qval_lowest'] = peptides['peptide_qval'].apply(min)
    peptides['peptide_PEP_lowest'] = peptides['peptide_PEP'].apply(min)
    peptides['peptide_score_highest'] = peptides['peptide_score'].apply(max)

    # Flag chimeric PSMs.
    peptides['#chimeric_PSMs'] = peptides['peptide'].map(chimera_count).fillna(0).astype(int)
    peptides['#chimeric_peptides'] = peptides['peptide'].map(chimera_pep).fillna('NA').astype(str)
    
    # Append the total number of PSMs.
    psms_total = psms.groupby(['peptide']).size().reset_index(name='PSMs_total')
    psms_total['PSMs_total'] = psms_total['PSMs_total'].astype(int)
    peptides = peptides.merge(psms_total,left_on='peptide',right_on='peptide',how='left')
    
    # Append the number of PSMs per run.
    psms_run = psms.groupby(['peptide', 'run']).size().unstack(fill_value=0)
    psms_run.columns = ['PSMs_run_' + str(col) for col in psms_run.columns]
    peptides = peptides.merge(psms_run,left_on='peptide',right_index=True,how='left')
    
    # Append the number of PSMs per condition.
    if 'condition' in psms.columns:
        psms_cond = psms.groupby(['peptide', 'condition']).size().unstack(fill_value=0)
        psms_cond.columns = ['PSMs_cond_' + str(col) for col in psms_cond.columns]
        peptides = peptides.merge(psms_cond,left_on='peptide',right_index=True,how='left')

    # Filter out the peptides only identified by chimeric PSMs!
    chimera_peptides = peptides[(peptides['PSMs_total'] <= 0) | (peptides['PSMs_total'].isna())]
    chimera_peptides.to_csv(f'{indir}/report/search_res/chimera_only_peptides.txt',sep='\t',index=False)
    peptides = peptides[peptides['PSMs_total'] > 0]
    logging.info(f'There were {len(peptides)} unique peptide sequences.')

    # Join quantification results if any.
    if os.path.exists(f'{indir}/report/quant/QuantifiedPeptides.tsv'):
        MSMS, MBR = read_quant(f'{indir}/report/quant/QuantifiedPeptides.tsv', runs)
        peptides = peptides.merge(MSMS,on='peptide', how='left')  #MSMS only quantifications
        if MBR.shape[1] > 1:  # No MBR ran if only one input file.
            peptides = peptides.merge(MBR,on='peptide', how='left')

    # Add PSM information to scans.
    logging.info('Appending PSM info the scan dataframe')
    scans = pd.DataFrame.from_dict(scans, orient='index')
    scans = scans.merge(psms_agg_engine[['scan_id','peptide','peptidoform','length','immuno','engine','protein_list',
        'species','gene','rescoring:spec_pearson','rescoring:rt_diff_best']],left_index=True,right_on='scan_id',how='left')
    
    return peptides, psms, scans

def histogram_plotter(df: pd.DataFrame, plot_folder: str, group: str) -> None:
    """
    Plots a histogram of peptide lengths from the given DataFrame.

    Parameters:
        df: pd.DataFrame
            DataFrame containing peptide data, which should include a 'length' column and an 'immuno' column.
        plot_folder: str
            Directory where the plots will be saved.
        group: str
            Label for the group of peptides being plotted, used in the title and filenames.

    Returns:
        None
    """

    # Prepare plotting data. 
    peptide_lengths = range(df['length'].min(), (df['length'].max()+1))  # Plot from 7 to 20.
    peptide_counts = df['length'].value_counts().reindex(peptide_lengths, fill_value=0) # Initiate counts for barplot.
    colors = ['firebrick' if df[df['length'] == length]['immuno'].any() else 'gray' for length in peptide_counts.index]    # Specify immuno in red.
    total_peptides = len(df)
    immuno_peptides = df['immuno'].sum()

    # Actual plot.
    sns.set(style="ticks")
    plt.figure(figsize=(8, 6))
    sns.barplot(x=peptide_counts.index, y=peptide_counts.values, hue=peptide_lengths, legend=False, palette=colors)
    plt.xlabel('Peptide Length')
    plt.ylabel('Peptides')
    plt.title(f'Length histogram for {group} peptides\n'
              f'{immuno_peptides}/{total_peptides} ({immuno_peptides/total_peptides:.2%}) were immunopeptides.')
    plt.xticks(rotation=45)

    # Save as PNG and SVG.
    plt.savefig(f'{plot_folder}/{group}_histogram.png',dpi=600)
    plt.savefig(f'{plot_folder}/{group}_histogram.svg',dpi=600)
    plt.close()

def histogram_plot(args, peptides: pd.DataFrame):
    """
    Histogram plots for all peptides and per species.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        peptides : pd.DataFrame
            DataFrame containing peptide data.

    Returns:
        None
    """

    logging.info('Plotting peptide length histogram plot.')

    # Plot all peptides regardless of species.
    plot_folder = f'{args.in_dir}/report/plots/len_histogram/'
    make_folder(plot_folder)
    histogram_plotter(peptides, plot_folder, 'all')

    # Plot for different species if present:
    if peptides['species'].nunique() > 1:
        for species in peptides['species'].unique():
            peptides_species = peptides[peptides['species'] == species]
            histogram_plotter(peptides_species, plot_folder, species)

def id_per_run_plotter(data: pd.DataFrame, ylabel: str, title: str, file_suffix: str, group: str, plot_folder: str) -> None:
    """
    Generates and saves a stacked bar plot from the given data.

    Parameters:
        data : pd.DataFrame
            DataFrame containing the data to be plotted. The index should represent 
            different runs and the columns should represent binding levels.
        ylabel : str
            Label for the y-axis of the plot.
        title : str
            Title for the plot.
        file_suffix : str
            Suffix to be appended to the output file names.
        plot_folder : str
            Directory path for saving the output plots and summary files.
        group : str
            Group identifier for naming the output files.

    Returns:
        None
    """

    colors = {'SB': 'firebrick', 'WB': 'dodgerblue', 'NB': 'gray', 'NA': 'gainsboro'}
    color_list = [colors[col] for col in data.columns]
    data.index = data.index.str.replace('PSMs_run_', '')
    data.to_csv(f'{plot_folder}/{file_suffix}.txt', sep='\t')
    data.plot(kind='bar', stacked=True, color=color_list, edgecolor='black')
    plt.xlabel('Runs')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(title='Binding Level')
    plt.tight_layout()
    plt.savefig(f'{plot_folder}/{group}_{file_suffix}.png', dpi=600)
    plt.savefig(f'{plot_folder}/{group}_{file_suffix}.svg', dpi=600)
    plt.clf()
    plt.close()

def id_per_run(df: pd.DataFrame, plot_folder: str, group: str) -> None:
    """
    Make stacked bar plots for PSMs and peptides identified per run, indicating binding levels.

    Parameters:
        df : pd.DataFrame
            DataFrame containing PSM data with binding levels and PSM counts per run.
        plot_folder : str
            Directory path for saving the output plots and summary files.
        group : str
            Group identifier for naming the output files.

    Returns:
        None
    """
    
    # If NetMHCpan was not ran: initiate here this column.
    df['bind_level'] = df.get('bind_level', 'NA')

    # Get PSM columns.
    psm_columns = [col for col in df.columns if col.startswith('PSMs_run')]

    # Count the number of PSMs and peptides per bind level.
    psm_data = pd.DataFrame(index=psm_columns, columns=['SB', 'WB', 'NB', 'NA'], data=0)
    pep_data = pd.DataFrame(index=psm_columns, columns=['SB', 'WB', 'NB', 'NA'], data=0)
    for index, row in df.iterrows():
        for psm_col in psm_columns:
            if row[psm_col] > 0:
                psm_data.loc[psm_col, row['bind_level']] += row[psm_col]    # Sum up PSMs.
                pep_data.loc[psm_col, row['bind_level']] += 1               # Add one as peptide is identified.

    # Plot PSMs and peptides
    id_per_run_plotter(psm_data, 'Number of PSMs', 'PSMs per run', 'psms_per_run', group, plot_folder)
    id_per_run_plotter(pep_data, 'Number of Peptides', 'Peptides identified per run', 'peps_per_run', group, plot_folder)

def numbers_per_run(args, peptides: pd.DataFrame) -> None:
    """
    Make stacked bar plots for PSMs and peptides identified per run, indicating binding levels.

    Parameters:
        args : argparse.Namespace
            Command-line arguments.
        peptides : pd.DataFrame
            The main peptide dataframe from the ms2rescore results.

    Returns:
        None
    """

    # Plot all PSMs/peptides per run regardless of species.
    logging.info('Plotting PSMs/peptides per run.')
    plot_folder = f'{args.in_dir}/report/plots/identifications_per_run/'
    make_folder(plot_folder)
    id_per_run(peptides, plot_folder, 'all')

    # Plot for different species if present:
    if peptides['species'].nunique() > 1:
        for species in peptides['species'].unique():
            peptides_species = peptides[peptides['species'] == species]
            id_per_run(peptides_species, plot_folder, species)

# Get best ranked binder
def get_best_binder(netmhcpan, netmhcpan_folder, mhc_binding, tool):

    # Get the best ranked HLA per peptide
    netmhcpan = netmhcpan.sort_values(by='rank',ascending=True)
    netmhcpan.to_csv(f'{netmhcpan_folder}parsed_output.txt',sep='\t',index=True)
    best_ranked = netmhcpan.groupby('peptide').agg({
        'mhc': 'first',
        'core': 'first',
        'rank': 'first',
        'bind_level': 'first'
        }).reset_index().fillna('NB')   #No binder.

    # Add tool information.
    if 'II' in tool: best_ranked['tool'] = f'{tool}-4.3'
    else: best_ranked['tool'] = f'{tool}-4.1'

    # Append to mhc_binding dictionary
    mhc_binding.update( best_ranked.set_index('peptide').T.to_dict() )

    return mhc_binding

# Run netMHCpan-4.1 or netMHCIIpan-4.3
def netMHCpan(args, peptides, mhc_binding, tool):
    netmhcpan_folder = f'{args.in_dir}/report/{tool}/'
    make_folder(netmhcpan_folder)

    # Get HLA alleles to check.
    if f'{tool}_{args.HLA}' not in hla_dict:
        logging.warning(f'Specified HLA flag "{args.HLA}" not in {tool} dictionary! Skipping {tool}!')
        return peptides
    else: HLA = hla_dict[f'{tool}_{args.HLA}']

    # Input peptide file - adjust length to HLA tool. Class I: 8 to 12, class II: 13 to 18.
    lengths = range(8,13) if tool == 'netMHCpan' else range(13,19)
    input_peptide = peptides.loc[ peptides['length'].isin(lengths), 'peptide' ]
    if len(input_peptide) == 0:
        logging.warning(f'No input peptides of length {lengths}! Skipping {tool}')
        return peptides
    else: input_peptide.to_csv(f'{netmhcpan_folder}input.txt',index=False,header=False)

    # Run
    if tool == 'netMHCpan':
        cmd = f'external_tools/netMHCpan-4.1/netMHCpan -p {netmhcpan_folder}input.txt -l 8,9,10,11,12 -s -a {HLA} > {netmhcpan_folder}output.txt'
    else:
        cmd = f'external_tools/netMHCIIpan-4.3/netMHCIIpan -inptype 1 -f {netmhcpan_folder}input.txt -s -a {HLA} > {netmhcpan_folder}output.txt'
    logging.info(f'Running {tool} for {len(input_peptide)} peptides ({lengths[0]} to {lengths[-1]} mers) for {HLA} alleles.')
    os.system(cmd)

    # Read and store netMHCpan output: need to filter specific lines.
    with open(f'{netmhcpan_folder}output.txt', 'r') as out_file: lines = out_file.readlines()
    if tool == 'netMHCpan':
        filtered_data = ''.join([line for line in lines if ('    PEPLIST' in line)])    # Actual results
        columns = ['Pos','mhc','peptide','core','Of','Gp','Gl','Ip','Il','Icore','Identity','score_EL','rank','arrow','bind_level']
    else:
        filtered_data = ''.join([line for line in lines if ('    Sequence  ' in line)])    # Actual results
        columns = ['Pos','mhc','peptide','Of','core','Core_Rel','Inverted','Identity','score_EL','rank','Exp_Bind','arrow','bind_level']

    # Get the best ranked HLA per peptide and update the binding dictionary
    netmhcpan = pd.read_csv(StringIO(filtered_data),header=None,names=columns, delim_whitespace=True)
    get_best_binder(netmhcpan, netmhcpan_folder, mhc_binding, tool)
    
    return mhc_binding

# Run NetMHCpan for MHC class I and/or MHC class II
def run_netMHCpan(args, peptides, scans):

    # MHC binding dict: store results here and append together
    mhc_binding = {}

    # MHC class I : netMHCpan-4.1
    mhc_binding = netMHCpan(args, peptides, mhc_binding, 'netMHCpan')
    
    # MHC class II : netMHCIIpan-4.3
    mhc_binding = netMHCpan(args, peptides, mhc_binding, 'netMHCIIpan')

    # Rare case if no binders, initiate 'NA' columns for downstream compatibility.
    if len(mhc_binding) > 0:
        for key in ['tool','mhc','rank','bind_level','core']:
            peptides[key] = peptides['peptide'].apply(lambda x: mhc_binding[x][key] if x in mhc_binding else 'NA')

    return peptides, scans

def pie_plot(column, name, ax, color_palette=None) -> None:
    unique, counts = np.unique(column, return_counts=True)
    # Check if colors passed
    if color_palette:
        colors = [color_palette.get(key, '#cccccc') for key in unique]
    else:
        colors = None
    ax.pie(counts, labels=unique, autopct='%1.1f%%',colors=colors)
    ax.set_title(f'{name} (N={len(column)})', fontsize="small")
    plt.tight_layout()

def netMHCpan_barchart(peptide_counts, ax, colors, title) -> None:
    peptide_counts.plot(kind='bar', stacked=True, color=[colors[level] for level in peptide_counts.columns],edgecolor='black',ax=ax)
    ax.set_xlabel('Peptide Length')
    ax.set_ylabel('Peptides')
    ax.set_title(title)
    plt.tight_layout()

# This will put a main overview of binding strength per peptide length.
def netMHCpan_overview(peptides, group, plot_folder) -> None:

    # Initial check.
    if 'mhc' not in peptides or len(peptides) < 5:
        logging.warning(f'Too few (<5) netMHCpan binding results found for {group} (species) peptides. Skipping MHC overview.')
        return
    
    # Count peptide lengths to be plotted.
    peptide_lengths = range(peptides['length'].min(), peptides['length'].max()+1)
    peptide_counts = peptides.groupby(['length', 'bind_level']).size().unstack(fill_value=0).reindex(peptide_lengths, fill_value=0)
    colors = {'SB': 'firebrick', 'WB': 'dodgerblue', 'NB': 'gray', 'NA': 'gainsboro'}

    # Adjust the order of the bars (SB -> WB -> NB, bottom to top)
    existing_levels = peptide_counts.columns.tolist()  # Get existing columns
    plot_levels = [level for level in ['SB', 'WB', 'NB','NA'] if level in existing_levels]  # Filter for existing levels
    peptide_counts = peptide_counts[plot_levels]

    # Stacked barplot on top, piecharts below for immunopeptides.
    sns.set(style="ticks")
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])  # First row, span across all columns
    ax2 = fig.add_subplot(gs[1, 0])  # Second row, first column
    ax3 = fig.add_subplot(gs[1, 1])  # Second row, second colum
   
    # Stacked barplot: binders per peptide length.
    netMHCpan_barchart(peptide_counts, ax1, colors, f'netMHCpan binding overview for {group} peptides')
        
    # Predicted binding strength for immunopeptides.
    df1 = peptides[peptides['immuno']]
    pie_plot(df1['bind_level'], f'Predicted binding strength immunopeptides', ax2)
        
    # Best ranked MHC for SB
    df2 = df1[df1['bind_level'].isin(['SB', 'WB'])]
    pie_plot(df2['mhc'], f'Best MHC for SB/WBs', ax3)
    
    # Save as PNG and SVG
    plt.tight_layout()
    plt.savefig(f'{plot_folder}/{group}_peptides.png', dpi=600)
    plt.savefig(f'{plot_folder}/{group}_peptides.svg', dpi=600)
    plt.clf()
    plt.close()

# Plot MHC binding logos for SB/WB binding level peptides.
def netMHCpan_logos(peptides, group, plot_folder):

    # First get the number of MHCs
    binders = peptides[peptides['bind_level'].isin(['SB', 'WB'])]
    mhcs = binders['mhc'].unique()

    # Plot per MHC a single plot with the sequence logo on the left and piechart of the peptide length(s) on the right.
    for mhc in mhcs:
        # Filter SB/WB for MHC and check if sufficient peptides for plot.
        df = binders[binders['mhc'] == mhc]
        if len(df) < 5:
            logging.warning(f'Too few (<5) netMHCpan SB/WB binding peptides for {mhc} logo plot.')
            continue

        # Initiate figure.
        fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1]})

        # Left subplot: MHC logo plot.
        logo_plot(df['core'].astype(str).tolist(), mhc, axes[0])

        # Right subplot: Bar chart of peptide length(s). 
        # Initiate peptide count and define color palette for binding strength.
        peptide_lengths = range(df['length'].min(),df['length'].max() + 1)
        peptide_counts = df.groupby(['length', 'bind_level']).size().unstack(fill_value=0).reindex(peptide_lengths, fill_value=0)
        colors = {'SB': 'firebrick', 'WB': 'dodgerblue'}

        # Adjust the order of the bars (SB -> WB -> NB, bottom to top)
        existing_levels = peptide_counts.columns.tolist()  # Get existing columns
        plot_levels = [level for level in ['SB', 'WB'] if level in existing_levels]  # Filter for existing levels
        peptide_counts = peptide_counts[plot_levels]
        
        # Make the barchart now.
        netMHCpan_barchart(peptide_counts, axes[1], colors, 'Length histogram')

        #Save figure.
        mhc = mhc.replace('*', '').replace(':', '_')   # Can't use asterisk or colon in filename, e.g. HLA-A*01:01.
        plt.tight_layout()
        plt.savefig(f'{plot_folder}/{mhc}.png', dpi=600)
        plt.savefig(f'{plot_folder}/{mhc}.svg', dpi=600)
        plt.clf()
        plt.close()

# Plot netMHCpan main overview and peptide sequence logos.
def netMHCpan_plot(peptides, args):
    logging.info('Plotting netMHCpan sequence motifs and overview plots.')

    # Folders
    main_folder = f'{args.in_dir}/report/plots/netMHCpan/'
    make_folder(main_folder)
    logo_folder = f'{args.in_dir}/report/plots/netMHCpan/mhc_logos/'
    make_folder(logo_folder)

    # Plot overview and logos for all peptides.
    netMHCpan_overview(peptides, 'all', main_folder)
    netMHCpan_logos(peptides, 'all', logo_folder)

    # Plot for different species if present:
    if peptides['species'].nunique() > 1:
        for species in peptides['species'].unique():
            peptides_species = peptides[peptides['species'] == species]
            netMHCpan_overview(peptides_species, species, main_folder)

# Plot KLD.
def KLD_plot(kld_file, args):
    
    # Make Gibbs plot folder.
    plot_folder = f'{args.in_dir}/report/plots/gibbs/'
    make_folder(plot_folder)

    # Plot the stacked bar plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot each column as a stack
    bottom = pd.Series([0]*len(kld_file), index=kld_file.index)
    for col in kld_file.columns:
        ax.bar(kld_file.index, kld_file[col], bottom=bottom, label=f'{col}', edgecolor='black')
        bottom += kld_file[col]

    # Add labels and title
    ax.set_xticks(kld_file.index)
    ax.set_xticklabels(kld_file.index)
    ax.set_xlabel('Number of Clusters')
    ax.set_ylabel('Kullbach Leibler Distance (KLD)')
    ax.set_title('GibbsCluster KLD plot')
    ax.legend(title='Clusters')
    plt.tight_layout()
    
    # Save as PNG and SVG
    plt.savefig(f'{plot_folder}/KLD_clusters.png', dpi=600)
    plt.savefig(f'{plot_folder}/KLD_clusters.svg', dpi=600)
    plt.clf()
    plt.close()

    return

# GibbsCluster2.0 is ran for immunopeptides
def run_Gibbs(args, peptides):

    # Prepare folders
    gibbs_folder = f'{args.in_dir}/report/gibbs/'
    make_folder(gibbs_folder)

    # Filter to set length
    immunopeptides = peptides[peptides['immuno']]
    peptide_lengths = immunopeptides['length'].unique()
    immunopeptides['peptide'].to_csv(f'{gibbs_folder}/peptide_input.txt',header=False,index=False)

    # Construct command - defaults see webserver, different for MHC class I (single len/multiple) and class II!
    gibbs_clusters = '[2-5]' if args.gibbs == 'auto' else f'[{clusters}]'
    cmd = (f'perl external_tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl -H /usr/bin/R -T -j 2 -l 9 -k 32 -S 5 -P GibbsCluster2'
          f' -g {gibbs_clusters}'                   # Gibbs clusters - user parameter
          f' -f {gibbs_folder}/peptide_input.txt'   # Peptide input file
          f' -R {gibbs_folder}')

    # MHC class I max length is 12.
    if max(peptide_lengths) < 13:
        cmd += ' -C'    # Make clustering moves at each iteration (not for MHC class II)
        if len(peptide_lengths) > 1:
            cmd += ' -D 4 -I 1' # Maximum deletion and insertion length, respectively.

    # Run
    logging.info(f'Running GibbsCluster2.0 for {str(len(immunopeptides))} peptides of specified immunolength(s) {args.len}.')
    os.system(cmd)

    # Move generated Gibbs folder (trailing integer) to output subfolder (replace if anything there)
    if os.path.exists(f'{gibbs_folder}/output/'): shutil.rmtree(f'{gibbs_folder}/output/')
    for folder in glob.glob(f'{gibbs_folder}/GibbsCluster2*'):
        shutil.move(folder,f'{gibbs_folder}/output/')

    # Read KLD info, print its plot and determine result cluster.
    kld_file = pd.read_csv(f'{gibbs_folder}/output/images/gibbs.KLDvsClusters.tab',sep='\t')
    KLD_plot(kld_file, args)
    row_sums = kld_file.sum(axis=1)     # Sum up KLD per row.
    final_cluster = row_sums.idxmax()   # Winning cluster = highest KLD sum.

    # Read clusters.
    logging.info(f'Will read results when using {final_cluster} clusters (KLD: {row_sums.loc[final_cluster]:.3f}).')
    with open(f'{gibbs_folder}/output/res/gibbs.{final_cluster}g.out', 'r') as out_file: lines = out_file.readlines()
    filtered_data = ''.join([line for line in lines if ('    Peplist' in line)])    # Actual results
    columns = ['G','cluster','nr0','peptide','gibbs_core','o','nr1','ip','nr2','il','nr3','dp','nr4','dl','nrx','peplist','sS','nr5','bgG','nr6','bgS','nr7','cS','nr8']

    # Get the best ranked HLA per peptide and update the binding dictionary
    gibbs = pd.read_csv(StringIO(filtered_data),header=None,names=columns, delim_whitespace=True)
    gibbs = gibbs[['peptide','cluster','gibbs_core']]
    peptides = peptides.merge(gibbs,on='peptide',how='left')
    peptides.fillna('NA', inplace=True)
    return peptides

# Make logo with logomaker for a peptide input (netMHCpan / Gibbs cores)
def logo_plot(peptide_list, name, ax):
    counts_mat = lm.alignment_to_matrix(peptide_list)
    info_mat = lm.transform_matrix(counts_mat,from_type='counts',to_type='information')
    logo = lm.Logo(info_mat,color_scheme="weblogo_protein",ax=ax)
    ax.set_xlabel('Position',fontsize=14)
    ax.set_ylabel("Bits", labelpad=-1,fontsize=14)
    ax.set_title(f'{name} (N={str(len(peptide_list))})')
    plt.tight_layout()
    
# Plot gibbs
def gibbs_plot(args, peptides) -> None:

    # Get the ones witihin a Gibbs cluster (only for immunolengths, non-trash clusters)
    clustered = peptides[peptides['cluster'] != 'NA']
    num_clusters = clustered['cluster'].nunique()
    rows = num_clusters
   
    # If HLA binding results present - we plot the % SB/WB and their MHC distribution (3 cols).
    if 'mhc' in clustered:
        fig, axes = plt.subplots(rows, 3, figsize=(15, 4 * rows), gridspec_kw={'width_ratios': [2, 1, 1]})
        unique_mhc_values = clustered['mhc'].dropna().unique()
        colors = sns.color_palette('hsv', len(unique_mhc_values)).as_hex()
        mhc_colors = dict(zip(unique_mhc_values, colors))
        bind_colors = {'SB': 'firebrick', 'WB': 'dodgerblue', 'NB': 'gray', 'NA': 'lightgrey'}
    else:
        fig, axes = plt.subplots(rows, 1, figsize=(10, 4 * rows))

    if rows == 1: axes = [axes]  # Ensure axes is always a list of lists, even if there's only one row

    # Get unique clusters in ascending manner
    unique_clusters = sorted(clustered['cluster'].astype(int).unique())

    # Plot them row by row on the graphic axes.
    for i, cluster in enumerate(unique_clusters):
        df = clustered[clustered['cluster'] == cluster]

        if 'mhc' in df.columns:
            if df['bind_level'].isin(['SB', 'WB']).any(): current_axes = axes[i]
            else: current_axes = axes[i][:2]

            logo_plot(df['gibbs_core'].astype(str).tolist(), f'Cluster {int(cluster)}', current_axes[0])
            pie_plot(df['bind_level'], f'Predicted binding strength', current_axes[1], bind_colors)

            # MHC distribution plot
            if df['bind_level'].isin(['SB', 'WB']).any():
                filtered_df = df[df['bind_level'] != 'NB']
                pie_plot(filtered_df['mhc'], f'MHC alleles SB/WB', current_axes[2], mhc_colors)
        
        # No HLA binding results.
        else:
            current_axes = [axes[i]]
            logo_plot(df['gibbs_core'].astype(str).tolist(), f'Cluster {int(cluster)}', current_axes[0])

    plt.tight_layout()
    plt.savefig(f'{args.in_dir}/report/plots/gibbs/All_Clusters.png', dpi=600)
    plt.savefig(f'{args.in_dir}/report/plots/gibbs/All_Clusters.svg', dpi=600)
    plt.clf()
    plt.close()

# Venn Diagram search engines
def run_Venn(engine_dict, plot_folder, name):
    
    # Need to run set for the search engine peptide list.
    for engine in engine_dict:
        engine_dict[engine] = set(engine_dict[engine])
    venny4py(sets=engine_dict,out=plot_folder)
    shutil.move(f'{plot_folder}/Venn_{str(len(engine_dict))}.png',f'{plot_folder}/{name}_peptides.png')

    plt.clf()
    plt.close()

# Parent function Venn diagrams.
def VennOverlap(peptides, args):
    logging.info('Plotting peptide overlap and barplot per search engine.')

    # Prepare folders.
    plot_folder = f'{args.in_dir}/report/plots/engine_overlaps/'
    make_folder(plot_folder)

    # Prepare the dictionaries - and check for more than single search engine.
    unique_engines = set(engine for sublist in peptides['engine'] for engine in sublist)
    if len(unique_engines) == 1:
        logging.info('Results only found for one search engine - skipping VennOverlap plots!')
        shutil.rmtree(plot_folder)
        return

    logging.info('Plotting Venn diagram of search engine results.')
    engine_dict_all = {engine: [] for engine in unique_engines}     # Overlap all peptides
    engine_dict_immuno = {engine: [] for engine in unique_engines}  # Overlap immunopeptides.
    engine_dict_SB = {engine: [] for engine in unique_engines}      # Overlap for immunopeptides - strong binders.
    
    # Iterate through the DataFrame and populate the dictionaries
    for index, row in peptides.iterrows():
        for engine in row['engine']:
            engine_dict_all[engine].append(row['peptide'])
            if row['immuno']:
                engine_dict_immuno[engine].append(row['peptide'])
                if 'bind_level' in row and row['bind_level'] == 'SB':
                    engine_dict_SB[engine].append(row['peptide'])

    # Make the Venn diagrams.
    run_Venn(engine_dict_all, plot_folder, 'all')
    run_Venn(engine_dict_immuno, plot_folder, 'immuno')
    if len(engine_dict_SB) > 0:
        run_Venn(engine_dict_SB, plot_folder, 'immuno_SB')

# BLASTP to extended host proteomes for double check sequence ambiguity (Ile <-> Leu).
def run_BLASTP(args, peptides):
    
    # Check whether host in pre-compiled BLAST databases.
    host = args.host.upper()
    blastdb = f'external_tools/blastdb/{host.lower()}'
    if not os.path.exists(f'{blastdb}.pin'):
        logging.warning(f'\nSpecied host species BLAST db {blastdb} not found!'
                        f'\nEither correct or generate BLAST protein database at path mentioned above. Skipping background check now')
        return

    # Write FASTA of non-host/contaminant peptides.
    nonhost_peptides = peptides.loc[
        (~peptides['species'].astype(str).str.contains(host)) &
        (~peptides['peptide'].astype(str).str.contains('contaminant')),
        'peptide'
    ].to_list()

    # Check number of non-host peptides
    if len(nonhost_peptides) == 0:
        logging.info('No non-host peptides found for BLASTP background check, skipping!')
        return peptides
    else:
        logging.info(f'Matching {str(len(nonhost_peptides))} peptide sequences to background {host} FASTA or contaminants.')
    
    # Write out tmp FASTA for query.
    with open('blastp_in.fasta', 'w') as outfile:
        for pep in nonhost_peptides:
            outfile.write(f'>{pep}\n{pep.replace("I","L")}\n')  # Replace I to L in queried peptide sequence!
    
    # Perform BLASTP
    cmd = f'blastp -task blastp-short -db {blastdb} -query blastp_in.fasta -out blastp_out.tsv -outfmt "6 qacc sacc qlen nident sseq"'
    os.system(cmd)

    # Read results and calculate+sort on the % identity - top hit will be stored per peptide.
    results = pd.read_csv('blastp_out.tsv',sep='\t',names=['peptide','BLASTP_match','qlen','nident','BLASTP_matchedSeq'])
    results['BLASTP_ident%'] = (results['nident'] / results['qlen'] * 100).round(2)
    results = results[['peptide', 'BLASTP_ident%', 'BLASTP_match', 'BLASTP_matchedSeq']] \
           .sort_values('BLASTP_ident%', ascending=False) \
           .groupby('peptide') \
           .first()

    # Add BLASTP output as novel columns to the peptides df.
    peptides = peptides.merge(results,left_on='peptide',right_index=True,how='left')
    peptides.fillna({'BLASTP_ident%': '0', 'BLASTP_match': 'NA', 'BLASTP_matchedSeq': 'NA'}, inplace=True)
    
    # Clean up
    os.remove('blastp_out.tsv')
    os.remove('blastp_in.fasta')

    # Adjust species if 100% identical match.
    peptides['species'] = peptides.apply(lambda row: row['species'] + f';{host}' if row['BLASTP_ident%'] == 100 else row['species'], axis=1)
    peptides['species'] = peptides['species'].apply( lambda x: ';'.join(sorted(x.split(';'))) )         # Sort species so compatible with other species rows!
    
    return peptides

# Transform intensities to % (divide maximum*100).
def transform_int(array) -> list:
    intensities = array.tolist()
    scaled = [(x/max(intensities))*100 for x in intensities]    # Relative intensity (%)
    return scaled

# Annotate ions
def annotate_ion(mz, fragdict):
    for ion, mz_ion in fragdict.items():
        if abs(mz_ion - mz) < 0.02:     # 0.02 Da tolerance.
            return ion
    return 'NA'

# Remove precursor ion with 1 Da tolerance.
def remove_precursor(mz_array,int_array,prec_mz):
    mask = np.abs(mz_array - prec_mz) > 1

    return mz_array[mask], int_array[mask]

# Return ion colors dependent on a/b/y and neutral loss or not.
def ion_color(ion):
    if ion == 'NA': return 'lightgrey'
    elif ion.startswith('a'):
        if 'H2O' in ion: return 'mediumseagreen'
        elif 'NH3' in ion: return 'yellowgreen'
        else: return 'green'
    elif ion.startswith('b'):
        if 'H2O' in ion: return 'darkturquoise'
        elif 'NH3' in ion: return 'steelblue'
        else: return 'darkblue'
    elif ion.startswith('y'):
        if 'H2O' in ion: return 'orangered'
        elif 'NH3' in ion: return 'palevioletred'
        else: return 'firebrick'

# Superscript ends.
def format_ion_text(ion_text):
    parts = ion_text.split('+', 1)
    if len(parts) == 1:
        # No '+' found, return the text as is
        return f"${ion_text}$"
    else:
        base = parts[0]
        superscript = '+' + parts[1] if parts[1] else '+'
        return f"${base}^{{{superscript}}}$"

# Plot spectra for non-host and non-contaminant peptides.
def plot_spectra(args, scans, instrument):
    
    # First check how many spectra and unique peptidoforms.
    scans['species'] = scans['species'].fillna('contaminant')       # The scans won't be returned - just to filter here next line.
    
    # To plot unambiguous assigned spectra.
    if args.plot_chimera:
        df = scans[
                (~scans['species'].str.contains(f'contaminant|{args.host.upper()}', regex=True)) |  # Non-host PSMs.
                (scans['peptide'].str.contains(';'))    # Also include chimeric PSMs.
                ]
    else:
        df = scans[~scans['species'].str.contains(f'contaminant|{args.host.upper()}', regex=True)]  # Non-host PSMs.
    
    if len(df) == 0:
        logging.warning('No non-host/contaminant matching spectra, will not plot annotated spectra.')
        return
    
    # Get peptidoforms and make folders
    all_pepforms = [pepform for sublist in df['peptidoform'] for pepform in sublist]
    uniq_pepforms = list(set(all_pepforms))
    logging.info(f'Plotting {len(df)} spectra matching {len(uniq_pepforms)} non-host/contaminant peptidoforms.')
    plot_folder = f'{args.in_dir}/report/plots/spectra/'
    make_folder(plot_folder)

    # Fragment dictionary per peptidoform for plotting - using pyteomics ProForma.
    fragdict = {}
    loss_masses = {
        'H2O': 18.01528,    # Water loss.
        'NH3': 17.02655     # Ammonia loss.
    }
    for pepform in uniq_pepforms:
        modpep = pepform.split('/')[0].replace('UNIMOD:','')
        fragdict.update({pepform: {}})            # Initiate main key (peptidoform).
        pep = proforma.ProForma.parse(modpep)   # Initiate proforma object.

        # Add all a/b/y ions singly and doubly charged.
        for ion_type in ['a','b','y']:
            for charge in [1,2]:
                fragments = pep.fragments(ion_type,charge=charge) 
                for series, mz in enumerate(fragments):
                    ion = f'{ion_type}{series+1}{charge*"+"}'
                    fragdict[pepform].update({ion: mz})

                    # Add water and ammonia neutral losses.
                    for loss in loss_masses:
                        loss_ion = f'{ion}{loss}'
                        fragdict[pepform].update({loss_ion: mz - loss_masses[loss]})

    # Determine which MS2PIP model to use.
    if 'TMT' in args.mod: model = 'TMT'
    elif instrument == 'timsTOF': model = 'timsTOF'
    elif instrument == 'orbitrap': model = 'Immuno-HCD'
    else:
        logging.warning('Undefined MS2PIP model - please check instrument/modification.')
        return
    
    # Run the MS2PIP prediction in batch mode.
    to_predict = pd.DataFrame({
        'peptidoform': uniq_pepforms,
        'spectrum_id': range(0,len(uniq_pepforms))
        })
    logging.info(f'Starting MS2PIP predictions using the {model} model.')
    to_predict.to_csv(f"{plot_folder}/to_predict.tsv",sep='\t',index=False)
    os.system(f"ms2pip predict-batch {plot_folder}/to_predict.tsv --model {model} > {plot_folder}/ms2pip.log 2>&1")
    
    # Read predictions.
    predictions = pd.read_csv(f'{plot_folder}/to_predict_predictions.tsv',sep='\t')
    predictions['intensity'] = predictions['predicted'].apply(lambda x: float((2 ** x) - 0.001))   # Undo log2 transformation.
    predictions['ion'] = predictions['ion_type'] + predictions['ion_number'].astype(str) + '+' 
    predictions = predictions.merge(to_predict,left_on='psm_index',right_on='spectrum_id',how='left')
    predictions['model'] = model
    predictions = predictions.drop(columns=['observed','predicted','rt','im','spectrum_id','ion_type','psm_index','ion_number'])
    predictions.to_csv(f'{plot_folder}/ms2pip_predictions.tsv',sep='\t',index=False)

    # Iterate over scans sample per sample - using the created MGFs.
    logging.info('Collected MS2PIP predictions - now plotting experimental/predicted mirror plots - can take a while!')
    for run in df['run'].unique():
        logging.info(f'Now plotting spectra for {run}.')
        df_run = df[df['run'] == run].copy()
        df_run['scanid'] = df_run['run'] + '_' + df_run['scan']
        df_run.set_index('scanid',inplace=True)
        mgf_file = f'{args.in_dir}/report/scans/{run}.mgf'
        with mgf.read(mgf_file) as reader:
            for spectrum in reader:
                title = spectrum.get('params', {}).get('title', 'N/A')
                match = re.search(r'scan=(\d+)_z=(\d+)', title)
                if match:
                    scan = match.group(1)
                    scanid = f'{run}_{scan}'
                    charge = match.group(2)
                    if scanid in df_run.index:
                        charge = match.group(2)
                        pepmass = spectrum.get('params', {}).get('pepmass', 'N/A')
                        mz_array, intensity_array = remove_precursor(spectrum['m/z array'],spectrum['intensity array'],pepmass[0])
                        intensity_array_scaled = transform_int(intensity_array)
                        
                        # It's possible multiple peptidoforms are found per spectrum by different search engines, plot all.
                        row = df_run.loc[f'{run}_{scan}']
                        peptidoforms, rt, spec_pearsons, spec_rt_diff_best = row[['peptidoform', 'rt', 'rescoring:spec_pearson', 'rescoring:rt_diff_best']]
                        pepform2pearson = dict(zip(peptidoforms,spec_pearsons))
                        pepform2rt = dict(zip(peptidoforms,spec_rt_diff_best))
                        for pepform in list(set(peptidoforms)):
                            if pepform in fragdict:
                                
                                # Initiate with the experimental spectrum.
                                exp_df = pd.DataFrame({'mz': mz_array,'int': intensity_array_scaled})
                                exp_df = exp_df.astype({'mz': 'float', 'int': 'float'})

                                # MS2PIP prediction.
                                pred_df = predictions[predictions['peptidoform'] == pepform].copy()
                                scaled_intensities = transform_int(pred_df['intensity'])
                                pred_df.loc[:, 'int'] = [x * -1 for x in scaled_intensities]    # Negative for mirror plot.
                                
                                # Annotate ions and their colors.
                                exp_df['ion'] = exp_df['mz'].apply(lambda x: annotate_ion(x, fragdict[pepform]))
                                
                                # Concatenate these
                                exp_df = pd.concat([exp_df, pred_df[['mz','int','ion']]], ignore_index=True) 

                                # Annotate ions.
                                exp_df['color'] = exp_df['ion'].apply(lambda x: ion_color(x))
                                exp_df['sort_order'] = exp_df['ion'].apply(lambda x: 0 if x == 'NA' else 1)
                                exp_df = exp_df.sort_values(by=['sort_order', 'mz'])

                                # Require at least 2% intensity
                                exp_df = exp_df[abs(exp_df['int']) > 2]
                                
                                # Make figure.
                                plt.figure(figsize=(10,6))
                                plt.vlines(exp_df['mz'], 0, exp_df['int'],colors=exp_df['color'], linewidth=1)
                                plt.xlabel('m/z')
                                plt.ylabel('Relative Intensity')
                                plt.axhline(y=0, linewidth=1, color='black')
                                plt.ylim(-110,110)
                                ticks = range(-100,125,25)  # Tick positions
                                labels = ['100%', '75%', '50%', '25%', '0%', '25%', '50%', '75%', '100%']
                                plt.yticks(ticks=ticks, labels=labels)
                                plt.grid(False)
                                plt.suptitle(f'{pepform}', fontsize=14)
                                plt.title(f'scan {scan} - {pepmass[0]} m/z - {float(rt):.3f} min\n'
                                        f'spec_pearson: {pepform2pearson[pepform]:.2f} - rt_diff_best: {pepform2rt[pepform]:.2f} min',fontsize=10, loc='center', pad=10)

                                # Show ion labels if higher than 5% intensity.
                                for index, row in exp_df.iterrows():
                                    if abs(row['int']) > 5 and row['ion'] != 'NA':
                                        plt.text(row['mz'], row['int']+0.05 if row['int'] > 0 else row['int']-0.05, format_ion_text(row['ion']),
                                                color=row['color'], ha='center', va='center')

                                # Save as PNG and SVG.
                                plain_pep = re.sub(r'\[.*?\]', '', pepform).translate(str.maketrans('', '', '0123456789/-'))
                                plt.savefig(f'{plot_folder}/{plain_pep}_{run}_{scan}.png',dpi=600)
                                plt.savefig(f'{plot_folder}/{plain_pep}_{run}_{scan}.svg',dpi=600)
                                plt.clf()
                                plt.close()
                            else:
                                logging.warning(f'No fragments were annotated for {pepform}. Skipping this one.')
                else:
                    logging.warning(f'MGF spectrum reading error for spectrum title {title}')

# Antigen histogram plotter.
def antigen_plotter(args, df, name, plot_folder):
    
    # Check if HLA binding info present.
    if 'bind_level' in df.columns:
        df.loc[~df['immuno'],'bind_level'] = 'NA'       # Only consider HLA info for specified immunolengths.
        agg_df = pd.crosstab(index=df['gene'], columns=df['bind_level'])
        agg_df = agg_df.reindex(columns=['SB','WB','NB','NA'])
        colors = ['firebrick', 'dodgerblue', 'grey', 'gainsboro']
    else:
        agg_df = pd.crosstab(index=df['gene'], columns=df['immuno'])
        colors = ['firebrick', 'dodgerblue']

    # Sort from highest to lowest number of peptides
    agg_df['total'] = agg_df.sum(axis=1)
    total_peptides = agg_df['total'].sum()
    agg_df = agg_df.sort_values(by=['total','SB','WB','NB'], ascending=[False,False,False,False]).drop(columns='total')
    agg_df.plot(kind='bar', stacked=True, figsize=(10, 6), color=colors, edgecolor='black')

    plt.title(f'Antigen histogram plot - {name}')
    plt.xlabel(f'Proteins (N={len(agg_df)})')
    plt.ylabel(f'Number of peptides (N={total_peptides})')
    plt.tight_layout()
    
    # Save as PNG and SVG.
    plt.savefig(f'{plot_folder}/{name}_epitopes.png',dpi=600)
    plt.savefig(f'{plot_folder}/{name}_epitopes.svg',dpi=600)
    plt.close()

# Plot number of epitopes per antigen.
def antigen_plot(args, peptides) -> None:
    # Filter for xeno-epitopes, prelimary check.
    df = peptides[~peptides['species'].str.contains(f'contaminant|{args.host.upper()}', regex=True)]
    if len(df) == 0:
        logging.warning('No non-host/contaminant proteins, skipping antigen histogram plot.')
        return
    
    logging.info(f'Found {len(df)} non-host/contaminant peptides pre-filtering, plotting epitope histogram for antigens.')
    plot_folder = f'{args.in_dir}/report/plots/antigen_histogram/'
    make_folder(plot_folder)

    # Plot first for all antigen peptides
    antigen_plotter(args, df, 'all', plot_folder)

    # Check for uninfected condition PSM count column.
    uninf_column = [x for x in df.columns if x.startswith('PSMs_cond_uninf')]
    if len(uninf_column) == 1:
        df = df[ (df[ uninf_column[0] ] == 0) & (df['#chimeric_PSMs'] == 0)]     # Reduce to peptides with no hits in uninfected condition and chimera.
        antigen_plotter(args, df, 'no_uninfected', plot_folder)

# Rescore plot per search engine.
def rescore_plot(args) -> None:
    logging.info('Plotting target-decoy mokapot score distributions per search engine before and after rescoring.')
    indir = args.in_dir
    plot_folder = f'{indir}/report/plots/rescore/'
    make_folder(plot_folder)

    # Append all PSMs to this dataframe
    columns = ['mokapot score','mokapot q-value','features','engine','tda']
    all_psms = pd.DataFrame(columns=columns)
    hits = {}

    for psm_file in glob.glob(f'{indir}/report/mokapot/*mokapot*psms.txt'):
        engine = os.path.basename(psm_file).split('_')[0]
        psm = pd.read_csv(psm_file,sep='\t')
        psm['engine'] = engine
        feature_set = 'ms2rescore' if '_ms2rescore.' in psm_file else 'search'
        psm['features'] = feature_set
        if 'decoy.psm' in psm_file:
            psm['tda'] = 'd'
        else:
            psm['tda'] = 't'
            below_1q = (psm['mokapot q-value'] < 0.01).sum()
            hits.update({f'{feature_set} | {engine}': str(below_1q)})
        all_psms = pd.concat([all_psms, psm[columns]], ignore_index=True)
        
    # Set color palette for the 'tda' column
    palette = {'d': 'red', 't': 'blue'}

    # Create a FacetGrid for separate ridge plots per search engine
    g = sns.FacetGrid(all_psms, row='features', col='engine', hue='tda', palette=palette, sharey=True, height=5,aspect=0.8)
    g.map(sns.kdeplot, 'mokapot score', fill=True, alpha=0.6, common_norm=False)

    # Adjust the appearance of each plot
    g.set_axis_labels('mokapot score', 'Density')
    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.add_legend(title='Target-Decoy')

    # Add custom text to each subplot based on the 'hits' column
    for ax in g.axes.flat:
        ax.set_xlim((-5,5))
        # Add the identified PSMs per search engine/set.
        ax.text(0.05, 0.95, f'{hits[ax.get_title()]} PSMs\n(q < 0.01)', transform=ax.transAxes, fontsize=12,
        verticalalignment='top', horizontalalignment='left',
        bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    # Adjust the layout and display the plot
    plt.subplots_adjust(top=0.9, hspace=0.3, wspace=0.3)
    plt.tight_layout()

    # Save as PNG and SVG.
    plt.savefig(f'{plot_folder}/rescore_ridgeplot.png',dpi=600)
    plt.savefig(f'{plot_folder}/rescore_ridgeplot.svg',dpi=600)
    plt.close()

# Barplot of identified peptide sequences.
def peptide_barplot(args, peptides) -> None:
    logging.info('Plotting peptide identification barplot for rescoring and engine aggregation effect.')

    # Peptides here are after re-scoring, just filter for the needed columns
    df = peptides[['peptide','engine']].copy()
    df['set'] = 'ms2rescore'

    # Now get the one prior re-scoring
    columns = ['peptide','engine']
    df_search = pd.DataFrame(columns=columns)
    peptide_files = glob.glob(f'{args.in_dir}/report/mokapot/*search*mokapot.peptides.txt')
    for mokapot_peptide in peptide_files:
        engine = os.path.basename(mokapot_peptide).split('_')[0]
        pepfile = pd.read_csv(mokapot_peptide,sep='\t')
        pepfile = pepfile[pepfile['mokapot q-value'] < 0.01]
        pepfile['engine'] = engine
        if engine == 'MSFragger':
            pepfile['Peptide'] = pepfile['Peptide'].apply(lambda x: x[2:-3])                  # e.g. R.SHYEEGPGKNLPFSVENKWS3.L
        elif engine == 'Comet':
            pepfile['Peptide'] = pepfile['Peptide'].apply(lambda x: x[2:-2])                  # e.g. R.SHYEEGPGKNLPFSVENKWS.L

        pepfile['peptide'] = pepfile['Peptide'].apply( lambda x: re.sub(r'[^A-Z]', '', re.sub(r'\[.*?\]', '', x)) )
        df_search = pd.concat([df_search, pepfile[['peptide','engine']]], ignore_index=True)

    # Aggregate per search engine.
    df_search = df_search.groupby('peptide')['engine'].agg(lambda x: list(set(x))).reset_index()
    df_search['set'] = 'search'

    # Function to count peptides per engine
    def count_peptides_per_engine(df):
        engine_list = [engine for engines in df['engine'] for engine in engines]
        engine_list_sorted = sorted(engine_list)
        return dict(Counter(engine_list_sorted))

    # Count peptides for each set and export the numbers.
    count_ms2rescore = count_peptides_per_engine(df)
    count_ms2rescore.update({'Combined': df['peptide'].nunique()})
    count_search = count_peptides_per_engine(df_search)
    count_search.update({'Combined': df_search['peptide'].nunique()})

    # Initiate figure plot and define colors per search engine.
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    engine_colors = {'Comet': 'cadetblue','Sage': 'green', 'PEAKS': 'darkorange', 'MSFragger': 'darkmagenta','Combined': 'firebrick'}
    
    # Left-hand side: peptide identifications after search (prior ms2rescore).
    search_plot = axes[0].bar(count_search.keys(), count_search.values(),
            edgecolor = 'black', color=[engine_colors.get(e, 'lightgray') for e in count_search])
    axes[0].set_title('search')
    axes[0].set_xlabel('Search engine')
    axes[0].set_ylabel('Identified peptides')
    for bar in search_plot:     # Add the numbers above the bar.
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width() / 2, height,
                 f'{height}',
                 ha='center', va='bottom', fontsize=10, color='black')

    # Right-hand side: peptide identifications after ms2rescore.
    ms2rescore_plot = axes[1].bar(count_ms2rescore.keys(), count_ms2rescore.values(),
            edgecolor = 'black',color=[engine_colors.get(e, 'lightgray') for e in count_ms2rescore])
    axes[1].set_title('ms2rescore')
    axes[1].set_xlabel('Search engine')
    for bar in ms2rescore_plot:     # Add the numbers above the bar.
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width() / 2, height,
                 f'{height}',
                 ha='center', va='bottom', fontsize=10, color='black')

    # Save as PNG and SVG.
    plt.tight_layout()
    plot_folder = f'{args.in_dir}/report/plots/rescore/'
    plt.savefig(f'{plot_folder}/peptide_barplot.png',dpi=600)
    plt.savefig(f'{plot_folder}/peptide_barplot.svg',dpi=600)
    plt.close()

def polygon_plot(args, scans) -> None:
    logging.info('Plotting polygons per run.')
    plot_folder = f'{args.in_dir}/report/plots/polygon/'
    make_folder(plot_folder)
    
    # Polygon plots per MS run separately.
    for run in scans['run'].unique():
        df_run = scans[scans['run'] == run].copy()             # Reduce per run.
        df_run[['mz','im']] = df_run[['mz','im']].astype(float) # Convert to float for plotting.

        # Split here non-identified, identified non-immunopeptides and identified immunopeptides.
        df_run['id_type'] = 'no_id' # Initialize.
        df_run.loc[df_run['peptide'].notna(), 'id_type'] = 'id'
        df_run.loc[df_run['immuno'] == True, 'id_type'] = 'immuno'
        non_id = df_run[df_run['id_type'] == 'no_id']
        id_noimmuno = df_run[df_run['id_type'] == 'id']
        id_immuno = df_run[df_run['id_type'] == 'immuno']

        # Initiate figure.
        figure, axis = plt.subplots(2, 2, figsize=(15, 10),width_ratios=[1, 3],height_ratios=[3, 1])
        max_im = df_run['im'].max()*1.05
        min_im = df_run['im'].min()*0.95
        max_mz = df_run['mz'].max()*1.05
        min_mz = df_run['mz'].min()*0.95

        figure.suptitle('Polygon plot for '+run,size="large")
        sns.set(style="ticks")
        sns.scatterplot(x="mz", y="im",data=non_id,color="grey",s=4,alpha=0.2,ax=axis[0,1])
        sns.scatterplot(x="mz", y="im",data=id_noimmuno,color="blue",s=4,ax=axis[0,1])
        sns.scatterplot(x="mz", y="im",data=id_immuno,color="red",s=4,ax=axis[0,1])
        axis[0,1].text(max_mz*0.2, max_im*0.96, f'{len(id_immuno)} PSMs {args.len} mers',color="red")
        axis[0,1].text(max_mz*0.2, max_im*0.93, f'{len(id_noimmuno)} PSMs other length',color="blue")
        axis[0,1].text(max_mz*0.2, max_im*0.90, f'{len(non_id)} non-ID spectra',color="grey")
        axis[0,1].set(xlabel=None)
        axis[0,1].set(ylabel=None)
        axis[0,1].set(xticklabels=[])
        axis[0,1].set(yticklabels=[])
        axis[1,1].set_xlabel("Observed m/z")
        axis[1,1].set_ylabel("Identified spectra")
        sns.histplot(data=df_run,x="mz", hue="id_type",multiple="layer",hue_order=['no_id','immuno','id'],palette=["grey", "red", "blue"],ax=axis[1,1])
        axis[0,0].set_xlabel("Identified spectra")
        axis[0,0].set_ylabel("Ion mobility 1/K0")
        sns.histplot(data=df_run,y="im", hue="id_type",multiple="layer",hue_order=['no_id','immuno','id'],palette=["grey", "red", "blue"],ax=axis[0,0])
        axis[1,0].axis('off')
        figure.tight_layout()
        figure.savefig(f'{plot_folder}/{run}.svg',dpi=600)
        figure.savefig(f'{plot_folder}/{run}.jpg',dpi=600)
        figure.clf()
        plt.close()

# Plot Pearson correlation (MS2PIP), Retention time (DeepLC) and ion mobility feature (IM2Deep) performance for PSMs < 1% qval.
def rescoring_features(args, peptides, instrument):
    logging.info('Plotting feature performance for all search engines.')
    indir = args.in_dir
    plot_folder = f'{indir}/report/plots/rescore/'

    # Parse from ms2rescore psm.tsv results.
    psm_files = glob.glob(f'{indir}/report/ms2rescore/*.psms.tsv')
    num_engines = len(psm_files)
    fig, axes = plt.subplots(2, num_engines, figsize=(5 * num_engines, 8))
    
    # If there's only one search engine, axes will be 1-dimensional
    if num_engines == 1:
        axes = axes[:, None]

    # Determine features to plot, if timsTOF we have CCS predictions as well.
    columns=['engine','scan_id','rescoring:spec_pearson','rescoring:observed_retention_time','rescoring:predicted_retention_time','rescoring:rt_diff_best']
    if instrument == 'timsTOF':
        fig, axes = plt.subplots(3, num_engines, figsize=(5 * num_engines, 8))
        columns = columns + ['rescoring:ccs_predicted_im2deep','rescoring:ccs_observed_im2deep']
    else:
        fig, axes = plt.subplots(2, num_engines, figsize=(5 * num_engines, 8))
    features = pd.DataFrame(columns=columns)
    
    # Plot the feature performance for all PSMs per search engine.
    for idx, psm_file in enumerate(psm_files):
        engine = os.path.basename(psm_file).split('.')[0]   # E.g. MSFragger.psms.tsv
        data = pd.read_csv(psm_file, sep='\t')              # Read psms.
        data = data[(data['meta:peptide_qvalue'] <= 0.01) & (data['qvalue'] <= 0.01)]
        data['scan_id'] = data['run'] + '_' + data['spectrum_id'].astype(str)
        data['engine'] = engine
        columns = ['scan_id', 'rescoring:spec_pearson', 'rescoring:observed_retention_time', 'rescoring:predicted_retention_time','engine']
        if instrument == 'timsTOF':
            columns = columns + ['rescoring:ccs_predicted_im2deep','rescoring:ccs_observed_im2deep']
        df = data[columns]
        features = pd.concat([features,df],ignore_index=True)
        
        # Histogram of Pearson correlation
        sns.set(style='ticks')
        sns.histplot(df['rescoring:spec_pearson'], ax=axes[0, idx], kde=False, bins=100, color='firebrick')
        median_pearson = df['rescoring:spec_pearson'].median()
        axes[0, idx].axvline(median_pearson, color='black', linestyle='--', label=f'Median: {median_pearson:.2f}')
        axes[0, idx].set_title(f'{engine} - MS2PIP')
        axes[0, idx].set_xlabel('Pearson Correlation')
        if idx == 0:
            axes[0, idx].set_ylabel('PSMs')
        else:
            axes[0, idx].set_ylabel('')
        axes[0, idx].legend(loc='upper center')

        # Scatterplot of observed vs predicted retention time
        sns.scatterplot(x='rescoring:observed_retention_time', y='rescoring:predicted_retention_time', data=df, ax=axes[1, idx],
                color='firebrick', alpha=0.05,edgecolor='firebrick')
        axes[1, idx].set_title(f'{engine} - DeepLC')
        axes[1, idx].set_xlabel('Observed RT')
        if idx == 0:
            axes[1, idx].set_ylabel('Predicted RT')
        else:
            axes[1, idx].set_ylabel('')
        # Add diagonal line on top
        min_retention = df[['rescoring:observed_retention_time', 'rescoring:predicted_retention_time']].min().min()
        max_retention = df[['rescoring:observed_retention_time', 'rescoring:predicted_retention_time']].max().max()
        axes[1, idx].plot([min_retention, max_retention], [min_retention + 1, max_retention + 1], color='black', linestyle='--')

        if instrument == 'timsTOF':
            # Scatterplot of observed vs predicted retention time
            sns.scatterplot(x='rescoring:ccs_observed_im2deep', y='rescoring:ccs_predicted_im2deep', data=df, ax=axes[2, idx],
                color='firebrick', alpha=0.05,edgecolor='firebrick')
            axes[2, idx].set_title(f'{engine} - IM2Deep')
            axes[2, idx].set_xlabel('Observed CCS')
            if idx == 0:
                axes[2, idx].set_ylabel('Predicted CCS')
            else:
                axes[2, idx].set_ylabel('')
            # Add diagonal line on top
            min_retention = df[['rescoring:ccs_observed_im2deep', 'rescoring:ccs_predicted_im2deep']].min().min()
            max_retention = df[['rescoring:ccs_observed_im2deep', 'rescoring:ccs_predicted_im2deep']].max().max()
            axes[2, idx].plot([min_retention, max_retention], [min_retention + 1, max_retention + 1], color='black', linestyle='--')

    # Adjust layout
    plt.tight_layout()
    plt.savefig(f'{plot_folder}/rescoring_features.svg',dpi=600)
    plt.savefig(f'{plot_folder}/rescoring_features.jpg',dpi=600)
    plt.clf()
    plt.close()

    # Now plot the feature performance per peptide.
    logging.info('Plotting best Pearson correlation for peptides in final reports.')
    df = peptides[['peptide','species','immuno','rescoring:spec_pearson','rescoring:rt_diff_best']].copy()
    df['bind_level'] = peptides.get('bind_level', 'NA')     # Add binding information if present, otherwise 'NA'.
    df.loc[~df['immuno'], 'bind_level'] = 'NA'              # Only consider binding predictions for peptides specified as immunopeptides.
    bind_level_order = ['SB', 'WB', 'NB', 'NA']             # Order from SB -> WB ->NB (no binder) - >NA (not immuno)
    df['bind_level'] = pd.Categorical(df['bind_level'], categories=bind_level_order, ordered=True)
    
    # Plot and add text with number of peptides
    palette = {'SB': 'firebrick', 'WB': 'dodgerblue', 'NB': 'grey', 'NA': 'gainsboro'}
    sns.set(style='ticks')
    plt.figure(figsize=(12, 8))
    ax = sns.boxplot(x='species', y='rescoring:spec_pearson', hue='bind_level', data=df, showfliers=False, palette=palette)
    
    # Calculate the number of data points for each group and add text annotations
    for i, species in enumerate(df['species'].unique()):
        for j, bind_level in enumerate(bind_level_order):
            peptides_plotted = df[(df['species'] == species) & (df['bind_level'] == bind_level)]
        
            # Calculate the position for annotation
            x_pos = i + 0.1  + (j - len(bind_level_order) / 2) * 0.2  # Adjust x position based on the order and spacing
            ax.text(x_pos, 1.02, len(peptides_plotted), ha='center', va='bottom', color='black', fontsize=10)

    plt.title('Spectrum Pearson correlation per species.')

    # Layout and save.
    plt.tight_layout()
    plt.savefig(f'{plot_folder}/peptides_pearson.svg',dpi=600)
    plt.savefig(f'{plot_folder}/peptides_pearson.jpg',dpi=600)
    plt.clf()
    plt.close()
    
