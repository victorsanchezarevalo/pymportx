import pandas as pd
import numpy as np
import os
import glob
import warnings
from .auxiliary_functions import *

def read_rsem(folders,
                tx_in =True,
                tx_out=True,
                tx2gene=None,
                countsFromAbundance='no',
                ignoreTxVersion=False,
                ignoreAfterBar=False):

    # Check if folders is a non-empty list of strings
    if not isinstance(folders, list) or not folders:
        raise ValueError("Folders must be a non-empty list of strings")
    
    # Check if tx2gene file exists
    if tx2gene is not None:
        if not os.path.exists(tx2gene):
            raise ValueError("'tx2gene' file does not exist")
        
    # Check the method to generate estimated counts using abundance estimates
    valid_countsFromAbundance = ["no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"]
    if countsFromAbundance not in valid_countsFromAbundance:
        raise ValueError(f"Invalid 'countsFromAbundance': {countsFromAbundance}")
    else:
        if countsFromAbundance != "no" and tx2gene is None:
            raise ValueError("tx2gene file must be provided to use 'countsFromAbundance'")
        if countsFromAbundance == "dtuScaledTPM" and not tx_out:
            raise ValueError("dtuScaledTPM can only be used with tx_out=True")

    if not tx_in and tx_out:
        raise ValueError("tx_out is only an option when transcript-level data is read in tx_in=True")
    
    if not tx_in and countsFromAbundance != "no":
        warnings.warn("countsFromAbundance is only used when tx_in=True (transcript level)")

    # Read isoforms.results.gz and genes.results.gz files 
    transcript_file_paths = []
    gene_file_paths = []

    for folder in folders:
        transcript_file_path = os.path.join(folder, '*isoforms.results.gz')
        transcript_files = glob.glob(transcript_file_path)
        gene_file_path = os.path.join(folder, '*genes.results.gz')
        gene_files = glob.glob(gene_file_path)

        if not transcript_files and not gene_files:
            raise ValueError(f"Folder {folder} does not contain 'isoforms.results.gz' file or 'genes.results.gz' file")
        else:
            transcript_file_paths.extend(transcript_files)
            gene_file_paths.extend(gene_files)
        

    # Create sample names
    sample_names = [os.path.basename(folder) for folder in folders]

    # Create variables to store the data
    abundanceMatTx_list = []
    countsMatTx_list = []
    lengthMatTx_list = []

    # Determine file paths and id_column based on the value of tx_in
    if tx_in:
        file_paths = transcript_file_paths
        id_column = 'transcript_id'
    else:
        file_paths = gene_file_paths
        id_column = 'gene_id'

    # Read sample data
    for i, file_path in enumerate(file_paths):
        print(f"Reading {sample_names[i]}")
        sample = pd.read_csv(file_path, 
                            sep='\t', 
                            dtype={'col1': str, 'col2': str, 'col3': float, 'col4': float, 'col5': float, 'col6': float, 'col7': float, 'col8': float})

        if i==0:
            names = sample[id_column].to_list()  # save the transcript names of the first sample
        else:
            if names != sample[id_column].to_list(): # check if transcript names are the same across samples
                raise ValueError("Transcript names are not the same across samples")
        
        lengthMatTx_list.append(sample['effective_length'].reset_index(drop=True)) 
        countsMatTx_list.append(sample['expected_count'].reset_index(drop=True))
        abundanceMatTx_list.append(sample['TPM'].reset_index(drop=True))
        
    abundanceMatTx = pd.concat(abundanceMatTx_list, axis=1)
    abundanceMatTx.columns = sample_names
    abundanceMatTx.index = names

    countsMatTx = pd.concat(countsMatTx_list, axis=1)
    countsMatTx.columns = sample_names
    countsMatTx.index = names

    lengthMatTx = pd.concat(lengthMatTx_list, axis=1)
    lengthMatTx.columns = sample_names
    lengthMatTx.index = names

    txi = {
        'abundance': abundanceMatTx,
        'counts': countsMatTx,
        'length': lengthMatTx,
        'countsFromAbundance': countsFromAbundance}

        
    if tx_out:
        if countsFromAbundance != "no":
            length4CFA = txi['length']
            if countsFromAbundance == "dtuScaledTPM":
                length4CFA = medianLengthOverIsoform(length4CFA, tx2gene, ignoreTxVersion, ignoreAfterBar)
                countsFromAbundance = "lengthScaledTPM"
            
            txi['counts'] = makeCountsFromAbundance(txi['counts'], txi['abundance'], length4CFA, countsFromAbundance)
        return txi
    else:
        if tx_in: 
            txiGene = summarizeToGene(txi, tx2gene, False, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)
            return txiGene
        else:
            return txi