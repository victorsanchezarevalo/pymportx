import pandas as pd
import numpy as np
import os
import warnings
import json
from .auxiliary_functions import *

def read_kallisto(folders,
                tx_out=True,
                tx2gene=None,
                countsFromAbundance='no',
                dropInfReps=False,
                varReduce=False,
                infRepStat=None,
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

    # Check if infRepStat method is valid
    valid_infRepStat = [None, "median"]
    if infRepStat not in valid_infRepStat:
        raise ValueError(f"Invalid infRepStat: {infRepStat}")
    
    # Read abundance.h5 or abundance.tsv.gz files 
    h5_file_paths = []
    tsv_file_paths = []

    for folder in folders:
        h5_file_path = os.path.join(folder, 'abundance.h5')
        tsv_file_path = os.path.join(folder, 'abundance.tsv.gz')
        run_info_path = os.path.join(folder, 'run_info.json')

        if os.path.exists(h5_file_path): 
            h5_file_paths.append(h5_file_path)
            print(f"Found h5 file in {folder}, using this file as it is faster.")
        elif os.path.exists(tsv_file_path):
            tsv_file_paths.append(tsv_file_path)        
        else:
            raise ValueError(f"Folder {folder} does not contain 'abundance.h5' or 'abundance.tsv.gz' file")   

        if os.path.exists(run_info_path):
            with open(run_info_path, 'r') as f:
                run_info = json.load(f)
                num_bootstrap = run_info.get('n_bootstraps', 0)
        else:
            raise ValueError(f"Folder {folder} does not contain 'run_info.json' file")
        
        # Check if inferential replicate information exists
        if dropInfReps==False and num_bootstrap == 0:
            warnings.warn(f"Folder {folder} does not contain inferential replicate infromation. 'dropInfReps' will be set to 'True'")
            dropInfReps = True

    # Check if dropInfReps is True, then infRepStat and varReduce should be None
    if dropInfReps:
        warnings.warn("Setting 'dropInfReps=True' will ignore 'infRepStat' and 'varReduce'")
        varReduce = False
        infRepStat = None
    
    # Process infRepStat
    infRepStat = process_infRepStat(infRepStat)

    # Create sample names
    sample_names = [os.path.basename(folder) for folder in folders]

    # Create variables to store the data
    abundanceMatTx_list = []
    countsMatTx_list = []
    lengthMatTx_list = []
    varMatTx_list = []
    infRepMatTx_dict = {}

    # Read sample data
    for i in range(len(folders)):
        print(f"Reading {sample_names[i]}")
        if h5_file_paths:
            sample = read_kallisto_h5(h5_file_paths[i])
        else:
            sample = pd.read_csv(tsv_file_paths[i], 
                                sep='\t', 
                                dtype={'col1': str, 'col2': int, 'col3': float, 'col4': float, 'col5': float})

        if i==0:
            tx_names = sample['target_id'].to_list()  # save the transcript names of the first sample
        else:
            if tx_names != sample['target_id'].to_list(): # check if transcript names are the same across samples
                raise ValueError("Transcript names are not the same across samples")
        
        lengthMatTx_list.append(sample['eff_length'].reset_index(drop=True)) 
        
        # Read replicates
        if dropInfReps:
            countsMatTx_list.append(sample['est_counts'].reset_index(drop=True))
            abundanceMatTx_list.append(sample['tpm'].reset_index(drop=True))
        else:
            # Read inferential replicate data
            boots_info = readInfRepKallisto(h5_file_paths[i])

            if infRepStat is not None:
                countsMatTx_list.append(pd.Series(infRepStat(boots_info['boots'])))
                tpm = np.array(countsMatTx_list[i]) / np.array(lengthMatTx_list[i])
                abundanceMatTx_list.append(pd.Series(tpm*1e6/np.sum(tpm)))
            else:
                countsMatTx_list.append(sample['est_counts'].reset_index(drop=True))
                abundanceMatTx_list.append(sample['tpm'].reset_index(drop=True))

            if varReduce and tx_out: 
                varMatTx_list.append(pd.Series(boots_info['boots_vars']))
            else:
                infRepMatTx_dict[sample_names[i]] = boots_info['boots']        

    abundanceMatTx = pd.concat(abundanceMatTx_list, axis=1)
    abundanceMatTx.columns = sample_names
    abundanceMatTx.index = tx_names

    countsMatTx = pd.concat(countsMatTx_list, axis=1)
    countsMatTx.columns = sample_names
    countsMatTx.index = tx_names

    lengthMatTx = pd.concat(lengthMatTx_list, axis=1)
    lengthMatTx.columns = sample_names
    lengthMatTx.index = tx_names

    txi = {
        'abundance': abundanceMatTx,
        'counts': countsMatTx,
        'length': lengthMatTx,
        'countsFromAbundance': countsFromAbundance}

    if not dropInfReps:
        if varReduce and tx_out:
            varMatTx = pd.concat(varMatTx_list, axis=1)
            varMatTx.columns = sample_names
            varMatTx.index = tx_names
            txi['variance'] = varMatTx
        elif varReduce and not tx_out:
            txi['infReps'] = infRepMatTx_dict
        if not varReduce:
            txi['infReps'] = infRepMatTx_dict

    
    if tx_out:
        if countsFromAbundance != "no":
            length4CFA = txi['length']
            if countsFromAbundance == "dtuScaledTPM":
                length4CFA = medianLengthOverIsoform(length4CFA, tx2gene, ignoreTxVersion, ignoreAfterBar)
                countsFromAbundance = "lengthScaledTPM"
            
            txi['counts'] = makeCountsFromAbundance(txi['counts'], txi['abundance'], length4CFA, countsFromAbundance)
        return txi
    else:
        txiGene = summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)       
        return txiGene