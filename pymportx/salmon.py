import pandas as pd
import numpy as np
import os
import warnings
from .auxiliary_functions import *
import anndata as ad

def read_salmon(folders,
                tx_out=True,
                tx2gene=None,
                countsFromAbundance='no',
                dropInfReps=False,
                varReduce=False,
                infRepStat=None,
                ignoreTxVersion=False,
                ignoreAfterBar=False):

    if not isinstance(folders, list) or not folders:
        raise ValueError("Folders must be a non-empty list of strings")
    
    if tx2gene is not None and not os.path.exists(tx2gene):
        raise ValueError("'tx2gene' file does not exist")
        
    valid_countsFromAbundance = ["no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"]
    if countsFromAbundance not in valid_countsFromAbundance:
        raise ValueError(f"Invalid 'countsFromAbundance': {countsFromAbundance}")
    elif countsFromAbundance == "dtuScaledTPM":
        if tx2gene is None:
            raise ValueError("tx2gene file must be provided to use 'dtuScaledTPM'")
        if not tx_out:
            raise ValueError("dtuScaledTPM can only be used with tx_out=True")

    # Check if infRepStat method is valid
    valid_infRepStat = [None, "median"]
    if infRepStat not in valid_infRepStat:
        raise ValueError(f"Invalid infRepStat: {infRepStat}")
    
    # Read quant.sf.gz files and meta_info.json files
    quant_file_paths = []
    meta_info_paths = []
    bootstraps_paths = []
    for folder in folders:
        quant_file_path = os.path.join(folder, 'quant.sf.gz')
        if not os.path.exists(quant_file_path):
            raise ValueError(f"Folder {folder} does not contain 'quant.sf.gz' file")       
        else:
            quant_file_paths.append(quant_file_path)

        meta_info_path = os.path.join(folder, 'aux_info/meta_info.json')
        if not os.path.exists(meta_info_path):
            raise ValueError(f"Folder {folder} does not contain 'meta_info.json' file")
        else:
            meta_info_paths.append(meta_info_path)
        
        # Check if bootstrap directory exists
        if dropInfReps==False:
            bootstraps_path = os.path.join(folder, 'aux_info/bootstrap/bootstraps.gz')
            if os.path.exists(bootstraps_path):
                bootstraps_paths.append(bootstraps_path)
            else:
                warnings.warn(f"Folder {folder} does not contain 'bootstraps.gz' file. 'dropInfReps' will be set to 'True'")
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
        sample = pd.read_csv(quant_file_paths[i], 
                             sep='\t', 
                             dtype={'col1': str, 'col2': int, 'col3': float, 'col4': float, 'col5': float})

        if i==0:
            tx_names = sample['Name'].to_list()  # save the transcript names of the first sample
        else:
            if tx_names != sample['Name'].to_list(): # check if transcript names are the same across samples
                raise ValueError("Transcript names are not the same across samples")
        
        lengthMatTx_list.append(sample['EffectiveLength'].reset_index(drop=True))
        
        # Read replicates
        if dropInfReps:
            countsMatTx_list.append(sample['NumReads'].reset_index(drop=True))
            abundanceMatTx_list.append(sample['TPM'].reset_index(drop=True))
        else:
            # Read bootstraps file
            boots_info = readInfRepFish('salmon', meta_info_paths[i], bootstraps_paths[i])

            if infRepStat is not None:
                countsMatTx_list.append(pd.Series(infRepStat(boots_info['boots'])))
                tpm = np.array(countsMatTx_list[i]) / np.array(lengthMatTx_list[i])
                abundanceMatTx_list.append(pd.Series(tpm*1e6/np.sum(tpm)))
            else:
                countsMatTx_list.append(sample['NumReads'].reset_index(drop=True))
                abundanceMatTx_list.append(sample['TPM'].reset_index(drop=True))

            if varReduce and tx_out:
                varMatTx_list.append(pd.Series(boots_info['boots_var']))
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
        #return txi
    else:
        txiGene = summarizeToGene(txi, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance)
        
        #return txiGene
        txi = txiGene
    
    # Crear el objeto AnnData
    adata = ad.AnnData(X=txi['counts'].T.astype(np.float32))
    
    # Agregar la abundancia como una capa
    adata.layers['abundance'] = txi['abundance'].T.astype(np.float32)
    
    # Agregar las longitudes como una capa
    adata.layers['length'] = txi['length'].T.astype(np.float32)
    
    # Agregar información de las variables (genes/transcritos)
    adata.var_names = txi['counts'].index
    adata.var['gene_ids'] = adata.var_names
    
    # Agregar información de las observaciones (muestras)
    adata.obs_names = txi['counts'].columns
    adata.obs['sample'] = adata.obs_names
    
    # Agregar información adicional a uns
    adata.uns['countsFromAbundance'] = txi['countsFromAbundance']
    
    # Si hay datos de varianza, agregarlos como una capa
    if 'variance' in txi:
        adata.layers['variance'] = txi['variance'].T.astype(np.float32)
    
    # Si hay réplicas inferenciales, agregarlas como un array multidimensional en uns
    if 'infReps' in txi:
        adata.uns['infReps'] = {sample: reps.astype(np.float32) for sample, reps in txi['infReps'].items()}

    return adata