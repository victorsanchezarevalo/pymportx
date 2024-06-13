import pandas as pd
import numpy as np
import json
import gzip
import h5py

# Function to read data from an h5 file
def read_kallisto_h5(fpath):
    with h5py.File(fpath, 'r') as f:
        counts = np.array(f['est_counts'])
        ids = np.array(f['aux/ids'], dtype=str)
        efflens = np.array(f['aux/eff_lengths'])
        
        # Check if the lengths match
        if not (len(counts) == len(ids) == len(efflens)):
            raise ValueError("Lengths of counts, ids, and eff_lengths do not match.")

        result = pd.DataFrame({
            'target_id': ids,
            'eff_length': efflens,
            'est_counts': counts,
        })

        # Calculate TPM
        normfac = (1e+06) / np.sum(result['est_counts'] / result['eff_length'])
        result['tpm'] = normfac * (result['est_counts'] / result['eff_length'])

    return result


# Function to read inferential replicates from a Kallisto h5 file
def readInfRepKallisto(fpath):
    with h5py.File(fpath,'r') as h5_file:
        boots = h5_file['bootstrap']
        numBoot = len(boots.keys())
        numTxp  = len(boots['bs0'])
        bootMat = np.empty((numTxp, numBoot))

        for bsn in range(numBoot):
            bootMat[:, bsn] = boots[f'bs{bsn}'][()]
        vars = np.var(bootMat, ddof=1, axis=1)
        
    return {"boots_vars": vars, "boots": bootMat}
        
# Function to read inferential replicates from a Salmon or Sailfish file
def readInfRepFish(fish_type, meta_info_path, bootstraps_path):
    # Read 'meta_info.json' file in auxiliary directory
    with open(meta_info_path, 'r') as meta_info_file:
        meta_info = json.load(meta_info_file)

    # Check Salmon version
    if fish_type == 'salmon':
        if 'salmon_version' in meta_info:
            required_version = (0, 8, 0)
            salmon_version = tuple(map(int, meta_info["salmon_version"].split('.')))
            if salmon_version < required_version:
                raise ValueError("Salmon version must be at least 0.8.0")
        else:
            raise ValueError("Salmon version not found in meta_info.json")
    # Check Sailfish version
    if fish_type == 'sailfish':
        if 'sailfish_version' in meta_info:
            required_version = (0, 9, 0)
            sailfish_version = tuple(map(int, meta_info["sailfish_version"].split('.')))
            if sailfish_version < required_version:
                raise ValueError("Sailfish version must be at least 0.9.0")
        else:
            raise ValueError("Sailfish version not found in meta_info.json")

    # Read bootstraps file
    try:
        with gzip.open(bootstraps_path, 'rb') as boot_con:
            boots = np.frombuffer(boot_con.read(), dtype=np.float64)
            if len(boots) == 0:
                boots = np.frombuffer(boot_con.read(), dtype=np.int32)           
    except Exception:
        raise ValueError(f"Error reading bootstraps file {bootstraps_path}")
 
    # Reshape the bootstraps
    boots = np.transpose(boots.reshape(meta_info["num_bootstraps"], meta_info["num_targets"]))

    # Calculate the variance of the bootstraps
    boots_var = np.var(boots, ddof=1, axis=1)

    return {'boots_var': boots_var, 'boots': boots}

# Function to summarize transcript-level data to gene-level data
def summarizeToGene(object, tx2gene, varReduce, ignoreTxVersion, ignoreAfterBar, countsFromAbundance):
    
    abundanceMatTx = object['abundance']
    countsMatTx = object['counts']
    lengthMatTx = object['length']

    txId = abundanceMatTx.index.tolist()

    # Check if all indices are the same in 'counts' and 'length' DataFrames
    if txId != countsMatTx.index.tolist():
        raise ValueError("Row indices of counts and abundance do not match!")

    if txId != lengthMatTx.index.tolist():
        raise ValueError("Row indices of length and abundance do not match!")

    # Need to associate tx to genes
    if tx2gene is not None:
        # Strip dots or bars and all remaining characters from the rownames of matrices
        if ignoreTxVersion:
            txId = [tx.split('.')[0] for tx in txId]
        elif ignoreAfterBar:
            txId = [tx.split('|')[0] for tx in txId]
        
    tx2gene = read_tx2gene(tx2gene)
    tx2gene = tx2gene[tx2gene['Tx'].isin(txId)]

    # Calculate the number of missing transcripts
    ntxmissing = sum(~pd.Series(txId).isin(tx2gene['Tx']))
    if ntxmissing > 0:
        print("Transcripts missing from tx2gene:", ntxmissing)

    txId = pd.Series(txId)
    sub_idx = np.array(txId.isin(tx2gene['Tx']))
    abundanceMatTx = abundanceMatTx[sub_idx]
    countsMatTx = countsMatTx[sub_idx]
    lengthMatTx = lengthMatTx[sub_idx]
   
    txId = txId[sub_idx]
    tx2gene_mapping = tx2gene.set_index('Tx')['Gene']
    
    # Map txId to gene using the dictionary
    geneId = tx2gene_mapping.loc[txId].tolist()

    # Summarize abundance
    print("Summarizing abundance")
    # argumento groups en groupby
    abundanceMat = abundanceMatTx.groupby(geneId).sum()
    # Summarize counts
    print("Summarizing counts")
    countsMat = countsMatTx.groupby(geneId).sum()
    print("Summarizing length")
    # Calculate the weighted average of transcript length
    weightedLength = (abundanceMatTx * lengthMatTx).groupby(geneId).sum()
    abundanceMatNaN = abundanceMat.replace(0, np.nan)
    weightedLength = weightedLength.replace(0, np.nan)

    # Normalize by abundance to remove length bias
    lengthMat = weightedLength.div(abundanceMatNaN, axis=0)

    # Calculate a simple average transcript length
    # Average the transcript lengths over samples
    aveLengthSamp = lengthMatTx.mean(axis=1)

    # Average the lengths within genes (not weighted by abundance)
    aveLengthSampGene = aveLengthSamp.groupby(geneId).mean()
    
    if isinstance(aveLengthSampGene, pd.DataFrame):
        if not all(aveLengthSampGene.columns == lengthMat.index):
                raise ValueError("Column names of aveLengthSampGene do not match row names of lengthMat")
        # If aveLengthSampGene is a Series
        elif isinstance(aveLengthSampGene, pd.Series):
            if not all(aveLengthSampGene.index == lengthMat.index):
                raise ValueError("Index of aveLengthSampGene does not match row names of lengthMat")
        else:
            raise ValueError("aveLengthSampGene must be either a DataFrame or a Series")
        
    lengthMat = replaceMissingLength(lengthMat, aveLengthSampGene)

    if countsFromAbundance != "no":
        countsMat = makeCountsFromAbundance(countsMat, abundanceMat, lengthMat, countsFromAbundance)
            
    out = {
    "abundance": abundanceMat,
    "counts": countsMat,
    "length": lengthMat,
    "countsFromAbundance": countsFromAbundance
    }

    if "infReps" in object:
        print("Summarizing inferential replicates")
        infReps = {}
        for key, value in object["infReps"].items(): 
            infRepMatTx = value[sub_idx]
            infRepMatdf = pd.DataFrame(data=infRepMatTx, index=abundanceMatTx.index)
            infReps[key] = infRepMatdf.groupby(geneId).sum()
            if varReduce:
                vars = np.var(infReps[key], ddof=1, axis=1)
                if "variance" not in out:
                    out["variance"] = pd.DataFrame(index=infReps[key].index)  # Initialize DataFrame
                out["variance"][key] = vars 
            else:
                out["infReps"] = infReps
    return out 

# Function to modify counts based on abundance and length
def makeCountsFromAbundance(countsMat, abundanceMat, lengthMat, countsFromAbundance):

    # Convert countsFromAbundance to lowercase
    countsFromAbundance = countsFromAbundance.lower()

    # Sum of counts per sample
    counts_sum = countsMat.sum(axis=0)

    if countsFromAbundance == "lengthscaledtpm":
        #new_counts = abundanceMat * lengthMat.mean(axis=1)
        aux_meanL = lengthMat.mean(axis=1)
        new_counts_array = abundanceMat.values * aux_meanL.values[:,np.newaxis]
        new_counts = pd.DataFrame(new_counts_array, columns = abundanceMat.columns, index = abundanceMat.index)
        #new_counts = abundanceMat.mul(aux_meanL, axis=1)

    elif countsFromAbundance == "scaledtpm":
        new_counts = abundanceMat

    else:
        raise ValueError("Expecting 'lengthScaledTPM' or 'scaledTPM'")

    new_sum = new_counts.sum(axis=0)

    countsMat = (new_counts.T.mul(counts_sum / new_sum, axis=0)).T

    return countsMat


#def replaceMissingLength(lengthMat, aveLengthSampGene):
#    nan_rows = np.where(np.any(np.isnan(lengthMat), axis=1))[0]
#    if len(nan_rows) > 0:
#        for i in nan_rows:
#            if np.all(np.isnan(lengthMat.iloc[i])):
#                # If all samples have 0 abundances for all tx, use the simple average
#                lengthMat.iloc[i] = aveLengthSampGene.iloc[i]
#            else:
#                # Otherwise use the geometric mean of the lengths from the other samples
#                idx = np.isnan(lengthMat.iloc[i])
#                lengthMat.loc[lengthMat.index[i], idx] = np.exp(np.mean(np.log(lengthMat.loc[lengthMat.index[i], np.logical_not(idx)]), axis=0))
#    return lengthMat


def replaceMissingLength(lengthMat, aveLengthSampGene):
    # Create a boolean mask for rows with NaN
    has_nan = lengthMat.isna().any(axis=1)
    
    # Create a boolean mask for rows with all NaN
    all_nan = lengthMat.isna().all(axis=1)
    
    # Replace rows with all NaN with aveLengthSampGene
    if all_nan.any():
        # Expand aveLengthSampGene to match the shape of lengthMat
        ave_length_expanded = np.tile(aveLengthSampGene.loc[all_nan].values[:, np.newaxis], (1, lengthMat.shape[1]))
        lengthMat.loc[all_nan] = ave_length_expanded
    
    # Create a boolean mask for rows with some NaN (but not all)
    some_nan = has_nan & ~all_nan
    
    if some_nan.any():
        # Get indices of non-NaN columns for each row
        valid_indices = lengthMat.notna().to_numpy()
        
        # Calculate the geometric mean of non-NaN values
        log_vals = np.log(lengthMat.to_numpy())
        log_means = np.sum(np.where(valid_indices, log_vals, 0), axis=1) / valid_indices.sum(axis=1)
        geom_means = np.exp(log_means)
        
        # Expand geom_means to match the shape of lengthMat
        geom_means_expanded = np.tile(geom_means[some_nan][:, np.newaxis], (1, lengthMat.shape[1]))
        
        # Replace NaN with the geometric mean
        lengthMat.loc[some_nan] = np.where(lengthMat.loc[some_nan].isna(), geom_means_expanded, lengthMat.loc[some_nan])
    
    return lengthMat


def read_tx2gene(tx2gene):
    # Read tx2gene data from CSV file
    tx2gene_df = pd.read_csv(tx2gene)

    # Check if the file has the shape (n, 2)
    if tx2gene_df.shape[1] != 2:
        raise ValueError("File must be a two-column CSV file with the mapping from transcript id to gene id. First column in the file is transcript id and second column is gene id")

    # Rename columns
    tx2gene_df.columns = ['Tx', 'Gene']

    # Remove duplicated transcript rows
    if tx2gene_df.duplicated('Tx').any():
        print("Removing duplicated transcript rows in 'tx2gene' file")
        tx2gene_df = tx2gene_df[~tx2gene_df.duplicated('Tx')]

    return tx2gene_df


def medianLengthOverIsoform(length, tx2gene, ignoreTxVersion, ignoreAfterBar):
    txId = length.index.tolist()
    
    if ignoreTxVersion:
        txId = [tx.split('.')[0] for tx in txId]
    elif ignoreAfterBar:
        txId = [tx.split('|')[0] for tx in txId]
    
    tx2gene = read_tx2gene(tx2gene)

    if not all(tx2gene['Tx'].isin(txId)):
        raise ValueError("Not all transcript IDs are present in 'tx2gene' file")
    
    tx2gene = tx2gene[tx2gene['Tx'].isin(txId)]
    
    tx2gene_mapping = tx2gene.set_index('Tx')['Gene']
    geneId = tx2gene_mapping.loc[txId].tolist()

    # average the lengths
    ave_len = length.mean(axis=1)
    
    # median over isoforms
    med_len = ave_len.groupby(geneId).median()
    
    one_sample = med_len.reindex(geneId).values
    
    return pd.DataFrame(np.tile(one_sample, (length.shape[1], 1)).T, index=length.index, columns=length.columns)


def process_infRepStat(infRepStat):
    # Change infRepStat name to its function. If None, keep it as None
    if infRepStat is None:
        pass
    elif infRepStat == "median":
        infRepStat = lambda x: np.median(x, axis=1)

    return infRepStat
