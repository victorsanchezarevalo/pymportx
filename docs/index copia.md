![Foto](Logo.png)

Here the main arguments of the package functions are described:

## Salmon, Sailfish or kallisto

 ``salmon.read_salmon()``, ``sailfish.read_sailfish()``, `kallisto.read_kallisto()`

This main functions have the following arguments:

| Argument | Description |
|:---------|:-----------:|
|`folders` | A list of strings with the folder path for each of the samples.|        
|`tx_out` | Boolean argument. Default is False for gene-level output. Set to True for transcript-level output.  |
|`tx2gene` | A two-column .csv file containing gene annotations: transcript ID in the first column and and gene ID in the second column.| 
|`countsFromAbundance` | Could be set to either: **" no "** (default), **" scaledTPM "**, **" lengthScaledTPM "**, **" dtuScaledTPM "**. See [**countsFromAbundance**](http://127.0.0.1:8000/index%20copia/#:~:text=*countsFromAbundance%3A) for more detail.|  
|`dropInfReps` | Whether to skip inferential replicates read or not (default is False).|
|`varReduce` | Whether to condense per-sample inferential replicated into a matrix displaying sample variances (default is False).|
|`infRepStat`| A predefined function to operate over rows of inferential replicates (default is median over rows).|
|`ignoreTxVersion` | Whether to ignore transcript isoforms by removing the version number from the transcriptID after the period ' . ' . |  
|`ignoreAfterBar` |Whether to ignore transcriptID characters after the bar ' / '. | 


### *countsFromAbundance:

This ``countsFromAbundance`` argument could be set to either:

* **" no "** (default) : Determining whether to produce estimated counts using abundance estimations.

* **" scaledTPM "** : Obtain estimated counts scaled up to library size.

* **" lengthScaledTPM "** : Adjusted by utilizing the mean transcript length across samples,    followed by the library size.

* **" dtuScaledTPM "** : Obtain estimated counts scaled by employing the median length of transcripts within gene isoforms, followed by the library size.

## RSEM

The arguments for ``rsem.read_rsem()`` are detailed in the following table:

| Argument | Description |
|:---------|:-----------:|
|`folders` | A list of strings with the quantification file path for each of the samples.|     
|`tx_in` | Boolean argument. Default is True for trasncript-level input. Set to False for transcript-level output. |     
|`tx_out` | Boolean argument. Default is False for gene-level output. Set to True for transcript-level output.  |
|`tx2gene` | A two-column .csv file containing gene annotations: transcript ID in the first column and and gene ID in the second column.| 
|`countsFromAbundance` | Could be set to either: **" no "** (default), **" scaledTPM "**, **" lengthScaledTPM "**, **" dtuScaledTPM "**. See [**countsFromAbundance**](http://127.0.0.1:8000/index%20copia/#:~:text=*countsFromAbundance%3A) for more detail.| 
|`ignoreTxVersion` | Whether to ignore transcript isoforms by removing the version number from the transcriptID after the period ' . ' . |  
|`ignoreAfterBar` |Whether to ignore transcriptID characters after the bar ' / '. | 

### *countsFromAbundance:

This ``countsFromAbundance`` argument could be set to either:

* **" no "** (default) : Determining whether to produce estimated counts using abundance estimations.

* **" scaledTPM "** : Obtain estimated counts scaled up to library size.

* **" lengthScaledTPM "** : Adjusted by utilizing the mean transcript length across samples,    followed by the library size.

* **" dtuScaledTPM "** : Obtain estimated counts scaled by employing the median length of transcripts within gene isoforms, followed by the library size.



