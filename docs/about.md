# Workflow
Here the different analysis pipelines, when using **`pymportx`**, are described.

## Previous Packages
---
As mentioned [before](https://pymportx.readthedocs.io/en/latest/#:~:text=Its%20upstream%20quantification%20methods%20(Salmon%2C%20Sailfish%2C%20Kallisto%2C%20and%20RSEM)), the **`pymportx`** package is compatible for fast gene count estimation with the following upstream trasncript quantification methods: Salmon, Sailfish, Kallisto, and RSEM. The package employment for each method is described in the [**Usage Tutorial**](https://pymportx.readthedocs.io/en/latest/#:~:text=as%20shown%20below%3A-,Usage%20tutorial,-Salmon).

## Downstream Analysis Packages
---

### PyDESeq2 and decoupleR

The `anndata` output from **`pymportx`** can be used to perform a Differential Expression Analysis using [PyDESeq2](https://pydeseq2.readthedocs.io/en/latest/). Other analyses, such as reads quality control or volcano plotting, can be also conducted using the [decoupleR](https://decoupler-py.readthedocs.io/en/latest/notebooks/bulk.html#Quality-control) package. Click on the links to access their documentation for more details.

Here is an example for Differential Expression Analysis using PyDESeq2:

```python
# Import PyDESeq2
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats
```
```python
# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(adata=adata,
    design_factors='condition',
    refit_cooks=True,
    inference=inference)
```

```pyhton
# Compute LFCs
dds.deseq2()
```
```python
# Extract contrast between COVID-19 vs normal
stat_res = DeseqStats(dds,
    contrast=["condition", 'treatment', 'control'],
    inference=inference)
```

```python
# Compute Wald test
stat_res.summary()
```

```python
# Extract results
results_df = stat_res.results_df
results_df
```


### PyWGCNA



```
wgcna = PyWGCNA.WGCNA(anndata=adata_cleaned)
 
wgcna.preprocess()
 
wgcna.findModules()
 
wgcna.analyseWGCNA()
```
