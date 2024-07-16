# pymportx
`pymportx` is a Python package designed for fast gene count estimation using transcript quantification files generated by pseudoalignment or quasi-mapping tools. It is a Python adaptation of the widely-used [`tximport`](https://bioconductor.org/packages/release/bioc/html/tximport.html) R package from Bioconductor.

## Documentation
Documentation is made available at [https://pymportx.readthedocs.io](https://pymportx.readthedocs.io/en/latest/). 

## Example
```python
from pymportx import salmon

results = salmon.read_salmon(salmon_folder_paths,
                             tx_out=False,
                             tx2gene=tx2gene_file_path,
                             countsFromAbundance='no')
```

## Citation

If you use `pymportx` in your research, please cite the following paper:

Pena Gonzalez, P., Lozano-Paredes, D., Rojo, J. L., Bote-Curiel, L., & Sanchez-Arevalo Lobo, V. J. (2024). Pymportx: Facilitating Next-Generation Transcriptomics Analysis in Python. *bioRxiv*. doi:10.1101/2024.07.12.598873

## License
This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact
For more information, you can contact us through victor.sanchezarevalo@ufv.es or luis.bote@urjc.es.
