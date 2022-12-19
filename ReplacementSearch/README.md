# ReplacementSearch
Code of generating a set of 16 protein-ligand complexes containing a bioisosteric hydrogen-fluorine substitution was obtained from the RCSB Protein Data Bank.

## Warning
Code were actual at February 2020
https://www.rcsb.org/news/5fc9176809ae2a096d081e28

## Installation
Our code requires: 
- [Python 3](https://www.python.org/downloads/)
- [Jupyter Notebook](https://jupyter.org/install)

Install the dependencies:

```sh
pip install -r requirements.txt
```

## Usage examples
1. Run Jupyter Notebook: 
```sh 
jupyter notebook Data_Analysis.ipynb
```
2. Run all cells in the notebook

## File organization

| Filename | Comment |
| ------ | ------ |
| `requirements.txt` | a file with requirements |
| `pdb2uniprot.py` | a python file additional functions for Data_Analysis.ipynb |
| `Data_Analysis.ipynb` | a Jupyter Notebook with procedure of generating a test set from RCSB Protein Data Bank. (METHODS)|