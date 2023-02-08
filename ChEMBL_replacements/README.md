# ChEMBL_replacements
The code used to match pairs of molecules with fluorine and their hydrogen analogues (the hydrogen analogue is a bioisioster in which the fluorine atom is replaced by a hydrogen atom).

For initial assessment of the effect of the bioisosteric replacement of hydrogen atom by a fluorine (further H → F) we evaluated the ratio between the values of biological activity constants (IC50, or EC50, or Kd, or Ki) of compounds with fluorine and the values of constants corresponding to their hydrogen-bearing counterparts.

In short, we used a data comparison algorithm that allows us to compare the biological activity of all compounds with fluorine available in ChEMBL database with their hydrogen analogues.

## Installation
Our code requires: 
- [Python 3](https://www.python.org/downloads/)
- [PostgreSQL](https://www.postgresql.org/download/)
- [CHEMBL28](https://chembl.gitbook.io/chembl-interface-documentation/downloads) for Postgres (`chembl_28_postgresql.tar.gz`)
- [Open Babel](https://openbabel.org/docs/dev/Installation/install.html) v3.1.0

Install the dependencies:

```sh
pip install -r requirements.txt
```

## Usage examples
1. Replace the data in the `chembl_credentials.py` file with your database connection parameters.
2. Run 3 files sequentially: 
```sh 
python3 ChEMBL_data_collection_v1.py 
```
```sh 
python3 ChEMBL_data_collection_v2.py 
```
```sh 
python3 ChEMBL_data_cleaning.py 
```

## File organization

| Filename | Comment |
| ------ | ------ |
| `requirements.txt` | a file with requirements for running algorithm files |
| `chembl_credentials.py` | a file with database connection parameters (port, host, user, password, database) |
| `ChEMBL_data_collection_v1.py` | a python script executing Algorithm 1 (METHODS)|
| `ChEMBL_data_collection_v2.py` | a python script executing Algorithm 2 (METHODS)|
| `ChEMBL_data_cleaning.py` | a python script that cleans the data and builds a reliable dataset (METHODS)|
| `f2h.exe` | a Windows program taking a MOL file of F-substituted ligand and producing an SDF file of H-substituted ligand; used in Algorithm_2.py  |
| `ReliableDataSet.csv` | a reliable dataset |
