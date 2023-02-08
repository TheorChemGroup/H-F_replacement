#!/usr/bin/env python
# coding: utf-8

# # Algorithm-1

# [Data base](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_28/)

# In[41]:


import re
import pandas as pd
pd.set_option('display.max_columns', None)
import subprocess
import chembl_credentials
import os
import psycopg2
import psycopg2.extras


CONNECTION = psycopg2.connect(database=chembl_credentials.database, 
                              user=chembl_credentials.user, 
                              password=chembl_credentials.password, 
                              host=chembl_credentials.host, 
                              port=chembl_credentials.port)


# In[42]:


os.makedirs("./alg1_files", exist_ok=True)


# In[43]:


QUERY_GET_MOLECULES = """
    SELECT 
      md.MOLREGNO,
      md.CHEMBL_ID,
      cs.CANONICAL_SMILES,
      array_agg(a.ASSAY_ID) assays
    FROM MOLECULE_DICTIONARY AS md
    JOIN COMPOUND_STRUCTURES AS cs
      ON md.MOLREGNO=cs.MOLREGNO
    JOIN ACTIVITIES AS a
      ON md.MOLREGNO=a.MOLREGNO
    WHERE a.ASSAY_ID IS NOT NULL
    AND cs.CANONICAL_SMILES ~ 'F([^emlr]|$)'
    GROUP BY md.MOLREGNO, cs.MOLREGNO
    ORDER BY md.MOLREGNO
    {}
    ;
"""


def get_molecules_with_fluorine(offset=0, limit='ALL'):
    query = QUERY_GET_MOLECULES.format(f'LIMIT {limit} OFFSET {offset}')
    
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
    cur.execute(query)
    rows = cur.fetchall()

    return rows


# In[9]:


REGEX_FLUORINE = re.compile(r'(\(F\)|F(?![emlr])|F$)')

def get_all_fluorine_indexes(smiles):
    finds = list(REGEX_FLUORINE.finditer(smiles))
    return finds


# In[10]:


def gen_next_smiles(finds, smiles):
    for match in finds:
        next_smiles = smiles[:match.start()] + smiles[match.end():]
        yield next_smiles


# In[11]:


def get_inchikey_by_smiles(smiles):
    list_files = subprocess.run(["obabel", "-:" + smiles, "-oinchikey"], stdout=subprocess.PIPE)
    return list_files.stdout.decode('utf-8').rstrip()


# In[12]:


QUERY_GET_MOLECULES_BY_INCHIKEY_AND_ASSAYS = """
    SELECT 
      md.MOLREGNO AS h_id,
      md.CHEMBL_ID AS h_chembl,
      cs.CANONICAL_SMILES h_smiles,
      a.ASSAY_ID AS assay_id,
      a.CHEMBL_ID AS assay_chembl,
      a.ASSAY_TYPE AS assay_type,
      td.TID AS target_id,
      td.CHEMBL_ID AS target_chembl,
      td.TARGET_TYPE AS target_type,
      ac.STANDARD_TYPE AS h_type,
      ac.STANDARD_RELATION AS h_relation,
      ac.STANDARD_VALUE AS h_value,
      ac.STANDARD_UNITS AS h_units
    FROM MOLECULE_DICTIONARY AS md
    JOIN COMPOUND_STRUCTURES AS cs
      ON md.MOLREGNO = cs.MOLREGNO
    JOIN ACTIVITIES AS ac
      ON md.MOLREGNO = ac.MOLREGNO
    JOIN ASSAYS AS a
      ON ac.ASSAY_ID = a.ASSAY_ID
    JOIN TARGET_DICTIONARY AS td
      ON a.TID = td.TID
    WHERE ac.ASSAY_ID IN %s
    AND cs.STANDARD_INCHI_KEY = %s
    ;
""" 


def get_molecule_by_inchikey_where_in_assays(inchikey, assays):
    query = QUERY_GET_MOLECULES_BY_INCHIKEY_AND_ASSAYS
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query, (tuple(assays), inchikey))
    molecule = cur.fetchall()

    return molecule


# In[13]:


QUERY_GET_MOLECULE_INFO = """
    SELECT 
      md.MOLREGNO AS f_id,
      md.CHEMBL_ID AS f_chembl,
      cs.CANONICAL_SMILES f_smiles,
      a.ASSAY_ID AS assay_id,
      a.CHEMBL_ID AS assay_chembl,
      a.ASSAY_TYPE AS assay_type,
      td.TID AS target_id,
      td.CHEMBL_ID AS target_chembl,
      td.TARGET_TYPE AS target_type,
      ac.STANDARD_TYPE AS f_type,
      ac.STANDARD_RELATION AS f_relation,
      ac.STANDARD_VALUE AS f_value,
      ac.STANDARD_UNITS AS f_units
    FROM MOLECULE_DICTIONARY AS md
    JOIN COMPOUND_STRUCTURES AS cs
      ON md.MOLREGNO = cs.MOLREGNO
    JOIN ACTIVITIES AS ac
      ON md.MOLREGNO = ac.MOLREGNO
    JOIN ASSAYS AS a
      ON ac.ASSAY_ID = a.ASSAY_ID
    JOIN TARGET_DICTIONARY AS td
      ON a.TID = td.TID
    WHERE cs.MOLREGNO = %s
    AND a.ASSAY_ID = %s
    AND ac.STANDARD_TYPE = %s
    ;
"""

def get_molecule_info(molregno, assay_id, standard_type):
    query = QUERY_GET_MOLECULE_INFO
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query, (molregno, assay_id, standard_type))
    row = cur.fetchall()

    return row


# In[ ]:


COLUMNS = ('f_id', 'f_chembl',
           'h_id', 'h_chembl',
           'f_smiles',
           'h_smiles',
           'assay_id', 'assay_chembl',
           'assay_type',
           'target_id', 'target_chembl',
           'target_type',
           'f_type',
           'f_relation',
           'f_value',
           'f_units',
           'h_type',
           'h_relation',
           'h_value',
           'h_units',
           'result'
           )

def get_molecules_with_other_smiles(molecules):
    molecules_info = pd.DataFrame(columns=COLUMNS)
    for i, molecule in enumerate(molecules):
        finds = get_all_fluorine_indexes(molecule.canonical_smiles)
        for other_smiles in gen_next_smiles(finds, molecule.canonical_smiles):
            other_molecules = get_molecule_by_inchikey_where_in_assays(
                get_inchikey_by_smiles(other_smiles),
                molecule.assays)
            for other_molecule in other_molecules:
                
                fs_info = get_molecule_info(molecule.molregno, other_molecule['assay_id'], other_molecule['h_type'])
                for f_info in fs_info:
                    
                    molecules_info = molecules_info.append(
                        pd.Series({**f_info, **other_molecule}, index=molecules_info.columns), ignore_index=True)

    return molecules_info


# In[ ]:


import csv
def to_csv(obj, filename):
    with open(filename,'w') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(obj[0]._fields)
        for row in obj:
            csv_out.writerow(row)


# In[ ]:


def main(molecules, offset=0):
    batch_size = 1000
    
    i = 0
    for i in range(offset, int(len(molecules) / batch_size)):
        start = i*batch_size
        end = start + batch_size
        
        mol = get_molecules_with_other_smiles(molecules[start:end])
        mol.to_csv(f'alg1_files/mol_{start}_{end-1}.csv')
        print(f'saved alg1_files/mol_{start}_{end-1}.csv')
        
    mol = get_molecules_with_other_smiles(molecules[i*batch_size:])
    mol.to_csv(f'alg1_files/mol_{i*batch_size}_{len(molecules)-1}.csv')


# In[ ]:


molecules = get_molecules_with_fluorine()


# In[ ]:


main(molecules)

