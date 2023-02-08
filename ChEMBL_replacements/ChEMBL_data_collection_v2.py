#!/usr/bin/env python
# coding: utf-8

# # Algorithm-2

# [Data base](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_28/)

# In[ ]:


import re
import pandas as pd
pd.set_option('display.max_columns', None)
import subprocess
import psycopg2
import psycopg2.extras
import glob
import os
import chembl_credentials


CONNECTION = psycopg2.connect(database=chembl_credentials.database, 
                              user=chembl_credentials.user, 
                              password=chembl_credentials.password, 
                              host=chembl_credentials.host, 
                              port=chembl_credentials.port)


# In[ ]:


os.makedirs("./molfiles", exist_ok=True)
os.makedirs("./sdf", exist_ok=True)
os.makedirs("./smi", exist_ok=True)
os.makedirs("./alg2_files", exist_ok=True)


# In[ ]:


QUERY_GET_ASSAYS = """
    SELECT
      ac.ASSAY_ID AS assay_id  
    FROM ACTIVITIES as ac
    WHERE ac.MOLREGNO IS NOT NULL
    ;
"""

def get_assays():
    query = QUERY_GET_ASSAYS
    cur = CONNECTION.cursor()
    cur.execute(query)
    rows = cur.fetchall()

    return rows


# In[ ]:


assays = list(set(get_assays()))


# In[ ]:


QUERY_GET_LIST_OF_MOLECULES_FROM_ASSAY = """
    SELECT
    ac.MOLREGNO
    FROM ACTIVITIES as ac
    WHERE ac.ASSAY_ID = %s
    ;
"""

def get_list_of_molecules_from_assay(assay_id):
    query = QUERY_GET_LIST_OF_MOLECULES_FROM_ASSAY
    cur = CONNECTION.cursor()
    cur.execute(query, (assay_id,))
    molecules = cur.fetchall()
    
    return molecules


# In[ ]:


QUERY_GET_MOLECULE_WITH_FLUORINE_IN_ASSAY = """
    SELECT 
      cs.MOLREGNO
    FROM COMPOUND_STRUCTURES as cs
    WHERE cs.MOLREGNO = %s
    AND cs.CANONICAL_SMILES ~ 'F([^emlr]|$)'
    ;
"""

def get_molecule_with_fluorine_in_assay(molregno):
    query = QUERY_GET_MOLECULE_WITH_FLUORINE_IN_ASSAY
    cur = CONNECTION.cursor()
    cur.execute(query, (molregno,))
    molecules = cur.fetchone()
    
    return molecules


# In[ ]:


QUERY_GET_MOLFILE_BY_MOLREGNO = """
    SELECT
      cs.MOLFILE
    FROM COMPOUND_STRUCTURES AS cs
    WHERE cs.MOLREGNO = %s
    ;
"""

def get_molfile_by_molregno(molregno):
    query = QUERY_GET_MOLFILE_BY_MOLREGNO
    cur = CONNECTION.cursor()
    cur.execute(query, (molregno,))
    row = cur.fetchone()
    
    return row


# In[ ]:


def save_as(molregno):
    with open(f'molfiles/{molregno}.mol','w') as file:
        print(get_molfile_by_molregno(molregno)[0], file=file)


# In[ ]:


def molfile_with_fluorine_in_all_assays(assays):
    molecules_with_fluorine = []
    
    for assay_id in assays:
        all_molecules = get_list_of_molecules_from_assay(assay_id)
        
        for molregno in all_molecules:
            with_fluorine = get_molecule_with_fluorine_in_assay(molregno[0])
            if with_fluorine:
                molecules_with_fluorine.append(with_fluorine[0])
                save_as(with_fluorine[0])     
                
    return molecules_with_fluorine


# In[ ]:


def all_molfiles(assays, offset=0):
    batch_size = 1000
    
    i = 0
    for i in range(offset, int(len(assays) / batch_size)):
        start = i * batch_size
        end = start + batch_size
        
        molfile_with_fluorine_in_all_assays(assays[start:end])
        print(f'saved {start}_{end-1}.mol')


# In[ ]:


all_molfiles(assays)


# In[ ]:


def get_sdf_by_mol(list_of_molregno):
    for i, molregno in enumerate(list_of_molregno):
        subprocess.run(['f2h.exe', 'molfiles/' + str(molregno) + '.mol', 'sdf/' + str(molregno) + '.sdf'])


# In[ ]:


list_of_molregno = []
for root, dirs, files in os.walk("./molfiles"):  
    for filename in files:
        if filename.count('.mol') == 1:
            list_of_molregno.append(filename[:-4])


# In[ ]:


len(list_of_molregno)


# In[ ]:


get_sdf_by_mol(list_of_molregno)


# In[ ]:


def gen_inchikey_by_sdf(sdf_file):
    list_files_1 = subprocess.run(["obabel", "sdf/" + sdf_file + ".sdf", "-O", "smi/" + sdf_file + ".smi"])
    list_files_2 = subprocess.run(["obabel", "smi/" + sdf_file + ".smi", "-oinchikey"], stdout=subprocess.PIPE)
    inchikey_list = list_files_2.stdout.decode('utf-8').rstrip().split('\n')
    
    return set(inchikey_list)


# In[ ]:


QUERY_GET_F_MOLECULE_BY_MOLREGNO = """
    SELECT 
      md.MOLREGNO AS f_id,
      md.CHEMBL_ID AS f_chembl,
      cs.CANONICAL_SMILES f_smiles,
      a.ASSAY_ID AS f_assay_id,
      a.CHEMBL_ID AS f_assay_chembl,
      a.ASSAY_TYPE AS f_assay_type,
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
    ;
"""


def get_f_molecule_by_molregno(molregno):
    query = QUERY_GET_F_MOLECULE_BY_MOLREGNO
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query, (molregno,))
    molecule = cur.fetchall()

    return molecule


# In[ ]:


QUERY_GET_MOLECULES_BY_INCHIKEY_AND_ASSAYS = """
    SELECT 
      md.MOLREGNO AS h_id,
      md.CHEMBL_ID AS h_chembl,
      cs.CANONICAL_SMILES h_smiles,
      a.ASSAY_ID AS h_assay_id,
      a.CHEMBL_ID AS h_assay_chembl,
      a.ASSAY_TYPE AS h_assay_type,
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
    WHERE a.ASSAY_ID = %s
    AND cs.STANDARD_INCHI_KEY = %s
    AND ac.STANDARD_TYPE = %s    
    ;
"""


def get_molecule_by_inchikey_where_in_assays(assay_id, inchikey, standard_type):
    query = QUERY_GET_MOLECULES_BY_INCHIKEY_AND_ASSAYS
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute(query, (assay_id, inchikey, standard_type))
    molecule = cur.fetchall()

    return molecule


# In[ ]:


COLUMNS = ('f_id', 'f_chembl',
           'h_id', 'h_chembl',
           'f_smiles',
           'h_smiles',
           'f_assay_id', 'f_assay_chembl',
           'f_assay_type',
           'h_assay_id', 'h_assay_chembl',
           'h_assay_type',
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


def new_alg(list_of_f_molregno):
    data = pd.DataFrame(columns=COLUMNS)
    
    for f_molregno in list_of_f_molregno:
        f_molecule_info = get_f_molecule_by_molregno(f_molregno) #список записей
        
        for f_info in f_molecule_info:
            f_standard_type = f_info['f_type']
            f_assay_id = f_info['f_assay_id']

            h_inchikey_list = gen_inchikey_by_sdf(f_molregno)

            for h_inchikey in h_inchikey_list:
                h_molecule_info = get_molecule_by_inchikey_where_in_assays(f_assay_id, 
                                                                           h_inchikey,
                                                                           f_standard_type)
                for h_info in h_molecule_info:
                    data.loc[len(data)] = pd.Series({**f_info, **h_info}, 
                                                     index=data.columns)

    return data


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
        
        mol = new_alg(molecules[start:end])
        mol.to_csv(f'alg2_files/mol_{start}_{end-1}.csv')
        print(f'saved alg2_files/mol_{start}_{end-1}.csv')
        
    mol = new_alg(molecules[i*batch_size:])
    mol.to_csv(f'alg2_files/mol_{i*batch_size}_{len(molecules)-1}.csv')
    print(f'saved alg2_files/mol_{i*batch_size}_{len(molecules)-1}.csv')


# In[ ]:


molecules_with_fluorine = []
for filename in glob.iglob('sdf/*.sdf'):
    molecules_with_fluorine.append(filename[4:-4])


# In[ ]:


main(molecules_with_fluorine)

