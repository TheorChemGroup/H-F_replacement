#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import glob

import psycopg2
import psycopg2.extras

import chembl_credentials


CONNECTION = psycopg2.connect(database=chembl_credentials.database, 
                              user=chembl_credentials.user, 
                              password=chembl_credentials.password, 
                              host=chembl_credentials.host, 
                              port=chembl_credentials.port)


# In[92]:


alg1 = pd.DataFrame()
for filename in glob.iglob('alg1_files/mol_*.csv'):
    
    alg1 = pd.concat([alg1, pd.read_csv(filename, index_col=0)])

alg1 = alg1.drop_duplicates()


# In[93]:


alg2 = pd.DataFrame()
for filename in glob.iglob('alg2_files/mol_*.csv'):
    
    alg2 = pd.concat([alg2, pd.read_csv(filename, index_col=0)])

alg2 = alg2.drop_duplicates()


# In[94]:


alg1 = alg1.rename(columns={'f_type': 'standard_type'})
alg2 = alg2.rename(columns={'f_type': 'standard_type'})
alg2 = alg2.rename(columns={'f_assay_id': 'assay_id'})

alg1 = alg1.drop('h_type', axis=1)
alg2 = alg2.drop('h_type', axis=1)
alg2 = alg2.drop('h_assay_id', axis=1)


# In[95]:


f_ids1 = alg1['f_id'].tolist()
h_ids1 = alg1['h_id'].tolist()
assays1 = alg1['assay_id'].tolist()
types1 = alg1['standard_type'].tolist()
f_values1 = alg1['f_value'].tolist()
h_values1 = alg1['h_value'].tolist()

keys1 = list(zip(f_ids1, h_ids1, assays1, types1, f_values1, h_values1))
alg1['key'] = keys1


# In[96]:


f_ids2 = alg2['f_id'].tolist()
h_ids2 = alg2['h_id'].tolist()
assays2 = alg2['assay_id'].tolist()
types2 = alg2['standard_type'].tolist()
f_values2 = alg2['f_value'].tolist()
h_values2 = alg2['h_value'].tolist()

keys2 = list(zip(f_ids2, h_ids2, assays2, types2, f_values2, h_values2))
alg2['key'] = keys2


# In[97]:


intersection = set(alg1['key'].tolist()).intersection(set(alg2['key'].tolist()))


# In[98]:


df = alg1.loc[alg1['key'].isin(intersection)]


# In[99]:


ok_types = ['Ki', 'Kd', 'IC50', 'EC50']
df = df.loc[df['standard_type'].isin(ok_types)]


# In[100]:


df = df.loc[(df['f_value'] > 0.1) & (df['f_value'] < 1000000)]
df = df.loc[(df['h_value'] > 0.1) & (df['h_value'] < 1000000)]


# In[101]:


df = df.loc[(df['f_units'] == 'nM') & (df['h_units'] == 'nM')]


# In[102]:


df = df.loc[(df['f_relation'] == '=') & (df['h_relation'] == '=')]


# In[103]:


df = df.loc[df['target_type'].isin(['SINGLE PROTEIN', 'PROTEIN COMPLEX'])]


# In[151]:


QUERY_GET_ADDITIONAL_PARAMETERS_BY_ASSAY_ID = """
    SELECT
      a.ASSAY_ID AS assay_id,
      a.ASSAY_STRAIN AS strain_id,
      a.ASSAY_TISSUE AS tissue_id,
      a.ASSAY_CELL_TYPE AS celltype,
      a.ASSAY_SUBCELLULAR_FRACTION AS fraction,
      a.CONFIDENCE_SCORE AS confidence,
      a.RELATIONSHIP_TYPE AS relationship,
      a.DOC_ID AS doc_id
    FROM ASSAYS AS a
    WHERE a.ASSAY_ID = %s
    ;
"""

def get_additional_parameters_by_assay_id(assay_id):
    query = QUERY_GET_ADDITIONAL_PARAMETERS_BY_ASSAY_ID
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
    cur.execute(query, (assay_id,))
    row = cur.fetchone()

    return row


# In[171]:


COLUMNS_ADDITIONAL = ['assay_id', 'strain_id', 'tissue_id', 'celltype', 'fraction', 
                      'confidence', 'relationship', 'doc_id']

def export_additional_parameters(list_of_assay_id):
    res = pd.DataFrame(columns=COLUMNS_ADDITIONAL)
    
    for assay_id in list_of_assay_id:
        additional_parameters = get_additional_parameters_by_assay_id(assay_id)
            
        res.loc[len(res)] = additional_parameters
    
    return res


# In[ ]:


assays = list(set(df['assay_id'].tolist()))


# In[156]:


additional_params_assay = export_additional_parameters(assays)


# In[161]:


df = df.merge(additional_params_assay, on='assay_id').drop_duplicates()


# In[167]:


QUERY_GET_ADDITIONAL_PARAMETERS_BY_DOC_ID = """
    SELECT
      d.DOC_ID AS doc_id,
      d.DOC_TYPE AS doc_type
    FROM DOCS AS d
    WHERE d.DOC_ID = %s
    ;
"""

def get_additional_parameters_by_doc_id(doc_id):
    query = QUERY_GET_ADDITIONAL_PARAMETERS_BY_DOC_ID
    cur = CONNECTION.cursor(cursor_factory=psycopg2.extras.NamedTupleCursor)
    cur.execute(query, (doc_id,))
    row = cur.fetchone()

    return row


# In[169]:


COLUMNS_ADDITIONAL_DOCS = ['doc_id', 'doc_type']

def export_additional_parameters_docs(list_of_doc_id):
    res = pd.DataFrame(columns=COLUMNS_ADDITIONAL_DOCS)
    
    for doc_id in list_of_doc_id:
        additional_parameters_docs = get_additional_parameters_by_doc_id(doc_id)
            
        res.loc[len(res)] = additional_parameters_docs
    
    return res


# In[ ]:


docs = list(set(df['doc_id'].tolist()))


# In[172]:


additional_params_docs = export_additional_parameters_docs(docs)


# In[175]:


df = df.merge(additional_params_docs, on='doc_id').drop_duplicates()


# In[177]:


df = df.loc[df['confidence'].isin([7, 8, 9])]


# In[179]:


df = df.loc[df['relationship'] == 'D']


# In[181]:


df = df.loc[df['doc_type'] == 'PUBLICATION']


# In[183]:


df = df.loc[df['tissue_id'].isna()]
df = df.loc[df['strain_id'].isna()]
df = df.loc[df['celltype'].isna()]
df = df.loc[df['fraction'].isna()]


# In[188]:


df['ratio'] = df['f_value'] / df['h_value']


# In[212]:


df.to_csv('ReliableDataSet.csv')


# In[2]:


df = pd.read_csv('ReliableDataSet.csv', index_col=0)


# In[7]:


print(df.shape[0])
print(df.loc[df['ratio'] > 4].shape[0], round((df.loc[df['ratio'] > 4].shape[0] / df.shape[0]* 100), 1))
print(df.loc[df['ratio'] < 0.25].shape[0], round((df.loc[df['ratio'] < 0.25].shape[0] / df.shape[0]* 100), 1))

