import requests, re
def get_uniprot(pdb_list):
    """
    Converts list of pdb to dict of uniprots {pdb: space_separated_uniprots}
    """
    result = {}
    exceptions = []

    for pdb in pdb_list:
        uniprots = []
        url2 = "https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=" + pdb
        try:
            uniprot = requests.get(url2)
        except requests.exceptions.RequestException:
            print(f"Exception for {pdb} (pdb -> uniprot)")
            exceptions.append(pdb)
        uniprots.extend(re.findall('dbAccessionId="([^"]{5,}?)"', uniprot.content.decode("utf8")))
        uniprots = list(set(uniprots))
        string = " ".join(uniprots)
        print(pdb, string)
        result[pdb] = string

    print(f"Warning! Uniprots not retrieved for {exceptions}")
    return result
    

def smiles2pdb(smiles_list):
    """
    Converts list of ligand smiles to dict with PDB of proteins which contains ligand {smiles: space_separated_pdbs}
    """
    exceptions = []

    result = {}
    for smi in smiles_list:

        url = "https://www.rcsb.org/pdb/rest/smilesQuery?smiles=" + smi + "&search_type=exact"
        try:
            pdbdir = requests.get(url)
        except requests.exceptions.RequestException:
            print(f"Exception for {smi} (smiles -> pdb)")
            exceptions.append(smi)
            continue
        pdbsF = re.findall('structureId="(.{4})"', pdbdir.content.decode("utf8"))  
        pdbsF = list(set(pdbsF))
        string = " ".join(pdbsF)
        print(string)
        result[smi] = string
    return result