import os
import pandas as pd
try:
    from commands import  getstatusoutput
except ImportError:
    from subprocess import getstatusoutput


class StupidError(Exception): 
    # Constructor or Initializer 
    def __init__(self, value): 
        self.value = value 
  
    # __str__ is to print() the value 
    def __str__(self): 
        return(repr(self.value)) 


def clusterProfiler(gene_list, id_type='SYMBOL'):
    """Gene ontology analysis of a specific gene list, API of clusterProfiler (http://amp.pharm.mssm.edu/Enrichr/)
    
    Paramters
    ---------
    gene_list : list,
        list of genes for analysis
    id_type : str,
        type of input list, ENSEMBL / SYMBOL / ENTREZID

    Returns
    -------
    pd.DataFrame :
        enriched data
    """
    from pathlib import Path
    from utils import generate_random_key
    token = generate_random_key(6).upper()

    src_file = '/tmp/clusterProfiler_{}.txt'.format(token)
    tar_file = '/tmp/clusterProfiler_{}_GO.txt'.format(token)

    with open(src_file, 'w') as out:
        for gene_id in gene_list:
            out.write('{}\n'.format(gene_id))

    rcmd = Path(__file__).parent.parent / 'R' / 'clusterProfiler.R'

    # Run clusterProfiler
    status, output = getstatusoutput('/usr/bin/Rscript {} --input {} --output {} --type {}'.format(rcmd, src_file, tar_file, id_type))
    if status != 0:
        raise StupidError(output)

    # Get results
    enriched_data = pd.read_csv(tar_file, index_col=0)

    return enriched_data


def enrichR(gene_list, input_libraries=None):
    """Gene ontology analysis of a specific gene list, API of enrichR (http://amp.pharm.mssm.edu/Enrichr/)
    
    Paramters
    ---------
    gene_list : list,
        list of genes for analysis
    input_libraries : list,
        Specific analysis libraries

    Returns
    -------
    pd.DataFrame :
        enriched data
    """
    import json
    import requests

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(gene_list)
    description = 'auto enrichr list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    
    user_id = data['userListId']
    enrichr_libraries = ['GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018']

    enriched_data = []
    for library in enrichr_libraries:
        library_data = _gene_set_analysis(user_id, library)
        library_data['Group'] = library
        enriched_data.append(library_data)
    enriched_data = pd.concat(enriched_data)
    return enriched_data
        
    
def _gene_set_analysis(user_list_id, gene_set_library):
    import json
    import requests
    
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library)
     )
    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    
    enriched_data = pd.DataFrame(data[gene_set_library], 
                             columns=['Rank', 'Term', 'p', 'z', 'CombinedScore', 'Genes', 'q', 'old-p', 'old-q'])
    enriched_data.index = enriched_data['Rank']
    enriched_data = enriched_data.drop('Rank', axis=1)
    return enriched_data


def revcomp(seq):
    """
    Convert sequence to reverse complementary
    """
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]

