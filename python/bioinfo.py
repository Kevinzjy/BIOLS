def enrichment_analysis(gene_list, input_libraries=None):
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
    import pandas as pd

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
    import pandas as pd
    
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

