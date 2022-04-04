from kegg_requests import *
import concurrent.futures
import numpy as np

def kegg_update(kegg_df, info):
    '''
    Description
    -----------
    This function simply uses the previous function 'kegg_dissect()' in a mulithread structure, in order to speed up the
    request process as much as possible. It's important to notice that we can't use many requests at the same times,
    otherwise the server will block the majority of them. For this reason it has been introduced a delay inside the
    function 'get_kegg_reference()' and 'get_kegg_dblinks()'.

    Parameters
    ----------
    kegg_df : Pandas.DataFrame
        DataFrame that holds kegg IDs. The id column must be labeled as 'KeggID'.
    info : str
        Specifies the type of research to perform in kegg database. It can be 'ref' for 'REFERENCES' or 'dblinks' for
        'DBLINKS'.

    Returns
    -------
    output : list
        List of the results obtained by 'kegg_dissect()' function, which are the adjacency lists obtained from each
        Kegg ID given.

    '''
    n_elements = len(kegg_df.index)
    res, rest = divmod(n_elements, 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers = 4) as executor:
        f1 = executor.submit(kegg_dissect, kegg_df, info, 0, res)
        f2 = executor.submit(kegg_dissect, kegg_df, info, res, 2 * res)
        f3 = executor.submit(kegg_dissect, kegg_df, info, 2 * res, 3 * res)
        f4 = executor.submit(kegg_dissect, kegg_df, info, 3 * res, 4 * res + rest)
        result_futures = [f1, f2, f3, f4]
        results = [f.result() for f in concurrent.futures.as_completed(result_futures)]
        return [element for elements in results for element in elements]
    
    def kegg_script():
    '''
    Description
    -----------
    At the moment this script update the full database.
    '''
    kegg_path = "./Databases/Kegg/"
    kegg_references_path = "./Databases/Kegg_References/"
    kegg_genomes = pd.read_csv("https://rest.kegg.jp/list/genome", sep = "(?<=\d)\t", engine = "python", header = None,
                                        names=["KeggID", "Description"])
    kegg_diseases = pd.read_csv("https://rest.kegg.jp/list/disease", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    kegg_drugs = pd.read_csv("https://rest.kegg.jp/list/drug", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    kegg_enzymes = pd.read_csv("https://rest.kegg.jp/list/enzyme", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    
    kegg_genomes.to_csv(kegg_path+'kegg_genomes.csv')
    kegg_diseases.to_csv(kegg_path+'kegg_diseases.csv')
    kegg_drugs.to_csv(kegg_path+'kegg_drugs.csv')
    kegg_enzymes.to_csv(kegg_path + 'kegg_enzymes.csv')
    
    gen = kegg_update(kegg_genomes, 'ref')
    dis = kegg_update(kegg_diseases, 'ref')
    dru = kegg_update(kegg_drugs, 'ref')
    enz = kegg_update(kegg_enzymes, 'ref')

    adjacency_gen = ['dummy', 'dummy']
    adjacency_dis = ['dummy', 'dummy']
    adjacency_dru = ['dummy', 'dummy']
    adjacency_enz = ['dummy', 'dummy']
    for arr1, arr2, arr3, arr4 in zip(gen, dis, dru, enz):
        adjacency_gen = np.vstack([adjacency_gen, arr1])
        adjacency_dis = np.vstack([adjacency_dis, arr2])
        adjacency_dru = np.vstack([adjacency_dru, arr3])
        adjacency_enz = np.vstack([adjacency_enz, arr])

    adjacency_gen = np.delete(adjacency_gen, (0), axis = 0)
    adjacency_dis = np.delete(adjacency_dis, (0), axis = 0)
    adjacency_dru = np.delete(adjacency_dru, (0), axis = 0)
    adjacency_enz = np.delete(adjacency_enz, (0), axis = 0)
    
    pd.DataFrame(adjacency_gen, columns = ['KeggID', 'REFERENCE'])).to_csv(kegg_references_path+'kegg_genomes_ref.csv')
    pd.DataFrame(adjacency_dis, columns = ['KeggID', 'REFERENCE']).to_csv(kegg_references_path+'kegg_diseases_ref.csv')
    pd.DataFrame(adjacency_dru, columns = ['KeggID', 'REFERENCE']).to_csv(kegg_references_path+'kegg_drugs_refe.csv')
    pd.DataFrame(adjacency_enz, columns = ['KeggID', 'REFERENCE']).to_csv(kegg_references_path+'kegg_enzymes_ref.csv')
    
