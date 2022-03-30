from kegg_requests import *
import concurrent.futures
import numpy as np

def kegg_update(kegg_df, info):
    n_elements = len(kegg_df.index)
    res, rest = divmod(n_elements, 4)
    res, rest = 5,0
    with concurrent.futures.ThreadPoolExecutor(max_workers = 4) as executor:
        f1 = executor.submit(kegg_dissect, kegg_df, info, 0, res)
        f2 = executor.submit(kegg_dissect, kegg_df, info, res, 2 * res)
        f3 = executor.submit(kegg_dissect, kegg_df, info, 2 * res, 3 * res)
        f4 = executor.submit(kegg_dissect, kegg_df, info, 3 * res, 4 * res + rest)
        result_futures = [f1, f2, f3, f4]
        results = [f.result() for f in concurrent.futures.as_completed(result_futures)]
        return [element for elements in results for element in elements]
    
    def kegg_script():
    db_path = "./Databases/Kegg"
    kegg_genomes = pd.read_csv("https://rest.kegg.jp/list/genome", sep = "(?<=\d)\t", engine = "python", header = None,
                                        names=["KeggID", "Description"])
    kegg_diseases = pd.read_csv("https://rest.kegg.jp/list/disease", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    kegg_drugs = pd.read_csv("https://rest.kegg.jp/list/drug", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    # kegg_genomes.to_csv(db_path+'kegg_genomes.csv')
    # kegg_diseases.to_csv(db_path+'kegg_diseases.csv')
    # kegg_drugs.to_csv(db_path+'kegg_drugs.csv')
    gen = kegg_update(kegg_genomes, 'ref')
    dis = kegg_update(kegg_diseases, 'ref')
    dru = kegg_update(kegg_drugs, 'ref')

    adjacency_gen = ['dummy', 'dummy']
    adjacency_dis = ['dummy', 'dummy']
    adjacency_dru = ['dummy', 'dummy']
    for array1, array2, array3 in zip(gen, dis, dru):
        adjacency_gen = np.vstack([adjacency_gen, array1])
        adjacency_dis = np.vstack([adjacency_dis, array2])
        adjacency_dru = np.vstack([adjacency_dru, array3])

    adjacency_gen = np.delete(adjacency_gen, (0), axis = 0)
    adjacency_dis = np.delete(adjacency_dis, (0), axis = 0)
    adjacency_dru = np.delete(adjacency_dru, (0), axis = 0)
    
