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
    function 'get_kegg_pmid()' and 'get_kegg_dblinks()'.

    Parameters
    ----------
    kegg_df : Pandas.DataFrame
        DataFrame that holds kegg IDs. The id column must be labeled as 'KeggID'.
    info : str
        Specifies the type of research to perform in kegg database. It can be 'pmid' for 'REFERENCES' or 'dblinks' for
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
    kegg_path = "./Databases/Kegg/"
    kegg_references_path = "./Databases/Kegg_PMID/"
    kegg_dblinks_path = "./Databases/Kegg_dblinks/"
    adjlists_path = "./Databases/ConnectionTables/"

    #Database used: ds, ec, dr, gn, cpd, genes(vg, vp, ag)
    #No Links: dr-ec, dr-gn, dr-dr, ds-ec, ds-ds, ec-ec, ec-gn, gn-gn, any genes

    print('\n\u2022 Loading Kegg self adjacency lists ... ', end='')

    kegg_dr_ds = pd.read_csv("https://rest.kegg.jp/link/dr/ds", sep = "(?<=\d)\t", engine = "python", header = None,
                                        names=["KeggID_1", "KeggID_2"])
    print('dr-ds '+ u'\u2713, ', end = '')

    kegg_dr_cpd = pd.read_csv("https://rest.kegg.jp/link/dr/cpd", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID_1", "KeggID_2"])
    print('dr-cpd '+ u'\u2713, ', end = '')

    kegg_ds_cpd = pd.read_csv("https://rest.kegg.jp/link/ds/cpd", sep = "(?<=\d)\t", engine = "python", header = None,
                              names = ["KeggID_1", "KeggID_2"])
    print('ds-cpd '+ u'\u2713, ', end = '')

    kegg_ds_gn = pd.read_csv("https://rest.kegg.jp/link/ds/gn", sep = "\t", engine = "python", header = None,
                              names = ["KeggID_1", "KeggID_2"])
    print('ds-gn '+ u'\u2713, ', end = '')

    kegg_ec_cpd = pd.read_csv("https://rest.kegg.jp/link/ec/cpd", sep = "(?<=\d)\t", engine = "python", header = None,
                              names = ["KeggID_1", "KeggID_2"])
    print('ec-cpd '+ u'\u2713, ', end = '')

    kegg_gn_cpd = pd.read_csv("https://rest.kegg.jp/link/gn/cpd", sep = "(?<=\d)\t", engine = "python", header = None,
                              names = ["KeggID_1", "KeggID_2"])
    print('gn-cpd '+ u'\u2713.', end = '')

    print('\n\u2022 Saving adjacency lists in \'' + adjlists_path + '\' ... ', end = '')
    kegg_dr_ds.to_csv(adjlists_path + 'kegg_dr_ds_adjlist.csv')
    kegg_dr_cpd.to_csv(adjlists_path + 'kegg_dr_cpd_adjlist.csv')
    kegg_ds_cpd.to_csv(adjlists_path + 'kegg_ds_cpd_adjlist.csv')
    kegg_ds_gn.to_csv(adjlists_path + 'kegg_ds_gn_adjlist.csv')
    kegg_ec_cpd.to_csv(adjlists_path + 'kegg_ec_cpd_adjlist.csv')
    kegg_gn_cpd.to_csv(adjlists_path + 'kegg_gn_cpd_adjlist.csv')
    print('Completed.', end = '')

    print('\n\n\u2022 Loading Kegg databases ... ', end = '')

    kegg_genomes = pd.read_csv("https://rest.kegg.jp/list/genome", sep = "(?<=\d)\t", engine = "python", header = None,
                                        names=["KeggID", "Description"])
    print('genomes ' + u'\u2713, ', end = '')

    kegg_diseases = pd.read_csv("https://rest.kegg.jp/list/disease", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    print('diseases ' + u'\u2713, ', end = '')

    kegg_drugs = pd.read_csv("https://rest.kegg.jp/list/drug", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    print('drugs ' + u'\u2713, ', end = '')

    kegg_enzymes = pd.read_csv("https://rest.kegg.jp/list/enzyme", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    print('enzymes ' + u'\u2713, ', end = '')

    kegg_compounds = pd.read_csv("https://rest.kegg.jp/list/compound", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    print('compounds ' + u'\u2713, ', end = '')

    kegg_genes_ag = pd.read_csv("https://rest.kegg.jp/list/ag", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    kegg_genes_vp = pd.read_csv("https://rest.kegg.jp/list/vp", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])
    kegg_genes_vg = pd.read_csv("https://rest.kegg.jp/list/vg", sep = "(?<=\d)\t", engine = "python", header = None,
                               names = ["KeggID", "Description"])

    kegg_genes = pd.concat([kegg_genes_ag, kegg_genes_vp, kegg_genes_vg], ignore_index = True)
    print('genes ' + u'\u2713 ', end = '')

    print('... Completed.', end = '')

    print('\n\u2022 Saving Kegg databases in \'' + kegg_path + '\' ... ', end = '')
    kegg_genomes.to_csv(kegg_path+'kegg_genomes.csv', index = False)
    kegg_diseases.to_csv(kegg_path+'kegg_diseases.csv', index = False)
    kegg_drugs.to_csv(kegg_path+'kegg_drugs.csv', index = False)
    kegg_enzymes.to_csv(kegg_path + 'kegg_enzymes.csv', index = False)
    kegg_compounds.to_csv(kegg_path + 'kegg_compounds.csv', index = False)
    kegg_genes.to_csv(kegg_path + 'kegg_genes.csv', index = False)
    print('Completed.', end = '')

    print('\n\n\u2022 Get Kegg-PMID adjacency lists ... ', end = '')

    gen = kegg_update(kegg_genomes, 'pmid')
    print('genomes-PMID ' + u'\u2713, ', end = '')

    dis = kegg_update(kegg_diseases, 'pmid')
    print('diseases-PMID ' + u'\u2713, ', end = '')

    dru = kegg_update(kegg_drugs, 'pmid')
    print('drugs-PMID ' + u'\u2713, ', end = '')

    enz = kegg_update(kegg_enzymes, 'pmid')
    print('enzymes-PMID ' + u'\u2713, ', end = '')

    gns = kegg_update(kegg_genes, 'pmid')
    print('genes-PMID ' + u'\u2713 ', end = '')

    print('... Completed.', end = '')

    print('\n\n\u2022 Get Kegg dblinks ... ', end = '')

    cpd = kegg_update(kegg_compounds, 'dblinks')
    print('compounds-dblinks ' + u'\u2713 ', end = '')

    print('... Completed.', end = '')

    adjacency_gen = ['dummy', 'dummy']
    adjacency_dis = ['dummy', 'dummy']
    adjacency_dru = ['dummy', 'dummy']
    adjacency_enz = ['dummy', 'dummy']
    adjacency_cpd = ['dummy', 'dummy']
    adjacency_gns = ['dummy', 'dummy']

    for arr1, arr2, arr3, arr4, arr5, arr6 in zip(gen, dis, dru, enz, cpd, gns):
        adjacency_gen = np.vstack([adjacency_gen, arr1])
        adjacency_dis = np.vstack([adjacency_dis, arr2])
        adjacency_dru = np.vstack([adjacency_dru, arr3])
        adjacency_enz = np.vstack([adjacency_enz, arr4])
        adjacency_cpd = np.vstack([adjacency_enz, arr5])
        adjacency_gns = np.vstack([adjacency_gns, arr6])

    print('\n\n\u2022 Saving Kegg-PMID adjacency lists in \'' + kegg_references_path + '\' ... ', end = '')

    pd.DataFrame(adjacency_gen, columns = ['KeggID', 'PMID']).to_csv(kegg_references_path+'kegg_genomes_ref.csv')
    pd.DataFrame(adjacency_dis, columns = ['KeggID', 'PMID']).to_csv(kegg_references_path+'kegg_diseases_ref.csv')
    pd.DataFrame(adjacency_dru, columns = ['KeggID', 'PMID']).to_csv(kegg_references_path+'kegg_drugs_ref.csv')
    pd.DataFrame(adjacency_enz, columns = ['KeggID', 'PMID']).to_csv(kegg_references_path+'kegg_enzymes_ref.csv')
    pd.DataFrame(adjacency_gns, columns = ['KeggID', 'PMID']).to_csv(kegg_references_path + 'kegg_genes_ref.csv')

    print('Completed.', end = '')

    print('\n\n\u2022 Saving Kegg-dblinks in \'' + kegg_dblinks_path + '\' ... ', end = '')

    pd.DataFrame(adjacency_cpd, columns = ['KeggID', 'DBLINKS']).to_csv(kegg_dblinks_path + 'kegg_compounds_dblinks.csv')
    
    print('Completed.', end = '')
