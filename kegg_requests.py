import numpy as np
import pandas as pd
import time as tm

def get_kegg_reference(keggid):
    """
    Description
    -----------
    Search in Kegg database the KeggID given and returns -- if they exist -- the associated REFERENCES in the form of
    an adjacency list.

    Parameters
    ----------
    keggid : str
        KeggID that will be searched (ex. 'gn:T00004').

    Note
    ----
    - The function uses a delay of 0.9 s to avoid server rejections to multiple requests in multithread.

    Returns
    -------
    output : ndarray
        Numpy array containing the adjacency list between the keggid given and the references found.
    """
    tm.sleep(0.9)
    df = pd.read_csv('http://rest.kegg.jp/get/' + keggid, names = ['filelines'], sep = "NoSeparatorRequired", engine = 'python')
    ref_positions = df['filelines'].str.startswith('REFERENCE').values
    if np.count_nonzero(ref_positions)>0:
        ref_list = df[ref_positions].values.tolist()
        pmid_list = [[i for i in [a.replace('[','').replace(']','') for a in pmid.split()]
                        if i.startswith('PMID')] for row in ref_list for pmid in row]
        pmid_list = [x for x in pmid_list if x != []]
        k_id = np.array([keggid] * len(pmid_list))
        if k_id.size == 0:
            return np.array([keggid, None])
        return np.array([k_id, np.array([item for sublist in pmid_list for item in sublist])]).T
    else:
        return np.array([keggid, None])

def get_kegg_dblinks(keggid, stopwords = ['REFERENCE','ATOM','BOND','///']):
    """
    Description
    -----------
    Search in Kegg database the KeggID given and returns -- if they exist -- the associated DBLINKS in the form of
    an adjacency list.

    Parameters
    ----------
    keggid : str
        KeggID that will be searched (ex. 'gn:T00004').
    stopwords : list, optional.
        Contains a list of possible words that are after the section of DBLINKS in Kegg database.

    Note
    ----
    - The function uses a delay of 0.9 s to avoid server rejections to multiple requests in multithread.
    - The parameter 'stopwords' is required to understand where is the possible end of the DBLINKS in Kegg server, which is usually different from database to database.
    - The default values of 'stopwords' are given with the following criterion:
        genome : REFERENCE
        drug : ATOM
        disease : REFERENCE
        enzyme : ///
        compound : ATOM, BOND

    Returns
    -------
    output : ndarray
        Numpy array containing the adjacency list between the keggid given and the dblinks found.
    """
    tm.sleep(0.9)
    df = pd.read_csv('http://rest.kegg.jp/get/' + keggid, names = ['filelines'], sep = "NoSeparatorRequired", engine = 'python')
    dblinks_bool = df['filelines'].str.startswith('DBLINKS').values
    end_indexes = []
    if np.any(dblinks_bool == True):
        beg = np.where(dblinks_bool == True)[0][0]
        for word in stopwords:
            if np.where(df['filelines'].str.startswith(word))[0].size != 0 and \
                    np.where(df['filelines'].str.startswith(word)==True)[0][0] > beg:
                end_indexes.append(np.where(df['filelines'].str.startswith(word).values[beg + 1:] == True)[0][0])
        end = np.min(end_indexes)
        res = list(df['filelines'].loc[beg+1:beg+end].values)
        res.append(' '.join(df['filelines'].loc[beg].split()[1:]))
        res = np.array(res)
        k_id = np.array([keggid] * len(res))
        return np.array([k_id, res]).T
    else:
        return np.array([keggid, None])

def kegg_dissect(df, info, beg, end):
    """
    Description
    -----------
    This function allows to slide in 'df' dataframe from 'beg' to 'end' and it performs for every object found the
    function specified by 'info'.

    Parameters
    ----------
    df : Pandas.DataFrame
        DataFrame that holds kegg IDs. The id column must be labeled as 'KeggID'.
    info : str
        Specifies the type of research to perform in kegg database. It can be 'ref' for 'REFERENCES' or 'dblinks' for
        'DBLINKS'.
    beg : int
        Start position in the dataframe.
    end : int
        End position in the dataframe.

    Note
    ----
    - The function is used in 'kegg_update', since it fits well a multithread computation.

    Returns
    -------
    output : list
        List of the adjacency lists obtained from each Kegg ID given.

    Examples
    --------
    >>> adjacency_list = ['dummy','dummy']
    >>> result = kegg_dissect(kegg_df, 'ref', 0, 5)
    >>> for array in result:
    >>>     adjacency_list = np.vstack([adjacency_list, array])
    >>> adjacency_list = np.delete(adjacency_list, (0), axis=0)
        [['gn:T00001' 'PMID:7542800']
         ['gn:T00002' 'PMID:7569993']
         ['gn:T00003' 'PMID:8688087']
         ['gn:T00004' 'PMID:8905231']
         ['gn:T00004' 'PMID:8590279']
         ['gn:T00005' 'PMID:8849441']]
    """
    if info == 'ref':
        return [get_kegg_reference(x) for x in df['KeggID'].iloc[beg:end].values]
    elif info == 'dblinks':
        return [get_kegg_dblinks(x) for x in df['KeggID'].iloc[beg:end].values]
    
