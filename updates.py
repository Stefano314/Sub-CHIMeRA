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