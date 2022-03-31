# CHIMeRA (v. 03/31/2022)
## Usage
Explanation of the functions and how they are intended to be used:

## Load Kegg Database into Pandas DataFrame
In order to load Kegg databases it is required to use a DataFrame that has that specific form, meaning two columns labeled as ```KeggID``` and ```Description```. The requests are based on [Kegg API](https://www.kegg.jp/kegg/rest/keggapi.html).
```python
kegg_df = pd.read_csv("https://rest.kegg.jp/list/genomes", sep = "(?<=\d)\t", engine = "python", header = None,
                  names=["KeggID", "Description"])
```

## get_kegg_reference()
This function searches in Kegg database the single Kegg Id given and it returns the information related to the voice "REFERENCE".
```python
get_kegg_reference('ds:H00047')
```
```
output: 
[['ds:H00047' 'PMID:15343276']
 ['ds:H00047' 'PMID:10526058']
 ['ds:H00047' 'PMID:11410326']]
```

## get_kegg_dblinks()
Similarly to the previous one, this function retrieves the information related to the voice "DBLINKS" of a single Kegg ID.
```python
get_kegg_dblinks('dr:D00001')
```
```
output:
[['dr:D00001' 'PubChem: 7847069']
 ['dr:D00001' 'ChEBI: 15377']
 ['dr:D00001' 'ChEMBL: CHEMBL1098659']
 ['dr:D00001' 'PDB-CCD: HOH O']
 ['dr:D00001' 'LigandBox: D00001']
 ['dr:D00001' 'NIKKAJI: J43.587B']
 ['dr:D00001' 'CAS: 7732-18-5']]
```

## kegg_dissect()
This function allows to slide the Kegg Ids stored inside a given dataframe -- in the form that we explained earlier -- from an initial position to an final one, both specified by the user. In addition, the user has to specify the function that will analyze each ID.
```python
kegg_dissect(kegg_df, 'ref', 0, 4)
```
```
output:
[array([['gn:T00001', 'PMID:7542800']], dtype='<U12'), 
array([['gn:T00002', 'PMID:7569993']], dtype='<U12'), 
array([['gn:T00003', 'PMID:8688087']], dtype='<U12'), 
array([['gn:T00004', 'PMID:8905231'],['gn:T00004', 'PMID:8590279']], dtype='<U12')]
```
The output of this function can be used to create a unique adjacency lists in this way:
```python
adjacency_gen = ['dummy', 'dummy']
for array in kegg_dissect(kegg_df, 'ref', 0, 4):
    adjacency_gen = np.vstack([adjacency_gen, array])
adjacency_gen = np.delete(adjacency_gen, (0), axis = 0)
```
```
output:
[['gn:T00001' 'PMID:7542800']
 ['gn:T00002' 'PMID:7569993']
 ['gn:T00003' 'PMID:8688087']
 ['gn:T00004' 'PMID:8905231']
 ['gn:T00004' 'PMID:8590279']]
```
## kegg_update()
This function simply uses the previous ```kegg_dissect()``` in a mulithread structure, in order to speed up the request process as much as possible. It's important to notice that we can't use many requests at the same times, otherwise the server will block the majority of them. For this reason it has been introduced a delay inside the function ```get_kegg_reference()``` and ```get_kegg_dblinks()```. In the example it is used only the first 4 entries (1 in 4 threads), but by default it uses all the entries in the dataframe.
```python
kegg_update(kegg_df, 'dblinks')
```
```
output:
       [array([['dr:D00004', 'PubChem: 7847072'],
       ['dr:D00004', 'ChEBI: 16526'],
       ['dr:D00004', 'ChEMBL: CHEMBL1231871'],
       ['dr:D00004', 'PDB-CCD: CO2'],
       ['dr:D00004', 'LigandBox: D00004'],
       ['dr:D00004', 'NIKKAJI: J43.600C'],
       ['dr:D00004', 'CAS: 124-38-9']], dtype='<U21'), 
       array([['dr:D00002', 'PubChem: 7847070'],
       ['dr:D00002', 'ChEBI: 15846'],
       ['dr:D00002', 'ChEMBL: CHEMBL1234613 CHEMBL1454168 CHEMBL3039307'],
       ['dr:D00002', 'DrugBank: DB01907'],
       ['dr:D00002', 'PDB-CCD: NAD NAJ'],
       ['dr:D00002', 'LigandBox: D00002'],
       ['dr:D00002', 'NIKKAJI: J136.554A'],
       ['dr:D00002', 'CAS: 53-84-9']], dtype='<U49'), 
       array([['dr:D00001', 'PubChem: 7847069'],
       ['dr:D00001', 'ChEBI: 15377'],
       ['dr:D00001', 'ChEMBL: CHEMBL1098659'],
       ['dr:D00001', 'PDB-CCD: HOH O'],
       ['dr:D00001', 'LigandBox: D00001'],
       ['dr:D00001', 'NIKKAJI: J43.587B'],
       ['dr:D00001', 'CAS: 7732-18-5']], dtype='<U21'), 
       array([['dr:D00003', 'PubChem: 7847071'],
       ['dr:D00003', 'ChEBI: 15379'],
       ['dr:D00003', 'ChEMBL: CHEMBL1234886'],
       ['dr:D00003', 'PDB-CCD: OXY'],
       ['dr:D00003', 'LigandBox: D00003'],
       ['dr:D00003', 'NIKKAJI: J44.420K'],
       ['dr:D00003', 'CAS: 7782-44-7']], dtype='<U21')]

```

## kegg_script()
This is simply a script that generates all the DataFrames that we need. It uses all the functions listed below. To perform a script that recovered the references from the complete ```genome``` and ```disease``` Kegg databases (10870 total kegg ids) it took ```time ~ 1h 17min```, using a delay of ```0.9 s``` and a four thread computation.
```python
kegg_script()
```
