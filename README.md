# CHIMeRA (v. 03/31/2022)
## Usage
Explanation of the functions and how it is intended to be used:

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
```

## get_kegg_dblinks()
Similarly to the previous one, this function retrieves the information related to the voice "DBLINKS" of a single Kegg ID.
```python
get_kegg_dblinks('dr:D00001')
```
```
output: 
```

## kegg_dissect()
This function allows to slide the Kegg Ids stored inside a given dataframe -- in the form that we explained earlier -- from an initial position to an final one, both specified by the user. In addition, the user has to specify the function that will analyze each ID.
```python
kegg_dissect(kegg_df, 'ref', 0, 5)
```
```
output: 
```

## kegg_update()
This function simply uses the previous ```kegg_dissect()``` in a mulithread structure, in order to speed up the request process as much as possible. It's important to notice that we can't use many requests at the same times, otherwise the server will block the majority of them. For this reason it has been introduced a delay inside the function ```get_kegg_reference()``` and ```get_kegg_dblinks()```.
```python
kegg_update(kegg_df, 'dblinks')
```
```
output: 
```

## kegg_script()
This is simply a script that generates all the DataFrames that we need. It uses all the functions listed below.
```python
kegg_script(' ')
```
```
output: 
```
