# URL update
**The PDC data download client will reject manifests older than one week (168 hours), as the download URLs will no longer be valid.**

So, we should update PDC file.

1. Place re-download PDC_file_manifest file from PDC database.
<center>
    <img src=".\example.png" width="400" height="300">
</center>

2. Set the file path.
```python
PDC_path = r"D:\tumor\PDC000436\PDC000436.sdrf.tsv"
file_manifest_path = r"D:\tumor\PDC000436\PDC_file_manifest_06032023_182700.csv"
```

3. Update URL
```python
def update_url(PDC_path,file_manifest_path):
    PDC = pd.read_csv(PDC_path,sep='\t')
    file_manifest = pd.read_csv(file_manifest_path)
    file_manifest = file_manifest[['File Name','File Download Link']]
    file_manifest.columns = ['comment[data file]','url']
    PDC = PDC.merge(file_manifest,on="comment[data file]",how='left')
    PDC.loc[:,'comment[file uri]'] = PDC['url']
    PDC.drop(columns='url',inplace=True)
    colname = PDC.columns.tolist()
    colname = [coln.split('.')[0] for coln in colname]
    PDC.to_csv(PDC_path,index=False,header=colname,sep='\t')

update_url(PDC_path,file_manifest_path)
```

**Pace refer update_url.ipynb**