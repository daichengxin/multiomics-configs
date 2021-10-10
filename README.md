# Multiomics configs

Multiomics configs is a Github repository that contains information and configuration files for reanalysis of proteomics, genomics and transcriptomics datasets. The current repository including the following sections: 

## Proteomics annotated datasets

Annotated datasets are divided into multiple categories, mainly **absolute** and **differential** expression. 

### Absolute expression datasets

[Absolute expression datasets](https://github.com/multiomics/multiomics-configs/tree/master/datasets/absolute-expression) are datasets that aim to quantified absolute protein expression. Currently, datasets are devided into multiple categories: 

- cell-lines 
- platelet
- tissues 
- tumor

The **factor value** for each absolute expression dataset should be the **organism part** (tissue) in tissue datasets or the cell-line code in the cell line datasets. In absolute expression profiles each dataset contains a global SDRF where all the samples are annotated (e.g. [PXD000612](https://github.com/multiomics/multiomics-configs/tree/master/datasets/absolute-expression/cell-lines/PXD000612)). This representation allows to analyze the experiment complete or divided by samples. 

### Differential expression datasets 

[Differential expression datasets](https://github.com/multiomics/multiomics-configs/tree/master/datasets/differential-datasets) are datasets that aim to quantified the diffenretial expressed proteins in specific diseases or conditions. Currently, datasets are divided into multiple categories taking to account the analytical method: 

- label-free 
- dia 
- tmt 






 
