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

The **factor value** for each absolute expression dataset should be the **organism part** (tissue) in tissue datasets or the cell-line code in the cell line datasets. In absolute expression profiles each dataset contains a global SDRF where all the samples are annotated (e.g. [PXD000612](https://github.com/multiomics/multiomics-configs/tree/master/datasets/absolute-expression/cell-lines/PXD000612)). This representation allows analyzing the experiment complete or divided by samples.

### Differential expression datasets

[Differential expression datasets](https://github.com/multiomics/multiomics-configs/tree/master/datasets/differential-datasets) are datasets that aim to quantify the differential expressed proteins in specific diseases or conditions. Currently, datasets are divided into multiple categories taking to account the analytical method:

- label-free
- dia
- tmt

**NOTE**: For each dataset only one factor value should be added for each SDRF. If more than one variable is studied, then multiple SDRFs should be added with the following structure PXD-{factor value}.sdrf.tsv

If possible, creates the IDF for each dataset.

## Databases

The Database folder contains multiple databases created for the reanalysis of data including UniProt, ENSEMBL, contaminants databases.


## Projects

### Detecting Non-canonical Peptides in Cell-lines and Tumor data

LFQ and TMT datasets are searched against non-canonical and tissue-specific variant databases generated with pgdb nextflow pipeline. The list of datasets can be found [here](https://github.com/multiomics/multiomics-configs/tree/master/projects/non-canonical)

## Citing

Please if you use this repo, cite the following manuscripts: 

- Dai C, Pfeuffer J, Wang H, Zheng P, Käll L, Sachsenberg T, Demichev V, Bai M, Kohlbacher O, Perez-Riverol Y. quantms: a cloud-based pipeline for quantitative proteomics enables the reanalysis of public proteomics data. Nature Methods. 2024 Sep;21(9):1603-7.
- Dai C, Füllgrabe A, Pfeuffer J, Solovyeva EM, Deng J, Moreno P, Kamatchinathan S, Kundu DJ, George N, Fexova S, Grüning B. A proteomics sample metadata representation for multiomics integration and big data analysis. Nature Communications. 2021 Oct 6;12(1):5854.

### Contributing

If you want to add a dataset to the repository, please create a Pull request with the annotations.





