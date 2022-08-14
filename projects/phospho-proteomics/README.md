## Phospho-proteomics reanalysis

The idea of this project is to reanalyze datasets to identified phospho-sites in public data. The following considerations are take into account to include projects: 

- TMT or LFQ experiments. 
- Human projects from Pride database


### Workflow 

Each dataset is annotated with an SDRF. Datasets are searched against human protein uniprot reviewed. The pipeline parameters for all experiments: 

- MSGF+ and Comet used
- PSM and protein FDR 1% filtered
- Luciphor LFR is computed for all the experiments where available for each msrun. **Note** It can happen for experiments with low number of phospho-sites per msrun that Luciphor can't be run, in those cases Luciophor is disabled. We are exploring the possibility of running AScore.


### Datasets

- PDC000326: CPTAC TMT (issues with protein inference never finished)
- PXD000759: 
- PXD001725:     
- PXD003531:    
- PXD005173: Done  
- PXD012255: 
- PXD014565: Error in the RAW files
- PDC000126: CPTAC TMT No annotated
- PDC000155: CPTAC TMT No annotated
- PXD000612:
- PXD001724:
- PXD001725: Run only with Comet
- PXD002255:
- PXD004452:
- PXD006482:
- PXD013923:
- PXD021865:
- PXD008012:
- PXD011175:
- PXD000497: 

