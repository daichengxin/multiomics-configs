# Differential expression analysis

- Differential expression datasets are Label-free, TMT or DIA datasets public in PRIDE.
- They can be run using multiple sdrf files depending on the variable under study (factor value).
- The analyses are run against the Human Swiss-Prot including isoforms (Homo-sapiens-uniprot-reviewed-isoforms-contaminants-decoy-202105.fasta). This is important because in other projects like the proteogenomics, we run against a bigger databases including variants etc. Using UNIPROT/Isoforms we aim to identified and quantified only canonical reviewed proteins.

- All the data is deposited in: https://ftp.pride.ebi.ac.uk/pride/data/proteomes/differential-expression/

### Differential Expression Label-Free

The datasets are run using the pipeline [proteomicsLFQ](https://github.com/nf-core/proteomicslfq).

- The parameters used for proteomicsLFQ are:

```bash
nextflow run /hps/nobackup2/proteomics/yperez_temp/proteomicslfq/main.nf -c /hps/nobackup2/proteomics/yperez_temp/proteomicslfq/nextflow.config -profile conda,lsf --root_folder **RAW_FILES** --database Homo-sapiens-uniprot-reviewed-isoforms-contaminants-decoy-202105.fasta --input SDRF --search_engines comet,msgf --protein_level_fdr_cutoff 0.01 --psm_pep_fdr_cutoff 0.05 --targeted_only false --outdir SDRF_OUTPUT --protein_inference bayesian --protein_quant shared_peptides --add_triqler_output -resume
```
