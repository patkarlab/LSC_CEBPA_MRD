# LSC CEBPA MRD workflow description

This repository describes the workflow for analysing MRD samples sequenced using single gene amplicon based assay for CEBPA. Background error calculation as mentioned in [article](https://doi.org/10.1182/blood-2018-02-829911) can be found in error_model.md 

## Usage
The following parameters need to be modified in the `params` section of the `nextflow.config` before executing the workflow :
- *genome* = Complete path to the human genome fasta file(hg19_all.fasta). Please ensure that the BWA index files (hg19_all.fasta.fai, hg19_all.fasta.amb, hg19_all.fasta.ann, hg19_all.fasta.bwt, hg19_all.fasta.pac, hg19_all.fasta.sa) are also present in the same genome folder. The assests folder currently contains placeholder genome and index files.

- *annovar_humandb* = Complete path to the humandb database folder for [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/ "annovar site")

- *gen_ref*     = Complete path to the gene_fullxref.txt file as downloaded from the [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/ "annovar site") site

- *bedfile* = This file contains the coordinates of the amplicons used for this assay

- *outdir* = Location to write the output folder

## Running the pipeline
1. Transfer the sample input files `*.fastq.gz` inside the `sequences/` folder.

2. Modify the `samplesheet.csv`. The sample_ids, without the file extension, should be mentioned in samplesheet in the following format - <br>
sample1
sample2
sample3
Please check for empty lines in the samplesheet before running the pipeline.

3. To execute the pipeline, use the following command
```bash
nextflow run cebpa_mrd_amplicon.nf -entry CEBPA_MRD -profile docker -resume -bg
```

## Output
Samplewise output folders are written to the folder name mentioned in the `outdir` param in the config file.

