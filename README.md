# ChromiumSingleCellPipeline
This repository contains scripts that can be used to produce files to operate
the 10x Genomics Cell Ranger pipeline.

## Usage
### Reference preparation
References can be prepared using the instructions as detailed on the [10x Genomics
website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references).

This repository includes [nucFasta2GTF.py](nucFasta2GTF.py), which creates a GTF
file based on FASTA sequences. This script makes two assumptions; first, the
FASTA sequence represents the nucleotide sequence of a single feature, and two,
the name of the gene is the same name used in the header of the FASTA sequence.
The script will work on single and multiple FASTA sequences in a file.

Usage:

```bash
python2.7 nucFasta2GTF.py -i input_fasta.fa -o output_gtf.gtf
```

### Generating PBS scripts
PBS scripts can be generated using the [scRNASeqPipe.py](scRNAseqPipe.py) script.

Usage:

```bash
python2.7 scRNAseqPipe.py --project ExampleProjectSheet.yaml --config pipeline_config.yaml
```

This script will generate PBS files for each stage of the pipeline. You may
specify which stages to run in the project YAML file.

If you have any questions, feel free to email me at a.senabouth@garvan.org.au.
