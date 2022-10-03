# abSENSE Workflow

This directory contains a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
workflow for performing analysis with abSENSE starting from a public URL or
local faa file.  Simply modify the `params.yaml` file to setup the analysis
and run `snakemake` to produce the analyzed results.

## Parameters

- `workdir`: The working directory, stores all data and results in
  subdirectories.  Will be created if does not exist.
- `protdist_replicates`: To gain a confidence interval of the genetic distance
  between species, the prodist estimate is repeated.  All passing BUSCO genes
  will be included in the first estimate.  Subsequent analyses draw N genes
  from that population with replacement, N is the total number of genes.  This
  attempts to correct bias due to placement of neighboring genes when generating
  the protdist input file.
- `validate_genes`:  Determine how many genes in the focal species to include
  for validation.  Setting to 0 will disable the validation step.  Setting to
  -1 will instead run all genes from the focal species.
- `species`: A dictionary of species names and parameters for analysis.
  - `focal`: When set to `true`, the species is taken as the focal species.
    When not included, assumed to be false by default.
  - `url`: A link to download the faa file from the web.  When set, will download
    the faa file with wget
  - `accession`: A genbank or refseq accession number.  When set, will download
    the faa file with `datasets` from ncbi.
  - `file`: A file path to a local faa file.  When set, will soft link the file
    into the faa directory in the working directory.  Absolute paths are valid
    and relative paths are relative to the `workdir`.  File must be uncompressed.

  When multiple file paths are provided, the following order is used:
    1) local file
    2) download accession with datasets
    3) download url with wget

  Only one method is required.

## Paths

The `paths.yaml` determines the organization of the output file structure.
Do not change the keys (e.g. `faa:`) or the name of wildcards (e.g. `{species}`).
You can safely change anything else in the string, though keeping each file in
a separate directory can simplify the snakemake preprocessing.  For example,
if you wanted to store the fastas in a directory called `fasta_aa`, you can
change line 4 to
```
  faa: "fasta_aa/{species}.faa"
```
Note that indentation is important in yaml files.

## Other directories
The `envs` and `scripts` directories contain descriptions for conda environments
and custom analysis scripts, respectively.  You don't need to modify them, but
can view the `envs` content for version numbers of any dependencies
