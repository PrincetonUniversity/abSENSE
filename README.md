# abSENSE: a method to interpret undetected homologs

[![Try Now!][pyodide-badge]][pyodide-link]
[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[pyodide-badge]:            https://img.shields.io/badge/Pyodide-Try%20It!-blue.svg
[pyodide-link]:             https://potential-waffle-8b9309b7.pages.github.io/
[actions-badge]:            https://github.com/PrincetonUniversity/abSENSE/workflows/CI/badge.svg
[actions-link]:             https://github.com/PrincetonUniversity/abSENSE/actions
[pypi-link]:                https://pypi.org/project/abSENSE/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/abSENSE
[pypi-version]:             https://badge.fury.io/py/abSENSE.svg
[rtd-badge]:                https://readthedocs.org/projects/abSENSE/badge/?version=latest
[rtd-link]:                 https://abSENSE.readthedocs.io/en/latest/?badge=latest

## INTRODUCTION
abSENSE is a method that calculates the probability that a homolog of a given
gene would fail to be detected by a homology search (using BLAST or a similar
method) in a given species, even if the homolog were present and evolving
normally.

The result of this calculation informs how one interprets the result of a
homology search failing to find homologs of a gene in some species. One
possibility to explain such a result is that the gene is _actually absent_ from
the genome in that species: a biological, and potentially interesting (e.g. if
due to a gene loss or the birth of a new gene), result.

A second explanation, often ignored, is that the homolog _is_ present in the
genome of that species, but that the homology search merely lacks statistical
power to detect it. Here, the apparent absense of the homolog is a
technical/statistical limitation, and does not reflect underlying biology.

By calculating the probability that your homology search would fail to detect a
homolog _even if one were present_ and _even if it were evolving normally_
(e.g. no rate accelerations on a specific branch, potentially suggestive of
biologically interesting changes), abSENSE informs the interpretation of a
negative homology search result. If abSENSE finds that there is a high
probability of a homolog being undetected even if present, you may not be as
inclined to invoke a biological explanation for the result: the null model of a
failure of the homology search is sufficient to explain what you observe.

The method is explained in complete detail in [Weisman CM, Murray AW, Eddy SR
(2020) Many, but not all, lineage-specific genes can be explained by homology
detection failure. PLoS Biol 18(11):
e3000862.](https://doi.org/10.1371/journal.pbio.3000862)

There, it is applied to the specific case of lineage-specific genes, for which
homologs appear absent in all species outside of a narrow lineage. The method
itself is applicable to any case in which a homolog appears absent (e.g. a
single species missing a homolog that one might interpret as a gene loss), and
likewise, this code is applicable to all such cases.

This repository is a rewrite of the [original implementation](https://github.com/caraweisman/abSENSE/tree/c355c458e83722a0ffdf7284d4ea1f6f29ce205f)

## Installation

```
pip install absense
```
Will install all dependencies and make the command `absense` available to run.
Snakemake workflows to generate input files and perform validation can be accessed
through this repository in the `workflows` directory.

## Usage
### abSENSE THE BASICS

#### Try it!

The visualization method can be run in your browser, [here][pyodide-link].  You
can replace the sample data with your own or upload files to analyze.

#### Quickstart: main analysis

The analysis `absense` calculates the probabilities of homologs of genes from
some "focal" species being undetected in a set of `N` other species, for an
arbitrary number of genes. It requires two input files:

- `--bitscores`: A tab separated text file containing the bitscores of homologs
  of each gene to be analyzed in at least three of the species.

- `--distances`: A tab separated text file containing the N evolutionary
  distances, in substitutions/site, between the focal species and the other
  species. The distance between the focal species and itself should be 0. (If
  you don't already have such distances, the snakemake workflow can produce them
  from `.faa` files).

Examples of both of these files for a subset of genes from S. cerevisiae and
their orthologs 11 other fungal species can be found
[here](https://github.com/caraweisman/abSENSE/tree/c355c458e83722a0ffdf7284d4ea1f6f29ce205f/Quickstart_Examples)
They exemplify the formatting required for abSENSE to run (explained in more
detail below).

To run abSENSE on a given bitscore and distance file:

```
absense \
    --distances <distance_file> \
    --bitscores <score_file>
```

For example, to run abSENSE on the example fungal genes, type:

```
absense \
    --distances Quickstart_Examples/Fungi_Distances \
    --bitscores Quickstart_Examples/Fungi_Example_Bitscores
```

For each gene in the input bitscore file, the following will be computed and 
saved in a directory with the start time of the analysis:

- The probabilities of a homolog being undetected in each species, in
  `failure_probabilities.tsv`.

- The expected bitscores of homologs in each species in
  `predicted_bitscores.tsv`.

- The 99% confidence interval around this bitscore in each species, low and
  high bounds listed in separate files: `99PI_lower_prediction.tsv` and
  `99PI_high_prediction.tsv`.

- The parameters of each exponential fit in `parameters.tsv`.

You can specify the output directory name with the `--out-dir` option.
Additional information on output files is below.


#### Visualization

Two options allow `absense` to produce plots in the `plots` directory in
addition to the previously described files.  `--plot-all` produces a plot for
each gene in the bitscore file.  Alternatively, `--plot` takes a comma separated
list of genes to plot.  All genes are still analyzed.

To run abSENSE on gene GENEID contained in a given bitscore with a given
distance file:

```
absense \
    --distances <distance_file> \
    --bitscores <score_file> \
    --plot <GENEID>
```

For example, to analyze the S. cerevisiae gene Uli1, listed in the bitscore
file under its RefSeq ID (NP\_116682.3), type:

```
absense \
    --distances Quickstart_Examples/Fungi_Distances \
    --bitscores Quickstart_Examples/Fungi_Example_Bitscores \
    --plot NP_116682.3
```


### All Options

You can specify additional options, found by running `absense --help`.

- `--out-dir`: The directory to which your results will be output.
  Default is the time at which you ran the analysis.

- `--e-value`: The E-value threshold to be used (above this value, homologs will
  be considered undetected). Default is 0.001 (fairly permissive).

- `--gene-lengths`: Allows you to specify a file containing the lengths (in aa)
  of all genes in the bitscore file to be analyzed. Default is 400 amino acids
  (~average protein size in many species) for all proteins.

  abSENSE predicts a bitscore, which is then converted to an E-value to determine
  detectability; this conversation technically requires knowledge of both the
  size of the database in which the search occurs (see below) and the length of
  the gene being searched. Because the conversion between these values and
  E-value is logarithmic, though, only fairly large changes in these values
  substantially affect results.

- `--db-lengths`: Allows you to specify a file containing the sizes (in aa) of
  the database in which the homology search for each of your N species is
  performed. Default is 400 amino acids * 20,000 amino acids / gene = 8,000,000
  amino acids (~average protein and proteome size in many species) for all
  species.

- `--predict-all`: When set, Causes abSENSE to calculate the
  probability of homologs being undetected, the expected bitscores, and 99\%
  confidence intervals not only in species in which homologs were actually
  undetected, but also for species in which homologs have been found. This is
  obviously not the main use case, and is especially uninformative when those
  homologs and their bitscores have been used in the prediction itself (see
  below). May potentially be useful to see if a homolog in one species,
  although detected, seems to be behaving anomalously compared to those in
  other species (e.g. due to rate acceleration).

- `--include-only`: A comma separated list of species to consider. Allows you
  to restrict the species whose bitscores are used in predicting bitscores in
  other species. Mainly useful to do control-type analyses, such as Figure 5 in
  Weisman et al 2020, to show that predictions made from only a subset of the
  data are nonetheless reliable. If not specified, abSENSE uses all available
  bitscores in the prediction.

- `--plot-all`: When set, produce a plot of each gene's fit in the `plots`
  directory.

- `--plot`: A comma separated list of genes to plot and save in the `plots`
  directory.

- `--validate`: When set, will to tabulate how well the current distances
  fit the bitscores.


For example, to run an analysis on all S. cerevisiae proteins in the selected
fungal species in which the lengths of each S. cerevisiae protein and the sizes
of each species' search database are specified:

```
absense \
    --distances Fungi_Data/Fungi_Distances \
    --bitscores Fungi_Data/Fungi_Bitscores \
    --gene-lengths Fungi_Data/S_cer_Protein_Lengths \
    --db-lengths Fungi_Data/Fungi_Database_Lengths
```


### Output Files

`absense` produces six output files.

#### `failure_probabilities.tsv`

The central output of the program. For each gene in the analysis, this contains
the predicted probability that a homolog in each species would be undetected at
the specified E-value by a homology search, even if the homolog were present.

By default, this is only calculated in species in which the gene was not
detected. Results for species in which homologs were detected are therefore
listed as "detected". The setting `--predict-all` will calculate this value for
all species, even those in which a homolog was in fact detected.

If not enough data for a gene was provided to generate a bitscore prediction
(bitscores of homologs from at least three species are needed), the results
will read "not\_enough\_data".

#### `predicted_bitscores.tsv`

For each gene in the analysis, this contains the predicted (maximum likelihood)
bitscore of a homolog in each species.

By default, bitscores are only predicted in species in which the gene was not
detected. Results for species in which homologs were detected are therefore
listed as "detected". The setting `--predict-all` will calculates this value
for all species, even those for which a homolog was in fact detected. Here, the
known bitscore (often used in the prediction process; see the option
`--include-only`) will be shown alongside the prediction. If the known
bitscore was used in the fitting process, of course, these will usually be
quite similar!

If not enough data for a gene was provided to generate a bitscore prediction
(bitscores of homologs from at least three species are needed), the results
will read "not\_enough\_data".

#### `99PI_{lower,high}_prediction.tsv`

For each gene in the analysis, these contain the lower and upper bound of the 99\%
confidence interval for the bitscore of a homolog in each species.

#### `parameters.tsv`

For each gene in the analysis, this contains the best-fit (maximum likelihood)
values of the a and b parameters. (See Weisman et al 2020 for full
explanation.)

These a and b parameters are calculated from bitscores of homologs in species
included in the prediction process.

#### `run_info.txt`

Contains information about the analysis, including names of input files,
options/settings used, and the analysis time.

### Input File Formats

#### Required files:

##### The bitscore file

For an analysis of M genes in N species (including the focal species), the
bitscore file should be a tab-delimited file of N+1 columns by M+1 rows.  The
first row should begin with a blank entry, and should be followed by N entries
containing the names of the N species in your analysis. These names should
match those in the distance file (below) exactly.  The remaining M rows should
each begin with the name/identifier of the gene from the focal species to be
analyzed, followed by the bitscore of that gene against its homolog in the
species indicated at the top of the given column. For species in which homologs
are undetected, this value should be `0`. For species in which a homolog is
detected, but the orthology is unclear and so you wish to exclude it from being
used in the fit (see Weisman et al 2020), this value should be `N/A`.

#### The distance file

For an analysis with N species (including the focal species), the distance file
should be a tab-delimited file of 2 or more columns by N rows.  Entries in the first
column should contain the name of each species in your analysis. These names
should match those in the bitscore file (above) exactly.  Entries in the remaining
columns should contain the evolutionary distance between each species in the
indicated column and the focal species. (The distance between the focal species
and itself should always be 0.)  If additional distance estimates are provided,
the uncertainty in the distance will be considered during analysis.

#### Optional files:

##### The gene length file

For an analysis of M genes, the gene length file should be a tab-delimited file
of 2 columns by M rows.  The first column should contain the names/identifiers
of the gene from the focal species to be analyzed. These should match exactly
the names in the bitscore file (above).  The second column should contain the
length in amino acids of that gene.

If you don't already have such a file, here is a command to make one (named
`OUTPUTFILENAME`), from a FASTA file (`SEQFILENAME`) containing all of the
sequences: (It requires the easel package, which comes with the [HMMER software](
http://hmmer.org/documentation.html).)

```
esl-seqstat -a (SEQFILENAME) | \
    awk '{print $2 "\t" $3}' | \
    tac | \
    sed -e '1,7d' | \
    tac > (OUTFILENAME)
```


##### The database length file

For an analysis of N species, the database length file should be a
tab-delimited file of 2 columns by N rows.  The first column should contain the
names of the N species. These should match exactly the names in the bitscore
and distance files (above).  The second column should contain the sizes, in
amino acids, of the database in which the homology search for each species is
performed. For example, if you are searching for a gene in each of the species'
annotated proteomes, it should be the size of that species' proteome in amino
acids. If instead you are searching against a pan-specific database, like for
example NR, for all species, it should be the size of that database in amino
acids.

If you don't already have such a file, you can make one easily if you have the
databases themselves in eg FASTA format (again requiring the easel package,
downloadable with HMMER as in c) above: just run the command `esl-seqstat
(FASTA FILE)` on each database; this will report the total length in aa of
each database file. You can then put these into a tab-delimited file manually.

If your BLAST search is done on a nucleotide genome of the outgroup species via
TBLASTN, the database length that you should use for each genome is 2N, where N
is the genome length in nucleotides. (For a genome of N nucleotides, there are
~N/3 codons in it, which can be read in each of 6 possible reading frames, for
a total of 6N/3 = 2N amino acids.)
