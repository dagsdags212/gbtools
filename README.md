# gbtools

A simple interface for downloading and parsing Genbank records from NCBI.

## Installation

Clone the git repository:
```
git clone https://github.com/dagsdags212/gbtools.git
cd gbtools
```

Install using pip:
```
pip install .
```

## Dependencies

```
requires-python = ">=3.10"
dependencies = [
    "biopython>=1.85",
    "click>=8.1.8",
    "ipython>=8.34.0",
    "marimo>=0.11.26",
    "polars>=1.26.0",
    "rich>=13.9.4",
    "setuptools>=78.0.2",
]
```

## Usage

### Downloading Records

Downlad a record from Genbank using an accession:
```bash
gbtools fetch AB050936.1
```

Specify an output file path:
```bash
gbtools fetch AB050936.1 -o AB050936.gbk
```

Or save it using output redirection:
```bash
gbtools fetch AB050936.1 > AB050936.gbk
```

### Retrieving Metadata

Get an overview of the contents of a Genbank record:
```bash
gbtools metadata AB050936.1
```

Get a list of annotated genes from a Genbank file:
```bash
# Download the Genbank record.
gbtools fetch AB050936.1 > AB050936.gbk

# Extract the gene features.
gbtools features AB050936.gbk

# Extract only coding regions.
gbtools features AB050936.gbk --cds

# Using tabs as separators.
gbtools features AB050936.gbk --cds --tsv
```

Get features without saving an intermediate gbk file:
```bash
gbtools fetch AB050936.1 | gbtools features > AB050936_features.csv
```

### Extract full sequences

Extract the full sequence in FASTA format:
```bash
gbtools seq AB050936.1
```

Pass in a Genbank filepath to avoid making an API request to the Entrez server:
```bash
# Save the sequence using the -o flag.
gbtools seq AB050936.gbk -o AB050936.fa

# Redirect output to a new file.
gbtools seq AB050936.gbk > AB050936.fa
```

### Extract sequences between annotated regions

In some cases, you would like to extract non-coding regions. This can be done for annotated features of a Genbank record using the following pattern:
```bash
gbtools between <gbk> <region_1> <region_2>
```

As an example, we could extract the intergenic region between _VP35_ and _VP40_ of `AB050936.gbk` by invoking:
```bash
gbtools between AB050936.gbk VP35 VP40
```

The order of the two regions should not matter. Specifying _VP40_ before _VP35_ would return the same sequence.

Pass in the `--translate` flag to return an amino acid sequence:
```bash
gbtools between AB050936.gbk VP35 VP40 --translate
```
