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

## Usage

### Downloading Records

Downlad a record from Genbank using an accession:
```bash
gbtools fetch AB050936.1
```

Specify a file path:
```bash
gbtools fetch AB050936.1 --out gb/AB050936.gbk
```

### Retrieving Metadata

Get an overview of the contents of a Genbank record:
```bash
gbtools info AB050936.1
```

Get a list of annotated genes in tabular format:
```bash
gbtools list-genes AB050936.1
```

Get a list of annotated genes in a tree hierarchy:
```bash
gbtools list-genes AB050936.1 --tree
```
