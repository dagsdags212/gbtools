import json
import os
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.pretty import pprint
from rich.table import Table
from typing_extensions import Annotated

from gbtools.utils import *

app = typer.Typer(no_args_is_help=True)

ECONFIG = EntrezConfig(
    email=os.environ.get("EMAIL", "jegsamson.dev@gmail.com"),
    api_key=os.environ.get("NCBI_API_KEY"),
)


@app.command()
def fetch(
    id: Annotated[str, typer.Argument()],
    outfile: Annotated[
        Path, typer.Option("--file", help="filepath for storing the Genbank record")
    ] = Path("record.gbk"),
    verbose: Annotated[
        bool, typer.Option("--verbose", help="increase verbosity of program")
    ] = False,
):
    """Retrieve a Genbank record from NCBI"""
    if outfile == "":
        handle = fetch_genbank_handle(id=id, config=ECONFIG, verbose=verbose)
        print(handle.read())
    else:
        outfile = Path(outfile)
        if not str(outfile).endswith(".gbk") and not str(outfile).endswith(".gb"):
            outfile = Path(str(outfile) + ".gbk")
        fetch_genbank_record(id=id, outfile=outfile, config=ECONFIG, verbose=verbose)


@app.command()
def metadata(
    id: Annotated[str, typer.Argument(default=...)],
    verbose: Annotated[
        bool, typer.Option("--verbose", help="increase verbosity of program")
    ] = False,
):
    """Display metadata from a Genbank accession"""
    record = fetch_genbank_record(id=id, config=ECONFIG, verbose=verbose)
    metadata = vars(record)
    metadata.pop("features")
    pprint(metadata, expand_all=True)


@app.command()
def extract(
    file: Annotated[str, typer.Option("--file", "-f", help="path to a Genbank file")],
    cds: Annotated[
        bool, typer.Option("--cds", "-p", help="extract protein sequences")
    ] = False,
    header: Annotated[
        bool, typer.Option("--header", "-H", help="include header line")
    ] = False,
    tsv: Annotated[
        bool, typer.Option("--tsv", "-t", help="use tabs for separating field columns")
    ] = False,
    verbose: Annotated[
        bool, typer.Option("--verbose", help="increase verbosity of program")
    ] = False,
):
    """Get all sequence features from a Genbank record"""
    gbk = load_gbk(file)
    if cds:
        features = gbk2features(gbk, "cds", verbose=verbose)
    else:
        features = gbk2features(gbk, "gene", verbose=verbose)

    sep = "\t" if tsv else ","

    if header:
        print(f"name{sep}pos_start{sep}pos_end{sep}strand")

    for name, info in features.items():
        print(
            f"{name}{sep}{info['pos_start']}{sep}{info['pos_end']}{sep}{info['strand']}"
        )


@app.command()
def fasta(
    id: Annotated[Optional[str], typer.Argument()] = None,
    file: Annotated[
        Optional[str], typer.Option("--file", "-f", help="path to a Genbank file")
    ] = None,
    outfile: Annotated[
        str,
        typer.Option("--output", "-o", help="filepath for storing the FASTA record"),
    ] = "",
    verbose: Annotated[
        bool, typer.Option("--verbose", help="increase verbosity of program")
    ] = False,
):
    """Extract the sequence from Genbank record in FASTA format"""
    if file:
        if outfile:
            gbk2fasta(infile=file, outfile=outfile, verbose=verbose, config=ECONFIG)
        else:
            fa_record = gbk2fasta(infile=file, verbose=verbose, config=ECONFIG)
            print_fasta(fa_record)

    if id:
        if outfile:
            gbk2fasta(id=id, outfile=outfile, verbose=verbose, config=ECONFIG)
        else:
            fa_record = gbk2fasta(id=id, verbose=verbose)
            print_fasta(fa_record)


@app.command()
def config(
    email: Annotated[
        str, typer.Option("--email", "-e", help="add your email address")
    ] = "",
    api_key: Annotated[
        str, typer.Option("--api-key", "-k", help="add an NCBI API key")
    ] = "",
    list: Annotated[
        bool, typer.Option("--list", "-l", help="display current configuration")
    ] = False,
) -> None:
    """Provide credentials (email, API key) to Entrez server"""
    curr_config = {}
    config_file = Path("/tmp/gbtools_auth.json")
    if config_file.exists() and config_file.is_file():
        curr_config = load_config().to_json()

    if email:
        curr_config.update({"EMAIL": email})

    if api_key:
        curr_config.update({"NCBI_API_KEY": api_key})

    with config_file.open("w") as fh:
        json.dump(curr_config, fh)

    if list:
        if len(curr_config) == 0:
            print("Config file not initialized")
            return
        for k, v in curr_config.items():
            print(f"{k}={v}")


if __name__ == "__main__":
    app()
