import json
from dataclasses import dataclass
from io import StringIO, TextIOWrapper
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass
class EntrezConfig:
    """Provides the Entrez server with user-related data

    Calling an instance of this class would provide Biopython with user metadata
    that restricts the rate-limit of the user depending if an API key is supplied.
    """

    email: Optional[str]
    api_key: Optional[str] = None
    max_tries: int = 3
    sleep_between_tries: int = 5

    def __post_init__(self) -> None:
        self.path: Path = Path('/tmp/gbtools_auth.json')
        if self.api_key:
            self.allowed_queries_per_second = 10
        else:
            self.allowed_queries_per_second = 3

    def __call__(self) -> None:
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        Entrez.max_tries = self.max_tries
        Entrez.sleep_between_tries = self.sleep_between_tries

    def save(self) -> None:
        """Write attributes to a config file."""
        data = dict(
            EMAIL=self.email,
            NCBI_API_KEY=self.api_key,
            max_tries=self.max_tries,
            sleep_between_tries=self.sleep_between_tries,
        )
        with self.path.open("w") as fh:
            json.dump(data, fh)

    def to_json(self) -> dict[str, Optional[str | int]]:
        """Return config file as a dictionary."""
        return dict(
            EMAIL=self.email,
            NCBI_API_KEY=self.api_key,
            max_tries=self.max_tries,
            sleep_between_tries=self.sleep_between_tries,
        )


def load_config() -> EntrezConfig:
    """Loads JSON data from `/tmp/gbtools_auth.json` as an EntrezConfig object.

    Parameters
    ----------
    None

    Returns
    -------
    EntrezConfig
        A dataclass containing auth data as attributes.
    """
    config_file = Path("/tmp/gbtools_auth.json")
    if not config_file.exists():
        raise FileNotFoundError("config file does not exist")
    if not config_file.is_file():
        raise IOError("cannot read a directory")

    with config_file.open("r") as fh:
        data = json.load(fh)
        return EntrezConfig(
            email=data.get("EMAIL", "jegsamson.dev@gmail.com"),
            api_key=data.get("NCBI_API_KEY"),
        )


def fetch_genbank_handle(
    id: str, config: Optional[EntrezConfig] = None, verbose: bool = False
) -> TextIOWrapper:
    """Retrieve a Genbank handle from NCBI using an accession identifier.

    Parameters
    ----------
    id : str
        An NCBI accession identifier.
    config : EntrezConfig, optional
        A dataclass containing user metadata for identifying to the Entrez server.
    verbose : bool, default=False
        Print additional logging messages to console.

    Returns
    -------
    TextIOWrapper
        A text handle for the SeqRecord object.
    """
    if config:
        config()
    if verbose:
        print(f"Fetching Genbank record with accession {id}...")

    return Entrez.efetch(db="nuccore", id=id, rettype="gb", retmode="text")


def fetch_genbank_record(
    id: str,
    outfile: Optional[str | Path] = None,
    config: Optional[EntrezConfig] = None,
    verbose: bool = False,
    **kwargs,
) -> SeqRecord:
    """Retrieve a single Genbank record from NCBI using an accession identifier.

    Parameters
    ----------
    id : str
        An NCBI accession identifier.
    outpath : str or Path, optional
        The target filepath for saving the Genbank record.
    config : EntrezConfig, optional
        A dataclass containing user metadata for identifying to the Entrez server.
    verbose : bool, default=False
        Print additional logging messages to console.
    Returns
    -------
    SeqRecord
        A Biopython class the holds sequence data and related annotations.
    """
    # Retrieve and parse gbk handle.
    handle = fetch_genbank_handle(id=id, config=config, verbose=verbose)
    record = SeqIO.read(handle, "gb")

    # Provide a record summary.
    if verbose:
        print(f"""
Record summary:
    accession   : {record.id}
    length      : {len(record.seq)}
    description : {record.description} 
        """)

    # Save record to specified file.
    if outfile:
        if verbose:
            print(f"Attempting to write Genbank record...")
        in_handle = Path(outfile).open("w")
        SeqIO.write(record, in_handle, "gb")
        in_handle.close()
        if verbose:
            print(f"Done! Genbank record saved to {outfile}")

    return record


def load_gbk(filepath: Path | str) -> SeqRecord:
    """Load Genbank file as a SeqRecord object.

    Parameters
    ----------
    filepath : Path or str
        Path to the Genbank file.

    Returns
    -------
    SeqRecord
        A Biopython object containing the Genbank record.
    """
    fp = Path(filepath)
    if not fp.exists():
        raise FileNotFoundError("Genbank file does not exist")
    if not fp.is_file():
        raise ValueError("path does not point to a file")
    return SeqIO.read(fp, "gb")


def gbk2fasta(
    id: Optional[str] = None,
    infile: Optional[str | Path] = None,
    outfile: Optional[str | Path] = None,
    config: Optional[EntrezConfig] = None,
    verbose: bool = False,
) -> SeqRecord:
    """Extracts the sequence from a Genbank record in FASTA format.

    If both `id` and `infile` are provided, prioritize the reading of the input file.

    Parameters
    ----------
    id : str
        An NCBI accession identifier.
    infile : str or Path, optional
        Path to a Genbank file.
    outfile : str or Path, optional
        Target path for saving the FASTA-formatted file.
    config : EntrezConfig, optional
        A dataclass containing user metadata for identifying to the Entrez server.
    verbose : bool, default=False
        Print additional logging messages to console.

    Returns
    -------
    SeqRecord
        A Biopython object containing the FASTA record.
    """
    fa_handle = StringIO("")
    # Extract sequence from local gbk file.
    if infile:
        if verbose:
            print("Input file detected...")

        # Raise error if filepath is invalid.
        infile = Path(infile)
        if not infile.is_file():
            raise ValueError("provided path is not a file")

        # TODO: add error handling using try-except clause
        n_records = SeqIO.convert(infile, "gb", fa_handle, "fasta")

    # Call Entrez API for gbk handle
    else:
        if id is None:
            raise ValueError("`id` not provided")
        gbk_handle = fetch_genbank_handle(id=id, config=config, verbose=verbose)
        n_records = SeqIO.convert(gbk_handle, "gb", fa_handle, "fasta")

    # Raise error if output FASTA is empty.
    if n_records == 0:
        raise IOError("FASTA file is empty")

    if outfile:
        if verbose:
            print("Attempting to write FASTA file...")
        with Path(outfile).open("w") as out_handle:
            out_handle.write(fa_handle.getvalue())
        if verbose:
            print(f"Done! FASTA file written to {outfile}...")

    # Write FASTA to temporary file.
    tmp = NamedTemporaryFile()
    with open(tmp.name, "w") as tmp_out:
        tmp_out.write(fa_handle.getvalue())

    with open(tmp.name, "r") as tmp_out:
        fa_record = SeqIO.read(tmp_out, "fasta")
        return fa_record


def gbk2features(gbk: SeqRecord, limit_to: str = "cds"):
    """Extract all sequence features from a Genbank record.

    Parameters
    ----------
    gbk : SeqRecord
        Parsed Genbank record loaded as a Biopython object.
    limit_to: {'cds', 'gene'}, default='cds'
        Limit extraction to this gene feature.

    Returns
    -------
    dict
        A dictionary mapping feature names to sequence information.
    """
    assert vars(gbk).get("features") is not None, "invalid file format"
    assert limit_to in ["cds", "gene"], "invalid value, expecting `CDS` or `gene`"

    features_dict = {}
    for f in gbk.features:
        f_name = f.qualifiers.get("gene", [None])[0]
        if f_name is None:
            continue
        # Extract nucleotide sequence and location
        f_data = dict(
            nuc_seq=str(f.extract(gbk.seq)),
            pos_start=int(f.location.start),
            pos_end=int(f.location.end),
            strand=f.location.strand,
        )
        if limit_to == "cds" and f.type == "CDS":
            f_data = dict(
                prot_seq=f.qualifiers.get("translation", [""])[0],
                prot_id=f.qualifiers.get("protein_id", [None])[0],
                prot_name=f.qualifiers.get("product", [None])[0],
                **f_data,
            )
        features_dict[f_name] = f_data

    return features_dict


def print_fasta(gbk: SeqRecord, linewidth: int = 60) -> None:
    """Print sequence from SeqRecord in FASTA format.

    Parameters
    ----------
    gbk : SeqRecord
        A Genbank recorded loaded into memory as a SeqRecord object.
    linewidth : int, default=60
        Number of characters per line in the FASTA file.

    Returns
    -------
    None
    """
    assert isinstance(gbk, SeqRecord), "first argument must be a SeqRecord object"
    print(f">{gbk.description}")
    for i in range(0, len(gbk.seq), linewidth):
        print(str(gbk.seq)[i : i + linewidth])
