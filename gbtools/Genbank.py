import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import polars as pl
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich import print as rprint
from rich.tree import Tree


@dataclass
class Gene:
    name: str
    source_id: str
    protein_id: str
    seq: Seq | str
    translation: Seq | str
    start: int
    end: int
    strand: int

    def __len__(self) -> int:
        return self.end - self.start

    def __str__(self) -> str:
        return f"Gene(name={self.name},source_id={self.source_id},start={self.start},end={self.end})"

    def to_dict(self) -> dict[str, str]:
        """Return attributes as a dictionary."""
        return dict(
            name=self.name,
            source_id=self.source_id,
            seq=self.seq,
            protein_id=self.protein_id,
            translation=self.translation,
            start=self.start,
            end=self.end,
            strand=self.strand,
        )

    def plot(self, context: int = 1) -> None:
        """Graph the gene using pygenomeviz."""
        raise NotImplementedError


class GenbankRecord:
    def __init__(self, gbk: Path | SeqRecord) -> None:
        if isinstance(gbk, (Path, str)):
            self.data = SeqIO.read(gbk, "genbank")
        elif isinstance(gbk, SeqRecord):
            self.data = gbk
        else:
            raise ValueError("format not supported")

    def __repr__(self) -> str:
        return f"GenbankRecord(id={self.id},description={self.description})"

    def __str__(self) -> str:
        return f"GenbankRecord(id={self.id},description={self.description})"

    @property
    def seq(self) -> Seq:
        return self.data.seq

    @property
    def id(self) -> str:
        return self.data.id

    @property
    def description(self) -> str:
        return self.data.description

    @property
    def annotations(self) -> dict[str, str]:
        return self.data.annotations

    def get_features(self, feat_type: str="CDS"):
        assert feat_type in ["CDS", "gene"], f"Error: invalid feature type, must be `CDS` or `gene` (feat_type={feat_type})"
        features = []
        for f in self.data.features:
            if f.type == feat_type:
                start, end = int(f.location.start), int(f.location.end)
                gene = Gene(
                    name=f.qualifiers.get("gene", [""])[0],
                    source_id=self.data.id,
                    seq=self.seq[start:end],
                    protein_id=f.qualifiers.get("protein_id", [""])[0],
                    translation=f.qualifiers.get("translation", [""])[0],
                    start=start,
                    end=end,
                    strand=f.location.strand,
                )
                features.append(gene)

        return features

    def extract_seq(self, gene: str, moltype: str="DNA") -> Optional[str]:
        """Extract the DNA or protein sequence of a gene if annotation is present in the record."""
        if moltype not in ["DNA", "PROTEIN"]:
            raise ValueError(f"Invalid moltype, expecting `DNA` or `PROTEIN` (moltype={moltype})")

        cds_table = self.tabulate("CDS")
        # Exit if gene is not in the table.
        if gene not in cds_table["name"]:
            print(f"Sequence not found (gene={gene})")
            return

        # Return the longer sequence if multiple hits were found.
        if moltype == "DNA":
            gene_hits = cds_table.filter(pl.col("name") == gene)["seq"]
            return str(max(gene_hits.to_list(), key=len))
        elif moltype == "PROTEIN":
            prot_hits = cds_table.filter(pl.col("name") == gene)["translation"]
            return max(prot_hits.to_list(), key=len)


    def tabulate(self, feat_type: str="CDS") -> pl.DataFrame:
        """Print a table for annotated genes"""
        assert feat_type in ["CDS", "gene"], f"Error: invalid feature type, must be `CDS` or `gene` (feat_type={feat_type})"
        cds = self.get_features(feat_type)
        data = {f: [] for f in cds[0].to_dict().keys()}
        for gene in cds:
            for k, v in vars(gene).items():
                data[k].append(v)
        return pl.DataFrame(data)

    def gene_tree(self) -> None:
        """Print a tree hierarchy for annotated genes"""
        root = Tree(self.id)
        for gene in self.get_features("CDS"):
            child = Tree(gene.name)
            child.add(f"id: {gene.source_id}")
            child.add(f"length: {gene.end - gene.start}")
            child.add(f"position: {gene.start} to {gene.end}")
            child.add(f"strand: {gene.strand}")
            child.add(f"protein_id: {gene.protein_id}")
            root.add(child)
        rprint(root)

    def to_dict(self) -> dict[str, Gene]:
        return {gene.name.replace(" CDS", ""): gene for gene in self.get_features("CDS")}

class Genbank:
    """
    A Python interface for downloading genbank files from NCBI.
    """

    @classmethod
    def authenticate(
        cls,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
    ) -> None:
        """Optionally provide Entrez server with user information."""
        if email:
            os.environ["EMAIL"] = email
            print(f"Email registered as {email}")
        if api_key:
            os.environ["NCBI_API_KEY"] = api_key
            print("API key has been set")

    @classmethod
    def info(cls, id: str) -> None:
        """Print metadata associated with a GenBank record."""
        # Authenticate with server.
        email = os.environ.get("EMAIL", "jegsamson.dev@gmail.com")
        api_key = os.environ.get("NCBI_API_KEY")
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        handle = Entrez.esummary(db="nuccore", id=id)
        summary = Entrez.read(handle)[0]

        json.dump(summary, sys.stdout)

    @classmethod
    def load(cls, id: str) -> SeqRecord:
        """Loads a Genbank record in memory."""
        email = os.environ.get("EMAIL", "jegsamson.dev@gmai.com")
        api_key = os.environ.get("NCBI_API_KEY")
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        # Retrieve and parse gbk file from NCBI.
        ehandle = Entrez.efetch(
            db="nuccore", id=id, rettype="genbank", retmode="text"
        )
        return SeqIO.read(ehandle, "genbank")


    @classmethod
    def fetch(cls, id: str, outpath: Optional[str | Path] = None) -> Path:
        """Retrieve a local copy of a GenBank record from NCBI."""
        # Authenticate with server.
        email = os.environ.get("EMAIL", "jegsamson.dev@gmail.com")
        api_key = os.environ.get("NCBI_API_KEY")
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

        # Retrieve and parse gbk file from NCBI.
        ehandle = Entrez.efetch(
            db="nuccore", id=id, rettype="genbank", retmode="text"
        )
        gbk = SeqIO.read(ehandle, "genbank")

        # Defaults to accession number.
        outpath = Path(outpath) if outpath else Path(f"{id}.gbk")

        if outpath.exists():
            overwrite = input(
                "Attempting to overwrite an existing file, continue? [Y/n]"
            )
            match overwrite:
                case "y" | "Y":
                    SeqIO.write(gbk, outpath, "genbank")
                case "n" | "N":
                    print("Exiting program.")
                    sys.exit(1)
                case _:
                    overwrite = input(
                        "Attempting to overwrite an existing file, continue? [Y/n]"
                    )
        else:
            SeqIO.write(gbk, outpath, "genbank")

        print(f"File written to {outpath.absolute()}")
        return outpath.absolute()
