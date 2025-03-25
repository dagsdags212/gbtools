import click

from gbtools.Genbank import Genbank, GenbankRecord


@click.group()
def cli():
    pass


@click.option('--id', help='genbank accession')
@click.option('--out', help='path for storing the GBK record')
@click.command()
def fetch(id, out):
    click.echo(f"Fetching {id} from Genbank")
    Genbank.fetch(id, out)

cli.add_command(fetch)

if __name__ == "__main__":
    cli()
