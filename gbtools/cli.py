import click
from rich.pretty import pprint

from gbtools.Genbank import Genbank, GenbankRecord


@click.group()
def cli():
    pass

@click.option('--email', type=str, help='email address')
@click.option('--api-key', type=str, default=None ,help='NCBI API key')
@click.command()
def auth(email, api_key):
    """Specify email and API key to be informed of API rate limits"""
    Genbank.authenticate(email=email, api_key=api_key)

@click.argument('id')
@click.option('--out', help='path for storing the GBK record')
@click.command()
def fetch(id: str, out):
    """Download a Genbank record from an accession"""
    click.echo(f"Fetching {id} from Genbank")
    Genbank.fetch(id, out)

@click.argument('id')
@click.command()
def info(id):
    """Retrieve metadata associated with a Genbank record"""
    Genbank.info(id)

@click.argument('id')
@click.option('--tree', is_flag=True, type=bool, show_default=False, default=False, 
              help="display in tree format")
@click.command()
def list_genes(id, tree):
    """List all annotated genes in tabular form"""
    gbk = GenbankRecord(Genbank.load(id))
    if tree:
        gbk.gene_tree()
    else:
        click.echo(gbk.tabulate())


cli.add_command(auth)
cli.add_command(fetch)
cli.add_command(info)
cli.add_command(list_genes)

if __name__ == "__main__":
    cli()
