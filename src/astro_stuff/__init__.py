import click


@click.group()
def main() -> None:
    """Astro-stuff: A CLI tool for celestial calculations."""
    pass


@main.command()
def hello() -> None:
    """Say hello from astro-stuff."""
    click.echo("Hello from astro-stuff!")
