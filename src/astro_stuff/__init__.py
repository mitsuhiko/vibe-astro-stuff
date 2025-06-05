import click


@click.group()
def main() -> None:
    """Astro-stuff: A CLI tool for celestial calculations."""
    pass


@main.command()
def hello() -> None:
    """Say hello from astro-stuff."""
    click.echo("Hello from astro-stuff!")


@main.command()
@click.option(
    "--latitude", "-lat", type=float, required=True, help="Observer latitude in degrees"
)
@click.option(
    "--longitude",
    "-lon",
    type=float,
    required=True,
    help="Observer longitude in degrees",
)
@click.option(
    "--elevation",
    "-elev",
    type=float,
    required=True,
    help="Observer elevation in meters",
)
@click.option(
    "--max-distance",
    "-d",
    type=float,
    default=None,
    help="Maximum preferred travel distance in km",
)
def iss_moon_transit(
    latitude: float, longitude: float, elevation: float, max_distance: float | None
) -> None:
    """Calculate when the ISS transits in front of the moon from your location."""
    click.echo(
        f"Calculating ISS-moon transits for location: {latitude}°, {longitude}°, {elevation}m"
    )
    if max_distance:
        click.echo(f"Maximum travel distance: {max_distance} km")

    # TODO: Implement ISS-moon transit calculations
    click.echo("Transit calculations not yet implemented.")
