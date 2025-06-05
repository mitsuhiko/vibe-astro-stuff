import click

from .iss_transit import calculate_iss_moon_transits


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
    transits = calculate_iss_moon_transits(latitude, longitude, elevation, max_distance)

    if not transits:
        click.echo("\nNo ISS-Moon transits or close passes found in the next year from your location.")
        return

    click.echo(f"\nFound {len(transits)} ISS-Moon transit(s) and close pass(es):\n")

    for i, transit in enumerate(transits, 1):
        # Determine event type
        sep_deg = transit['separation_degrees']
        moon_radius = transit['moon_radius_degrees']
        
        if sep_deg < moon_radius:
            event_type = "Transit"
        else:
            event_type = "Close pass"
            
        # Convert angular separation to arcminutes and arcseconds
        sep_arcmin = sep_deg * 60
        sep_deg_int = int(sep_deg)
        sep_arcmin_int = int((sep_deg - sep_deg_int) * 60)
        sep_arcsec = ((sep_deg - sep_deg_int) * 60 - sep_arcmin_int) * 60
        
        click.echo(f"{event_type} #{i}:")
        click.echo(f"  Date/Time (UTC): {transit['time'].strftime('%a %Y-%m-%d, %H:%M:%S.%f')[:-4]}")
        click.echo(f"  ISS angular size: {transit['iss_angular_size_arcmin'] * 60:.2f}″; distance: {transit['iss_distance_km']:.2f} km")
        click.echo(f"  Angular separation: {sep_deg_int}° {sep_arcmin_int}′ {sep_arcsec:.0f}″; azimuth: {transit['iss_azimuth']:.1f}°; altitude: {transit['iss_altitude']:.1f}°")
        click.echo(f"  Moon illumination: {transit['moon_illumination']:.1f}%")
        click.echo()
