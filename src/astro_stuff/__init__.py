import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

import click
import numpy as np
import pytz
import requests
from skyfield import almanac
from skyfield.api import EarthSatellite, wgs84, Loader
from skyfield.positionlib import ICRF


# Constants
MOON_RADIUS_KM = 1737.4
ISS_WIDTH_METERS = 109
EARTH_RADIUS_KM = 6371.0
COARSE_SEARCH_THRESHOLD_DEGREES = 30.0
FINE_SEARCH_THRESHOLD_DEGREES = 3.0
MIN_ALTITUDE_DEGREES = 5.0

# ISS orbital constants
ISS_ORBITAL_PERIOD_MINUTES = 93.0  # ISS orbital period in minutes
ISS_ORBITS_PER_DAY = 1440 / ISS_ORBITAL_PERIOD_MINUTES  # ~15.5 orbits/day
MOON_ANGULAR_SPEED_DEG_PER_HOUR = 0.5  # Moon moves ~0.5° per hour


class ISSMoonTransitCalculator:
    """Calculates ISS-Moon transit events from a given observer location."""

    TLE_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle"
    TLE_CACHE_DIR = Path.home() / ".cache" / "astro-stuff"
    TLE_CACHE_FILE = TLE_CACHE_DIR / "iss_tle.json"
    TLE_MAX_AGE_DAYS = 7

    def __init__(self, latitude: float, longitude: float, elevation: float):
        """Initialize calculator with observer location."""
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation

        # Create observer location using wgs84
        self.observer = wgs84.latlon(latitude, longitude, elevation_m=elevation)

        # Load Skyfield data with local caching
        loader = Loader(".")  # Use current directory for cache
        self.ts = loader.timescale()
        self.ephemeris = loader("de421.bsp")
        self.earth = self.ephemeris["earth"]
        self.moon = self.ephemeris["moon"]

        # Load ISS TLE data
        self.iss_satellite = self._load_iss_tle()

    def _fetch_fresh_tle(self) -> Tuple[str, str]:
        """Fetch fresh TLE data from Celestrak."""
        response = requests.get(self.TLE_URL, timeout=10)
        response.raise_for_status()

        lines = response.text.strip().split("\n")
        # Find ISS (ZARYA) in the TLE data
        for i in range(0, len(lines), 3):
            if "ISS (ZARYA)" in lines[i]:
                return lines[i + 1], lines[i + 2]

        raise ValueError("ISS TLE not found in Celestrak data")

    def _load_iss_tle(self) -> EarthSatellite:
        """Load ISS TLE data, using cache if available and fresh."""
        self.TLE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

        # Check if cached TLE exists and is fresh
        if self.TLE_CACHE_FILE.exists():
            with open(self.TLE_CACHE_FILE, "r") as f:
                cache_data = json.load(f)

            cache_time = datetime.fromisoformat(cache_data["timestamp"])
            if datetime.now() - cache_time < timedelta(days=self.TLE_MAX_AGE_DAYS):
                # Cache is fresh, use it
                line1 = cache_data["line1"]
                line2 = cache_data["line2"]
                print(
                    f"Using cached ISS TLE data from {cache_time.strftime('%Y-%m-%d %H:%M:%S UTC')}"
                )
            else:
                # Cache is stale
                print("Cached TLE data is stale, fetching fresh data...")
                line1, line2 = self._fetch_fresh_tle()
                self._save_tle_cache(line1, line2)
        else:
            # No cache, fetch fresh
            print("Fetching ISS TLE data from Celestrak...")
            line1, line2 = self._fetch_fresh_tle()
            self._save_tle_cache(line1, line2)

        return EarthSatellite(line1, line2, "ISS (ZARYA)", self.ts)

    def _save_tle_cache(self, line1: str, line2: str) -> None:
        """Save TLE data to cache."""
        cache_data = {
            "timestamp": datetime.now().isoformat(),
            "line1": line1,
            "line2": line2,
        }
        with open(self.TLE_CACHE_FILE, "w") as f:
            json.dump(cache_data, f)

    def _angular_separation(self, pos1: ICRF, pos2: ICRF) -> float:
        """Calculate angular separation between two positions in degrees."""
        # Get unit vectors
        r1 = pos1.position.au / np.linalg.norm(pos1.position.au)
        r2 = pos2.position.au / np.linalg.norm(pos2.position.au)

        # Calculate angle using dot product
        cos_angle = np.clip(np.dot(r1, r2), -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)

    def _moon_radius_degrees(self, t) -> float:
        """Calculate apparent radius of moon in degrees at given time."""
        moon_distance = (self.moon - self.earth).at(t).distance().km
        angular_radius_rad = np.arctan(MOON_RADIUS_KM / moon_distance)
        return np.degrees(angular_radius_rad)

    def _get_positions(self, t):
        """Get ISS and Moon positions at given time."""
        iss_relative = (self.iss_satellite - self.observer).at(t)
        moon_astrometric = (self.earth + self.observer).at(t).observe(self.moon)
        moon_relative = moon_astrometric.apparent()
        return iss_relative, moon_relative

    def find_transits(
        self,
        start_date: datetime,
        end_date: datetime,
        max_distance_km: Optional[float] = None,
        verbose: bool = False,
    ) -> List[Dict[str, Any]]:
        """Find ISS-Moon transit events using fast orbital intersection algorithm."""
        print("Using fast orbital intersection algorithm...")
        
        # Convert to Skyfield time
        t0 = self.ts.from_datetime(start_date.replace(tzinfo=pytz.UTC))
        t1 = self.ts.from_datetime(end_date.replace(tzinfo=pytz.UTC))
        
        # Fast orbital-based search
        return self._fast_orbital_search(t0, t1, verbose)
    
    def _fast_orbital_search(self, t0, t1, verbose: bool) -> List[Dict[str, Any]]:
        """Fast search based on ISS orbital periods and Moon trajectory prediction."""
        
        # Step 1: Sample ISS positions at orbital intervals
        orbital_period_days = ISS_ORBITAL_PERIOD_MINUTES / (24 * 60)
        total_days = t1.tt - t0.tt
        num_orbits = int(total_days / orbital_period_days)
        
        print(f"Analyzing {num_orbits} ISS orbits over {total_days:.1f} days...")
        
        # Sample key points in each orbit (4 points per orbit: quadrants)
        iss_samples = []
        for orbit in range(num_orbits):
            for phase in [0.0, 0.25, 0.5, 0.75]:  # Orbit quadrants
                t_orbit = self.ts.tt_jd(t0.tt + (orbit + phase) * orbital_period_days)
                iss_samples.append(t_orbit)
        
        if verbose:
            print(f"Generated {len(iss_samples)} ISS sample points")
        
        # Step 2: Pre-compute Moon trajectory checkpoints
        moon_checkpoints = []
        checkpoint_interval_hours = 3  # Moon position every 3 hours
        checkpoint_interval_days = checkpoint_interval_hours / 24.0
        
        num_checkpoints = int(total_days / checkpoint_interval_days)
        for i in range(num_checkpoints + 1):
            t_checkpoint = self.ts.tt_jd(t0.tt + i * checkpoint_interval_days)
            moon_astrometric = (self.earth + self.observer).at(t_checkpoint).observe(self.moon)
            moon_relative = moon_astrometric.apparent()
            moon_alt, moon_az, _ = moon_relative.altaz()
            
            moon_checkpoints.append({
                'time': t_checkpoint,
                'alt': moon_alt.degrees,
                'az': moon_az.degrees,
                'position': moon_relative
            })
        
        if verbose:
            print(f"Pre-computed {len(moon_checkpoints)} Moon trajectory points")
        
        # Step 3: Fast intersection detection
        promising_times = []
        
        for iss_time in iss_samples:
            # Get ISS position
            iss_relative = (self.iss_satellite - self.observer).at(iss_time)
            iss_alt, iss_az, _ = iss_relative.altaz()
            
            # Skip if ISS below horizon
            if iss_alt.degrees < MIN_ALTITUDE_DEGREES:
                continue
            
            # Find nearest Moon checkpoint and interpolate
            moon_pos = self._interpolate_moon_position(iss_time, moon_checkpoints)
            if moon_pos is None or moon_pos['alt'] < MIN_ALTITUDE_DEGREES:
                continue
            
            # Quick angular separation check
            alt_diff = abs(iss_alt.degrees - moon_pos['alt'])
            az_diff = abs(iss_az.degrees - moon_pos['az'])
            if az_diff > 180:
                az_diff = 360 - az_diff
            
            # Rough angular distance (faster than full spherical calculation)
            rough_separation = np.sqrt(alt_diff**2 + (az_diff * np.cos(np.radians(moon_pos['alt'])))**2)
            
            if rough_separation < COARSE_SEARCH_THRESHOLD_DEGREES:
                promising_times.append(iss_time)
                if verbose:
                    print(f"  Promising: {iss_time.utc_datetime().strftime('%m-%d %H:%M')} - rough sep {rough_separation:.1f}°")
        
        print(f"Fast scan found {len(promising_times)} promising orbital positions")
        
        if not promising_times:
            return []
        
        # Step 4: Precise refinement only around promising times
        print("Refining promising candidates...")
        return self._refine_candidates(promising_times, verbose)
    
    def _interpolate_moon_position(self, target_time, moon_checkpoints):
        """Interpolate Moon position at target time from pre-computed checkpoints."""
        if not moon_checkpoints:
            return None
        
        # Find surrounding checkpoints
        before = None
        after = None
        
        for checkpoint in moon_checkpoints:
            if checkpoint['time'].tt <= target_time.tt:
                before = checkpoint
            elif checkpoint['time'].tt > target_time.tt:
                after = checkpoint
                break
        
        if before is None:
            return moon_checkpoints[0]
        if after is None:
            return moon_checkpoints[-1]
        
        # Linear interpolation
        time_fraction = (target_time.tt - before['time'].tt) / (after['time'].tt - before['time'].tt)
        
        # Interpolate altitude and azimuth
        alt_interp = before['alt'] + time_fraction * (after['alt'] - before['alt'])
        
        # Handle azimuth wraparound
        az_diff = after['az'] - before['az']
        if az_diff > 180:
            az_diff -= 360
        elif az_diff < -180:
            az_diff += 360
        
        az_interp = before['az'] + time_fraction * az_diff
        if az_interp < 0:
            az_interp += 360
        elif az_interp >= 360:
            az_interp -= 360
        
        return {
            'alt': alt_interp,
            'az': az_interp,
            'time': target_time
        }
    
    def _refine_candidates(self, promising_times, verbose: bool) -> List[Dict[str, Any]]:
        """Precisely refine promising orbital positions."""
        close_approaches = []
        
        for promising_time in promising_times:
            # Search ±10 minutes around this orbital position with 30-second resolution
            search_window_minutes = 10
            search_window_days = search_window_minutes / (24 * 60)
            
            start_refine = promising_time.tt - search_window_days
            end_refine = promising_time.tt + search_window_days
            
            refine_step = 30.0 / 86400.0  # 30 seconds in days
            refine_samples = int((end_refine - start_refine) / refine_step)
            
            best_separation = float('inf')
            best_time = None
            
            for j in range(refine_samples):
                t = self.ts.tt_jd(start_refine + j * refine_step)
                
                iss_relative, moon_relative = self._get_positions(t)
                iss_alt, _, _ = iss_relative.altaz()
                moon_alt, _, _ = moon_relative.altaz()
                
                if (iss_alt.degrees < MIN_ALTITUDE_DEGREES or 
                    moon_alt.degrees < MIN_ALTITUDE_DEGREES):
                    continue
                
                separation = self._angular_separation(iss_relative, moon_relative)
                
                if separation < best_separation:
                    best_separation = separation
                    best_time = t
                
                if separation < FINE_SEARCH_THRESHOLD_DEGREES:
                    close_approaches.append({
                        "time": t,
                        "separation": separation,
                        "iss_alt": iss_alt.degrees,
                        "moon_alt": moon_alt.degrees,
                    })
            
            if verbose and best_time is not None:
                print(f"  Best in window: {best_time.utc_datetime().strftime('%m-%d %H:%M:%S')} - {best_separation:.2f}°")
        
        print(f"Refinement found {len(close_approaches)} close approach points")
        
        # Find local minima (actual transit events)
        if close_approaches:
            return self._find_transit_events(close_approaches)
        
        return []

    def _find_transit_events(
        self, close_approaches: List[Dict]
    ) -> List[Dict[str, Any]]:
        """Find actual transit events from close approach data points."""
        if not close_approaches:
            return []

        # Sort by time
        close_approaches.sort(key=lambda x: x["time"].tt)

        events = []
        i = 0

        while i < len(close_approaches):
            # Start a new event group
            group_start = i
            best_idx = i
            best_separation = close_approaches[i]["separation"]

            # Find all points within 2 hours of this point
            while (
                i + 1 < len(close_approaches)
                and (
                    close_approaches[i + 1]["time"].tt
                    - close_approaches[group_start]["time"].tt
                )
                * 24
                < 2
            ):
                i += 1
                # Track the best (minimum separation) in this group
                if close_approaches[i]["separation"] < best_separation:
                    best_separation = close_approaches[i]["separation"]
                    best_idx = i

            # Create transit event for the best point in this group
            best_time = close_approaches[best_idx]["time"]
            events.append(self._create_transit_event(best_time, best_separation))

            i += 1

        return events

    def _create_transit_event(self, t, separation: float) -> Dict[str, Any]:
        """Create a transit event dictionary with all relevant information."""
        # Get positions at transit time
        iss_relative, moon_relative = self._get_positions(t)

        # Get altitude and azimuth
        iss_alt, iss_az, _ = iss_relative.altaz()
        moon_alt, moon_az, _ = moon_relative.altaz()

        # Calculate ISS angular size
        iss_distance = iss_relative.distance().km
        iss_size_arcsec = (
            ISS_WIDTH_METERS / (iss_distance * 1000) * 206265
        )  # Convert km to m for calculation

        # Get Moon phase
        phase_angle = almanac.moon_phase(self.ephemeris, t)

        return {
            "time": t.utc_datetime(),
            "separation_degrees": separation,
            "moon_radius_degrees": self._moon_radius_degrees(t),
            "iss_altitude": iss_alt.degrees,
            "iss_azimuth": iss_az.degrees,
            "moon_altitude": moon_alt.degrees,
            "moon_azimuth": moon_az.degrees,
            "iss_distance_km": iss_distance,
            "iss_angular_size_arcsec": iss_size_arcsec,
            "moon_phase_degrees": phase_angle.degrees,
            "moon_illumination": (1 + np.cos(np.radians(phase_angle.degrees)))
            / 2
            * 100,
        }


def calculate_iss_moon_transits(
    latitude: float,
    longitude: float,
    elevation: float,
    max_distance_km: Optional[float] = None,
    verbose: bool = False,
    start_date: Optional[datetime] = None,
) -> List[Dict[str, Any]]:
    """Calculate ISS-Moon transits for the next year from given location."""
    calculator = ISSMoonTransitCalculator(latitude, longitude, elevation)

    # Search for transits in the next 30 days (start from beginning of today)
    if start_date is None:
        start_date = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
    end_date = start_date + timedelta(days=30)

    print(
        f"\nSearching for ISS-Moon transits from {latitude}°, {longitude}°, {elevation}m"
    )
    print(
        f"Time range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}"
    )

    if max_distance_km:
        print(f"Maximum travel distance: {max_distance_km} km")
        print("(Travel distance feature not yet implemented)")

    transits = calculator.find_transits(start_date, end_date, max_distance_km, verbose)

    return transits


@click.command()
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
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show verbose debug information during search",
)
def main(
    latitude: float,
    longitude: float,
    elevation: float,
    max_distance: float | None,
    verbose: bool,
) -> None:
    """Calculate when the ISS transits in front of the moon from your location."""
    transits = calculate_iss_moon_transits(
        latitude, longitude, elevation, max_distance, verbose
    )

    if not transits:
        click.echo(
            "\nNo ISS-Moon transits or close passes found in the next year from your location."
        )
        return

    click.echo(f"\nFound {len(transits)} ISS-Moon transit(s) and close pass(es):\n")

    for i, transit in enumerate(transits, 1):
        # Determine event type
        sep_deg = transit["separation_degrees"]
        moon_radius = transit["moon_radius_degrees"]

        if sep_deg < moon_radius:
            event_type = "Transit"
        else:
            event_type = "Close pass"

        # Convert angular separation to arcminutes and arcseconds
        sep_deg_int = int(sep_deg)
        sep_arcmin_int = int((sep_deg - sep_deg_int) * 60)
        sep_arcsec = ((sep_deg - sep_deg_int) * 60 - sep_arcmin_int) * 60

        click.echo(f"{event_type} #{i}:")
        click.echo(
            f"  Date/Time (UTC): {transit['time'].strftime('%a %Y-%m-%d, %H:%M:%S.%f')[:-4]}"
        )
        click.echo(
            f"  ISS angular size: {transit['iss_angular_size_arcsec']:.2f}″; distance: {transit['iss_distance_km']:.2f} km"
        )
        click.echo(
            f"  Angular separation: {sep_deg_int}° {sep_arcmin_int}′ {sep_arcsec:.0f}″; azimuth: {transit['iss_azimuth']:.1f}°; altitude: {transit['iss_altitude']:.1f}°"
        )
        click.echo(f"  Moon illumination: {transit['moon_illumination']:.1f}%")
        click.echo()
