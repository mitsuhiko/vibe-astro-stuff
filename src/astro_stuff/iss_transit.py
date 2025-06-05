import os
import json
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any

import numpy as np
import pytz
import requests
from skyfield import almanac
from skyfield.api import Topos, load, EarthSatellite, wgs84
from skyfield.positionlib import ICRF
from skyfield.units import Angle


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
        from skyfield.api import Loader
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
        moon_radius_km = 1737.4  # Moon's radius in km
        angular_radius_rad = np.arctan(moon_radius_km / moon_distance)
        return np.degrees(angular_radius_rad)

    def _is_iss_sunlit(self, t) -> bool:
        """Check if ISS is sunlit at given time."""
        sun = self.ephemeris["sun"]
        iss_position = self.iss_satellite.at(t)
        sun_position = (sun - self.earth).at(t)

        # Get ISS position relative to Earth center
        iss_vector = iss_position.position.km
        sun_vector = sun_position.position.au * 149597870.7  # Convert AU to km

        # Check if ISS is in Earth's shadow
        # Simplified check: if angle between ISS and Sun vectors > 90 degrees
        # and ISS is closer than Earth's shadow cone
        dot_product = np.dot(iss_vector, sun_vector)
        if dot_product < 0:  # ISS is on opposite side of Earth from Sun
            # Calculate if ISS is within Earth's shadow cone
            iss_distance = np.linalg.norm(iss_vector)
            earth_radius = 6371.0  # km
            if iss_distance < earth_radius * 1.01:  # Add small margin
                return False
        return True

    def find_transits(
        self, start_date: datetime, end_date: datetime, max_distance_km: Optional[float] = None, verbose: bool = False
    ) -> List[Dict[str, Any]]:
        """Find ISS-Moon transit events within the given time range."""
        transits = []

        # Convert to Skyfield time
        t0 = self.ts.from_datetime(start_date.replace(tzinfo=pytz.UTC))
        t1 = self.ts.from_datetime(end_date.replace(tzinfo=pytz.UTC))

        # Coarse-to-fine search for speed
        print("Phase 1: Coarse search (1-hour intervals)...")
        
        # Phase 1: Coarse search every hour
        coarse_step = 3600.0 / 86400.0  # 1 hour in days
        coarse_samples = int((t1.tt - t0.tt) / coarse_step)
        
        promising_times = []
        
        above_horizon_count = 0
        
        for i in range(coarse_samples):
            t = self.ts.tt_jd(t0.tt + i * coarse_step)
            
            # Quick check: calculate separation
            difference = self.iss_satellite - self.observer
            iss_relative = difference.at(t)
            moon_astrometric = (self.earth + self.observer).at(t).observe(self.moon)
            moon_relative = moon_astrometric.apparent()
            
            # Check horizon
            iss_alt, _, _ = iss_relative.altaz()
            moon_alt, _, _ = moon_relative.altaz()
            
            separation = self._angular_separation(iss_relative, moon_relative)
            
            # Count when both are above horizon for reference
            if iss_alt.degrees > 0 and moon_alt.degrees > 0:
                above_horizon_count += 1
            
            # Debug: print some samples regardless of horizon
            if verbose and i % 24 == 0:  # Every day
                print(f"  {t.utc_datetime().strftime('%m-%d %H:%M')}: ISS {iss_alt.degrees:.1f}° Moon {moon_alt.degrees:.1f}° Sep {separation:.1f}°")
            
            # Mark as promising if separation is small (regardless of horizon for coarse search)
            if separation < 30.0:
                promising_times.append(t)
                if verbose:
                    print(f"  Promising: {t.utc_datetime().strftime('%m-%d %H:%M')} - {separation:.1f}° (ISS {iss_alt.degrees:.1f}° Moon {moon_alt.degrees:.1f}°)")
        
        if verbose:
            print(f"  Checked {coarse_samples} hours, {above_horizon_count} with both above horizon")
        
        print(f"Found {len(promising_times)} promising hour(s)")
        
        if not promising_times:
            return []
        
        # Phase 2: Fine search around promising times
        print("Phase 2: Fine search (2-minute intervals around promising times)...")
        close_approaches = []
        
        for promising_time in promising_times:
            # Search ±2 hours around this time with 2-minute resolution
            start_fine = promising_time.tt - (2.0 / 24.0)  # 2 hours before
            end_fine = promising_time.tt + (2.0 / 24.0)    # 2 hours after
            
            fine_step = 120.0 / 86400.0  # 2 minutes in days
            fine_samples = int((end_fine - start_fine) / fine_step)
            
            for j in range(fine_samples):
                t = self.ts.tt_jd(start_fine + j * fine_step)
                
                difference = self.iss_satellite - self.observer
                iss_relative = difference.at(t)
                moon_astrometric = (self.earth + self.observer).at(t).observe(self.moon)
                moon_relative = moon_astrometric.apparent()
                
                iss_alt, _, _ = iss_relative.altaz()
                moon_alt, _, _ = moon_relative.altaz()
                
                if iss_alt.degrees < 5 or moon_alt.degrees < 5:
                    continue
                
                separation = self._angular_separation(iss_relative, moon_relative)
                
                if separation < 3.0:
                    close_approaches.append({
                        'time': t,
                        'separation': separation,
                        'iss_alt': iss_alt.degrees,
                        'moon_alt': moon_alt.degrees
                    })
        
        print(f"Found {len(close_approaches)} close approach time points")
        
        # Find local minima (actual transit events)
        transits = []
        if close_approaches:
            transits = self._find_transit_events(close_approaches)
        
        print(f"\nSearch complete! Found {len(transits)} event(s)")
        return transits

    def _find_transit_events(self, close_approaches: List[Dict]) -> List[Dict[str, Any]]:
        """Find actual transit events from close approach data points."""
        if not close_approaches:
            return []
        
        # Sort by time
        close_approaches.sort(key=lambda x: x['time'].tt)
        
        events = []
        i = 0
        
        while i < len(close_approaches):
            # Start a new event group
            group_start = i
            best_idx = i
            best_separation = close_approaches[i]['separation']
            
            # Find all points within 2 hours of this point
            while (i + 1 < len(close_approaches) and 
                   (close_approaches[i + 1]['time'].tt - close_approaches[group_start]['time'].tt) * 24 < 2):
                i += 1
                # Track the best (minimum separation) in this group
                if close_approaches[i]['separation'] < best_separation:
                    best_separation = close_approaches[i]['separation']
                    best_idx = i
            
            # Create transit event for the best point in this group
            best_time = close_approaches[best_idx]['time']
            events.append(self._create_transit_event(best_time, best_separation))
            
            i += 1
        
        return events

    def _create_transit_event(self, t, separation: float) -> Dict[str, Any]:
        """Create a transit event dictionary with all relevant information."""
        # Get positions at transit time
        difference = self.iss_satellite - self.observer
        iss_relative = difference.at(t)
        
        moon_astrometric = (self.earth + self.observer).at(t).observe(self.moon)
        moon_relative = moon_astrometric.apparent()

        # Get altitude and azimuth
        iss_alt, iss_az, _ = iss_relative.altaz()
        moon_alt, moon_az, _ = moon_relative.altaz()

        # Calculate ISS angular size
        iss_distance = iss_relative.distance().km
        iss_size_arcsec = 109 / (iss_distance * 1000) * 206265  # ISS is ~109m wide, convert km to m

        # Get Moon phase
        sun = self.ephemeris["sun"]
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
            "moon_illumination": (1 + np.cos(np.radians(phase_angle.degrees))) / 2 * 100,
        }


def calculate_iss_moon_transits(
    latitude: float,
    longitude: float,
    elevation: float,
    max_distance_km: Optional[float] = None,
    verbose: bool = False,
) -> List[Dict[str, Any]]:
    """Calculate ISS-Moon transits for the next year from given location."""
    calculator = ISSMoonTransitCalculator(latitude, longitude, elevation)

    # Search for transits in the next 30 days (start from beginning of today)
    start_date = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
    end_date = start_date + timedelta(days=30)

    print(f"\nSearching for ISS-Moon transits from {latitude}°, {longitude}°, {elevation}m")
    print(f"Time range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")

    if max_distance_km:
        print(f"Maximum travel distance: {max_distance_km} km")
        print("(Travel distance feature not yet implemented)")

    transits = calculator.find_transits(start_date, end_date, max_distance_km, verbose)

    return transits