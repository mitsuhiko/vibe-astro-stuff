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
        self, start_date: datetime, end_date: datetime, max_distance_km: Optional[float] = None
    ) -> List[Dict[str, Any]]:
        """Find ISS-Moon transit events within the given time range."""
        transits = []

        # Convert to Skyfield time
        t0 = self.ts.from_datetime(start_date.replace(tzinfo=pytz.UTC))
        t1 = self.ts.from_datetime(end_date.replace(tzinfo=pytz.UTC))

        # Sample every 1 minute for maximum precision
        time_step = 60.0 / 86400.0  # 1 minute in days
        num_samples = int((t1.tt - t0.tt) / time_step)

        print(f"Searching {num_samples} time points for transits...")
        print(f"Detection threshold: 3.0° angular separation")
        
        # Quick test of ISS positioning
        test_time = self.ts.now()
        try:
            test_iss_pos = self.iss_satellite - self.observer
            test_iss_relative = test_iss_pos.at(test_time)
            test_alt, test_az, _ = test_iss_relative.altaz()
            print(f"ISS test position: {test_alt.degrees:.1f}° altitude, {test_az.degrees:.1f}° azimuth")
        except Exception as e:
            print(f"ERROR in ISS positioning: {e}")

        # Track when we're close to a transit
        near_transit = False
        min_separation = float("inf")
        best_time = None

        for i in range(num_samples):
            t = self.ts.tt_jd(t0.tt + i * time_step)

            # Get positions relative to observer on Earth
            difference = self.iss_satellite - self.observer
            iss_relative = difference.at(t)
            
            moon_astrometric = (self.earth + self.observer).at(t).observe(self.moon)
            moon_relative = moon_astrometric.apparent()

            # Skip if ISS or Moon is below horizon  
            iss_alt, iss_az, _ = iss_relative.altaz()
            moon_alt, moon_az, _ = moon_relative.altaz()

            if iss_alt.degrees < 5 or moon_alt.degrees < 5:
                if near_transit and best_time is not None:
                    # We were tracking a close pass that has now ended
                    near_transit = False
                    if min_separation < 3.0:
                        transits.append(self._create_transit_event(best_time, min_separation))
                    min_separation = float("inf")
                    best_time = None
                continue

            # Calculate angular separation
            separation = self._angular_separation(iss_relative, moon_relative)
            moon_radius = self._moon_radius_degrees(t)

            # Check if we're within 5 degrees (close pass zone)
            if separation < 5.0:
                near_transit = True
                if separation < min_separation:
                    min_separation = separation
                    best_time = t
            elif near_transit and best_time is not None:
                # We've moved away from the close pass zone
                near_transit = False
                # Record all close passes under 3 degrees
                if min_separation < 3.0:
                    # Check if ISS is sunlit for actual transits
                    if min_separation < moon_radius or self._is_iss_sunlit(best_time):
                        transits.append(self._create_transit_event(best_time, min_separation))
                min_separation = float("inf")
                best_time = None

            # Progress indicator with debugging on June 21st
            if i % 1000 == 0 and i > 0:
                progress = i / num_samples * 100
                transits_found = len(transits)
                current_time = t.utc_datetime()
                if current_time.day == 21 and current_time.month == 6:
                    print(f"\nDEBUG {current_time.strftime('%Y-%m-%d %H:%M')}: ISS {iss_alt.degrees:.1f}° Moon {moon_alt.degrees:.1f}° Sep {separation:.2f}°")
                print(f"Progress: {progress:.1f}% - Found {transits_found} transit(s)", end="\r")

        print("\nSearch complete!")
        return transits

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
        iss_size_arcmin = 109 / iss_distance * 3437.75  # ISS is ~109m wide

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
            "iss_angular_size_arcmin": iss_size_arcmin,
            "moon_phase_degrees": phase_angle.degrees,
            "moon_illumination": (1 + np.cos(np.radians(phase_angle.degrees))) / 2 * 100,
        }


def calculate_iss_moon_transits(
    latitude: float,
    longitude: float,
    elevation: float,
    max_distance_km: Optional[float] = None,
) -> List[Dict[str, Any]]:
    """Calculate ISS-Moon transits for the next year from given location."""
    calculator = ISSMoonTransitCalculator(latitude, longitude, elevation)

    # Search for transits in the next 30 days
    start_date = datetime.now()
    end_date = start_date + timedelta(days=30)

    print(f"\nSearching for ISS-Moon transits from {latitude}°, {longitude}°, {elevation}m")
    print(f"Time range: {start_date.strftime('%Y-%m-%d')} to {end_date.strftime('%Y-%m-%d')}")

    if max_distance_km:
        print(f"Maximum travel distance: {max_distance_km} km")
        print("(Travel distance feature not yet implemented)")

    transits = calculator.find_transits(start_date, end_date, max_distance_km)

    return transits