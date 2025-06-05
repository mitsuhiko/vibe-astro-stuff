# Fixed test data for deterministic testing

# Sample TLE data for ISS (fixed point in time for consistent testing)
SAMPLE_ISS_TLE = {
    "line1": "1 25544U 98067A   24173.12345678  .00001234  00000-0  12345-4 0  9999",
    "line2": "2 25544  51.6400 123.4567  0001234  12.3456 347.8901 15.48900000123456",
}

# Test location: Kansas (matching user's previous test data)
TEST_LOCATION = {
    "latitude": 39.7392,
    "longitude": -104.9903,
    "elevation": 1609,  # Denver elevation in meters
}

# Mock transit event for output testing
MOCK_TRANSIT_EVENT = {
    "time": "2024-06-21T12:30:45.123456",
    "separation_degrees": 0.25,
    "moon_radius_degrees": 0.26,
    "iss_altitude": 45.7,
    "iss_azimuth": 180.5,
    "moon_altitude": 50.2,
    "moon_azimuth": 175.3,
    "iss_distance_km": 408.2,
    "iss_angular_size_arcsec": 20.5,
    "moon_phase_degrees": 90.0,
    "moon_illumination": 75.3,
}
