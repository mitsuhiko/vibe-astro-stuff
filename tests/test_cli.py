from datetime import datetime
from unittest.mock import patch, MagicMock
import pytest
from click.testing import CliRunner

from astro_stuff import main, ISSMoonTransitCalculator


class TestCLI:
    def test_missing_required_params(self):
        """Test the command fails with missing required parameters."""
        runner = CliRunner()
        result = runner.invoke(main, [])
        assert result.exit_code != 0
        assert "Missing option" in result.output

    def test_short_flags(self):
        """Test the command with short flag versions."""
        runner = CliRunner()
        with patch("astro_stuff.calculate_iss_moon_transits") as mock_calc:
            mock_calc.return_value = []
            result = runner.invoke(
                main,
                [
                    "-lat",
                    "40.7128",
                    "-lon",
                    "-74.0060",
                    "-elev",
                    "10",
                    "-d",
                    "50",
                ],
            )
            assert result.exit_code == 0
            mock_calc.assert_called_once_with(40.7128, -74.006, 10.0, 50.0, False)

    def test_verbose_flag(self):
        """Test the verbose flag is passed through."""
        runner = CliRunner()
        with patch("astro_stuff.calculate_iss_moon_transits") as mock_calc:
            mock_calc.return_value = []
            result = runner.invoke(
                main,
                [
                    "--latitude",
                    "37.7749",
                    "--longitude",
                    "-122.4194",
                    "--elevation",
                    "16",
                    "--verbose",
                ],
            )
            assert result.exit_code == 0
            mock_calc.assert_called_once_with(37.7749, -122.4194, 16.0, None, True)

    def test_no_transits_found(self):
        """Test output when no transits are found."""
        runner = CliRunner()
        with patch("astro_stuff.calculate_iss_moon_transits") as mock_calc:
            mock_calc.return_value = []
            result = runner.invoke(
                main,
                [
                    "--latitude",
                    "37.7749",
                    "--longitude",
                    "-122.4194",
                    "--elevation",
                    "16",
                ],
            )
            assert result.exit_code == 0
            assert "No ISS-Moon transits or close passes found" in result.output

    def test_transit_output_formatting(self):
        """Test the output formatting for transit events."""
        runner = CliRunner()
        mock_transit = {
            "time": datetime(2024, 6, 21, 12, 30, 45, 123456),
            "separation_degrees": 0.25,
            "moon_radius_degrees": 0.26,
            "iss_angular_size_arcsec": 20.5,
            "iss_distance_km": 408.2,
            "iss_azimuth": 180.5,
            "iss_altitude": 45.7,
            "moon_illumination": 75.3,
        }

        with patch("astro_stuff.calculate_iss_moon_transits") as mock_calc:
            mock_calc.return_value = [mock_transit]
            result = runner.invoke(
                main,
                [
                    "--latitude",
                    "37.7749",
                    "--longitude",
                    "-122.4194",
                    "--elevation",
                    "16",
                ],
            )
            assert result.exit_code == 0
            assert "Found 1 ISS-Moon transit(s)" in result.output
            assert "Transit #1:" in result.output
            assert "Fri 2024-06-21, 12:30:45.12" in result.output
            assert "20.50″" in result.output
            assert "408.20 km" in result.output
            assert "75.3%" in result.output


class TestISSMoonTransitCalculator:
    @pytest.fixture
    def mock_skyfield_data(self):
        """Mock Skyfield components for testing."""
        with patch("astro_stuff.Loader") as mock_loader:
            mock_ts = MagicMock()
            mock_ephemeris = MagicMock()
            mock_earth = MagicMock()
            mock_moon = MagicMock()

            mock_loader_instance = MagicMock()
            mock_loader.return_value = mock_loader_instance
            mock_loader_instance.timescale.return_value = mock_ts
            mock_loader_instance.return_value = mock_ephemeris
            mock_ephemeris.__getitem__ = lambda self, key: {
                "earth": mock_earth,
                "moon": mock_moon,
            }[key]

            yield {
                "loader": mock_loader,
                "ts": mock_ts,
                "ephemeris": mock_ephemeris,
                "earth": mock_earth,
                "moon": mock_moon,
            }

    def test_calculator_initialization(self, mock_skyfield_data):
        """Test calculator initializes with correct parameters."""
        with patch("astro_stuff.ISSMoonTransitCalculator._load_iss_tle"):
            calc = ISSMoonTransitCalculator(39.7392, -104.9903, 1609)
            assert calc.latitude == 39.7392
            assert calc.longitude == -104.9903
            assert calc.elevation == 1609

    def test_angular_separation_calculation(self, mock_skyfield_data):
        """Test angular separation calculation between two positions."""
        with patch("astro_stuff.ISSMoonTransitCalculator._load_iss_tle"):
            calc = ISSMoonTransitCalculator(39.7392, -104.9903, 1609)

            # Mock two positions with known angular separation
            pos1 = MagicMock()
            pos2 = MagicMock()
            pos1.position.au = [1.0, 0.0, 0.0]  # Unit vector along x-axis
            pos2.position.au = [0.0, 1.0, 0.0]  # Unit vector along y-axis

            # These should be 90 degrees apart
            separation = calc._angular_separation(pos1, pos2)
            assert abs(separation - 90.0) < 0.001

    def test_moon_radius_calculation(self, mock_skyfield_data):
        """Test moon radius calculation."""
        with patch("astro_stuff.ISSMoonTransitCalculator._load_iss_tle"):
            calc = ISSMoonTransitCalculator(39.7392, -104.9903, 1609)

            # Mock time and distance
            mock_time = MagicMock()
            mock_distance = MagicMock()
            mock_distance.km = 384400  # Average Earth-Moon distance in km

            with patch.object(calc, "moon") as mock_moon, patch.object(calc, "earth"):
                mock_moon_earth = MagicMock()
                mock_moon.__sub__.return_value = mock_moon_earth
                mock_moon_earth.at.return_value.distance.return_value = mock_distance

                radius = calc._moon_radius_degrees(mock_time)
                # Moon's angular radius should be approximately 0.26 degrees
                assert 0.25 < radius < 0.27
