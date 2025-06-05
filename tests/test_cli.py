from click.testing import CliRunner

from astro_stuff import main


class TestCLI:
    def test_hello_command(self):
        """Test the hello command works."""
        runner = CliRunner()
        result = runner.invoke(main, ["hello"])
        assert result.exit_code == 0
        assert "Hello from astro-stuff!" in result.output

    def test_iss_moon_transit_command_basic(self):
        """Test the iss-moon-transit command with required parameters."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "iss-moon-transit",
                "--latitude",
                "37.7749",
                "--longitude",
                "-122.4194",
                "--elevation",
                "16",
            ],
        )
        assert result.exit_code == 0
        assert (
            "Calculating ISS-moon transits for location: 37.7749°, -122.4194°, 16.0m"
            in result.output
        )
        assert "Transit calculations not yet implemented." in result.output

    def test_iss_moon_transit_command_with_max_distance(self):
        """Test the iss-moon-transit command with optional max distance parameter."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "iss-moon-transit",
                "--latitude",
                "37.7749",
                "--longitude",
                "-122.4194",
                "--elevation",
                "16",
                "--max-distance",
                "100",
            ],
        )
        assert result.exit_code == 0
        assert "Maximum travel distance: 100.0 km" in result.output

    def test_iss_moon_transit_command_missing_required_params(self):
        """Test the iss-moon-transit command fails with missing required parameters."""
        runner = CliRunner()
        result = runner.invoke(main, ["iss-moon-transit"])
        assert result.exit_code != 0
        assert "Missing option" in result.output

    def test_iss_moon_transit_command_short_flags(self):
        """Test the iss-moon-transit command with short flag versions."""
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "iss-moon-transit",
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
        assert (
            "Calculating ISS-moon transits for location: 40.7128°, -74.006°, 10.0m"
            in result.output
        )
        assert "Maximum travel distance: 50.0 km" in result.output
