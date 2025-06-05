# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python CLI application called "astro-stuff" for ISS-Moon transit calculations. It predicts when the International Space Station (ISS) will transit in front of the Moon from a given observer location. The project uses modern Python packaging with pyproject.toml and uv for dependency management.

## Development Commands

- **Install dependencies**: `uv sync` (creates/updates .venv with dependencies)
- **Install dev dependencies**: `uv sync --extra dev` (includes pytest, pytest-cov)
- **Run the CLI**: `uv run astro-stuff --latitude LAT --longitude LON --elevation ELEV`
- **Run tests**: `make test` or `uv run python -m pytest`
- **Format code**: `make format`
- **Lint code**: `make lint`
- **Clean cache**: `make clean`

## Architecture

- **Single command CLI**: `src/astro_stuff/__init__.py` contains all functionality in one file
- **CLI structure**: Uses `@click.command()` with a single focused command (no groups)
- **Project script**: Configured in pyproject.toml as `astro-stuff = "astro_stuff:main"`
- **Core algorithm**: Fast orbital intersection algorithm using ISS orbital mechanics
- **Caching**: TLE data cached in `~/.cache/astro-stuff/` for 7 days

## Algorithm Details

The ISS-Moon transit detection uses an intelligent orbital intersection algorithm:
- **Orbital sampling**: Samples ISS positions at 4 key points per 93-minute orbit
- **Moon trajectory prediction**: Pre-computes Moon positions every 3 hours with interpolation
- **Fast geometric screening**: Uses approximate angular distance before expensive calculations
- **Targeted refinement**: High-resolution search only around promising candidates
- **Performance**: ~40x faster than brute-force, completing in ~1 second

## Dependencies

### Core Dependencies
- **Click**: CLI framework
- **Skyfield**: Astronomical calculations and ephemeris data
- **NumPy**: Numerical computations and vector operations
- **Requests**: Fetching TLE data from Celestrak
- **PyTZ**: Timezone handling

### Development Dependencies
- **pytest**: Testing framework
- **pytest-cov**: Test coverage reporting
- **ruff**: Code formatting and linting

### System Requirements
- **Python**: Requires >= 3.12.1 (specified in .python-version)
- **Build system**: Uses Hatchling as the build backend

## Testing

The project includes comprehensive tests with mocking for deterministic results:
- Unit tests for CLI interface and parameter validation
- Tests for core calculation algorithms with mock data
- Performance validation and regression testing
- All tests use fixed test data to avoid time-dependent failures

## Performance Notes

- Initial run downloads ephemeris data (~15MB de421.bsp file)
- TLE data is cached and refreshed weekly
- Algorithm completes 30-day transit search in ~1 second
- Memory usage is minimal due to strategic sampling approach