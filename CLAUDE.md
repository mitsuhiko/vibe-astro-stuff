# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python CLI application called "astro-stuff" for celestial calculations, built using the Click framework. The project uses modern Python packaging with pyproject.toml and uv for dependency management.

## Development Commands

- **Install dependencies**: `uv sync` (creates/updates .venv with dependencies)
- **Run the CLI**: `uv run astro-stuff`

## Architecture

- **Entry point**: `src/astro_stuff/__init__.py` contains the main Click group and commands
- **CLI structure**: Uses Click's group/command pattern with `@click.group()` as the main entry point
- **Project script**: Configured in pyproject.toml as `astro-stuff = "astro_stuff:main"`

## Dependencies

- **Click**: Primary framework for CLI interface
- **Python**: Requires >= 3.12.1 (specified in .python-version)
- **Build system**: Uses Hatchling as the build backend