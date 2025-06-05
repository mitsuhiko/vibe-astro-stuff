.PHONY: test install dev clean format check-format lint

# Run tests
test:
	uv run pytest tests

# Install dependencies
install:
	uv sync

# Install with dev dependencies
dev:
	uv sync --extra dev

# Format code with ruff
format:
	uvx ruff format .

# Check code formatting
check-format:
	uvx ruff format --check .

# Lint code with ruff
lint:
	uvx ruff check .

# Clean up cache and build artifacts
clean:
	rm -rf .pytest_cache/
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
