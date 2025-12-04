.PHONY: all build install clean docs docs-open test test-cov format lint help

help:
	@echo "SCRIBE Makefile"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build C++ module and install Python package (default)"
	@echo "  build     - Build only the C++ simulation engine (TICG)"
	@echo "  install   - Install only the Python package (scribe)"
	@echo "  clean     - Remove build artifacts"
	@echo "  docs      - Build Sphinx documentation"
	@echo "  docs-open - Build and open documentation in browser"
	@echo "  test      - Run fast unit tests (skips slow/integration tests)"
	@echo "  test-all  - Run all unit tests (including slow/integration tests)"
	@echo "  test-cov  - Run fast tests with coverage report"
	@echo "  format    - Format code with black and ruff"
	@echo "  lint      - Check code style without modifying"
	@echo ""
	@echo "Quick start:"
	@echo "  conda env create -f environment.yml"
	@echo "  conda activate scribe"
	@echo "  make all"

all: build install

build: 
	(cd src && make pybind && mv scribe_engine* ../scribe)

install: 
	python setup.py bdist_wheel
	python -m pip install dist/* --force-reinstall

clean:
	rm -rf build dist *.egg-info scribe/*.so scribe/scribe_engine*

docs:
	cd docs && make html

docs-open: docs
	open docs/build/html/index.html

test:
	pytest tests/ -m 'not slow'

test-all:
	pytest tests/

test-cov:
	pytest tests/ --cov=scribe --cov-report=term-missing --cov-report=html

format:
	black scribe/ tests/ examples/
	ruff check --fix scribe/ tests/ examples/

lint:
	black --check scribe/ tests/ examples/
	ruff check scribe/ tests/ examples/
