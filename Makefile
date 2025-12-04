.PHONY: all build install clean docs docs-open test test-cov format lint help

help:
	@echo "SCRIBE Makefile"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build C++ module and install Python package (default)"
	@echo "  build     - Build only the C++ simulation engine (TICG)"
	@echo "  install   - Install only the Python package (pylib)"
	@echo "  clean     - Remove build artifacts"
	@echo "  docs      - Build Sphinx documentation"
	@echo "  docs-open - Build and open documentation in browser"
	@echo "  test      - Run unit tests"
	@echo "  test-cov  - Run tests with coverage report"
	@echo "  format    - Format code with black and ruff"
	@echo "  lint      - Check code style without modifying"
	@echo ""
	@echo "Quick start:"
	@echo "  conda env create -f environment.yml"
	@echo "  conda activate scribe"
	@echo "  make all"

all: build install

build: 
	(cd src && make pybind && mv pyticg* ../pylib)

install: 
	python setup.py bdist_wheel
	python -m pip install dist/* --force-reinstall

clean:
	rm -rf build dist *.egg-info pylib/*.so pylib/pyticg*

docs:
	cd docs && make html

docs-open: docs
	open docs/build/html/index.html

test:
	pytest tests/

test-cov:
	pytest tests/ --cov=pylib --cov-report=term-missing --cov-report=html

format:
	black pylib/ tests/ examples/
	ruff check --fix pylib/ tests/ examples/

lint:
	black --check pylib/ tests/ examples/
	ruff check pylib/ tests/ examples/
