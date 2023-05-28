.PHONY: install
install:
	pip install -e .[dev,docs]

.PHONY: test
test:
	pytest

.PHONY: clean
clean: clean-pyc clean-build

.PHONY: clean-pyc
clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f  {} +

.PHONY: clean-build
clean-build:
	rm -rf build/
	rm -rf dist/

.PHONY: build
build: clean-build
	python -m build --sdist --wheel .

.PHONY: docs
docs:
	cd docs && make html

.PHONY: publish
publish: build
	twine upload dist/*

.PHONY: publish-test
publish-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*
