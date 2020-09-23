.PHONY: install clean-pyc clean-build build test publish docs-init docs

install:
	pip install -r requirements-dev.txt
	pip install -e .

test:
	pytest

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f  {} +

clean-build:
	rm -rf build/
	rm -rf dist/

clean: clean-pyc clean-build

build: clean-build
	python setup.py sdist
	python setup.py bdist_wheel

docs-init:
	pip install -r docs/requirements_doc.txt

docs:
	cd docs && make html

publish: build
	twine upload dist/*

publish-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*
