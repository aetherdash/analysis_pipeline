.PHONY: all
all: clean init build

clean:
	git clean -Xfd

init:
	pip install --upgrade pipenv
	pipenv install --dev

test:
	pipenv run py.test -s

build:
	pipenv run python setup.py bdist_wheel

deploy:
	pipenv run python setup.py bdist_wheel upload -r http://pypi.aetherbio.com/

.PHONY: static_analysis
static_analysis: flake8 isort mypy pycodestyle pylint coverage

black:
	pipenv run black federation tests setup.py

flake8:
	pipenv run flake8

isort:
	pipenv run isort --apply

mypy:
	pipenv run mypy federation tests

pycodestyle:
	pipenv run pycodestyle federation tests setup.py

pylint:
	pipenv run pylint federation setup.py tests

coverage:
	pipenv run py.test --cov-config .coveragerc --verbose --cov-report term --cov-report xml --cov-report html --cov-fail-under=95 --cov=federation tests
