.PHONY: flake8 test coverage

flake8:
	flake8 --max-line-length=100 --count --statistics --exit-zero kt_simul/

test:
	nosetests kt_simul/mecabio -v

coverage:
	nosetests kt_simul/mecabio --with-coverage --cover-package=kt_simul.mecabio -v

clean:
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "__pycache__" -exec rm -rf {} \;
