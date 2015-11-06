.PHONY: flake8 test coverage doc
export PATH := bin:$(PATH)

flake8:
	flake8 --max-line-length=100 --count --statistics --exit-zero kt_simul/

test:
	nosetests kt_simul/mecabio -v

coverage:
	nosetests kt_simul/mecabio --with-coverage --cover-package=kt_simul.mecabio -v

clean:
	find . -name "*.pyc" -exec rm -rf {} \;
	find . -name "__pycache__" -exec rm -rf {} \;

doc:
	cd doc/ && pandoc simu_spindle_1d.md -o simu_spindle_1d.pdf
	cd doc/ && pandoc simu_spindle_3d.md -o simu_spindle_3d.pdf
	cd doc/ && pandoc simu_spindle_three_attachments.md -o simu_spindle_three_attachments.pdf
