test:
	python -m nose -x
debug:
	python -m nose --pdb
coverage:
	python -m nose -v --with-coverage --cover-package ledges --cover-html
