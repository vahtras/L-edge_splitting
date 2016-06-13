test:
	python -m nose -x
pytest:
	python -m pytest -x tests
debug:
	python -m nose -vx --pdb
coverage:
	python -m nose -v --with-coverage --cover-package ledges --cover-html
