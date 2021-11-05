
SHELL := /bin/bash

all : clean test

clean:
	rm -rf output
	#rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash tests/run_test_0.sh
	bash tests/run_test_1.sh
	bash tests/run_test_2.sh
	bash tests/run_test_3.sh
	bash tests/run_test_4.sh
	bash tests/run_test_5.sh
