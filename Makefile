
SHELL := /bin/bash

all : clean test

.PHONY: all test clean

clean:
	rm -rf output
	rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash test/scripts/run_test_0.sh
	bash test/scripts/run_test_1.sh
	bash test/scripts/run_test_2.sh
	bash test/scripts/run_test_3.sh
	bash test/scripts/run_test_4.sh
	bash test/scripts/run_test_5.sh
	bash test/scripts/run_test_6.sh
	bash test/scripts/run_test_7.sh
	bash test/scripts/run_test_8.sh
	bash test/scripts/run_test_10.sh
	bash test/scripts/run_test_11.sh
	bash test/scripts/run_test_12.sh
	bash test/scripts/run_test_13.sh
	bash test/scripts/run_test_14.sh
