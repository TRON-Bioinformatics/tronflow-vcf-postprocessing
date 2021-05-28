
all : clean test check

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
	nextflow main.nf -profile test,conda --output output/test1 --input_files test_data/test_input.txt

check:
	test -s output/test1/sample1/sample1.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample2/sample2.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample3/sample4.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test1/sample4/sample4.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
