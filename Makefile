
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
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --output output/test1
	nextflow main.nf -profile test,conda --output output/test2 --filter PASS,MNV
	nextflow main.nf -profile test,conda --output output/test3 --skip_decompose_complex
	nextflow main.nf -profile test,conda --input_files test_data/test_input_no_ad.txt --output output/test4
	nextflow main.nf -profile test,conda --output output/test5 --input_files false --input_vcf test_data/test_single_sample.vcf

check:
	test -s output/test1/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test3/tumor_normal/tumor_normal.normalized.vcf || { echo "Missing test 3 output file!"; exit 1; }
	test -s output/test1/single_sample/single_sample.normalized.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/single_sample/single_sample.normalized.vcf || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test3/single_sample/single_sample.normalized.vcf || { echo "Missing test 3 output file!"; exit 1; }
	test -s output/test4/sample_no_ad/sample_no_ad.normalized.vcf || { echo "Missing test 4 output file!"; exit 1; }
	test -s output/test5/test_single_sample/test_single_sample.normalized.vcf || { echo "Missing test 5 output file!"; exit 1; }