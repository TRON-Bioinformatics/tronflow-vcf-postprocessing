from unittest import TestCase
import pkg_resources
import subprocess
import os
from cyvcf2 import VCF


class IntegrationTests(TestCase):

    def setUp(self):
        self.test_vcf = pkg_resources.resource_filename(__name__, "resources/test.vcf")
        self.test2_vcf = pkg_resources.resource_filename(__name__, "resources/test2.vcf")
        self.input_files = pkg_resources.resource_filename(__name__, "resources/input_files.txt")
        self.sample1 = 'test'
        self.sample2 = 'test2'
        self._write_config_file(self.input_files, entries=[
            [self.sample1, self.test_vcf],
            [self.sample2, self.test2_vcf]
        ])

        self.test_vcf_with_duplicates = pkg_resources.resource_filename(__name__, "resources/test_with_duplicates.vcf")
        self.input_files_with_duplicates = pkg_resources.resource_filename(__name__, "resources/input_files_with_duplicates.txt")
        self.sample_with_duplicates = 'test_with_duplicates'
        self._write_config_file(self.input_files_with_duplicates, entries=[
            [self.sample_with_duplicates, self.test_vcf_with_duplicates]
        ])

        self.test_vcf_with_mnvs = pkg_resources.resource_filename(__name__, "resources/test_with_mnvs.vcf")
        self.input_files_with_mnvs = pkg_resources.resource_filename(__name__, "resources/input_files_with_mnvs.txt")
        self.sample_with_mnvs = 'test_with_mnvs'
        self._write_config_file(self.input_files_with_mnvs, entries=[
            [self.sample_with_mnvs, self.test_vcf_with_mnvs]
        ])

        self.module_settings_cmd = 'module purge; module load java/11.0.1;'
        workflow = pkg_resources.resource_filename(__name__, "../main.nf")
        self.nextflow_cmd = "nextflow {workflow} -profile test ".format(workflow=workflow)

    def test_help(self):
        self._run_nextflow(options='--help')

    def test_default_settings(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_default_settings")
        options = "--input_files {input_files} --output {output}".format(input_files=self.input_files, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample1, input_vcf=self.test_vcf)
        self._check_results(output_folder=output, sample=self.sample2, input_vcf=self.test2_vcf)

    def test_filter_pass(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_filter")
        options = "--input_files {input_files} --output {output} --filter PASS".format(
            input_files=self.input_files, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample1, input_vcf=self.test_vcf, filter=4)
        self._check_results(output_folder=output, sample=self.sample2, input_vcf=self.test2_vcf)

    def test_skip_split_vcf(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_skip_split_vcf")
        options = "--input_files {input_files} --output {output} --skip_split_vcf_by_type".format(
            input_files=self.input_files, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample1, input_vcf=self.test_vcf, split_by_type=False)
        self._check_results(output_folder=output, sample=self.sample2, input_vcf=self.test2_vcf, split_by_type=False)

    def test_skip_duplication_removal(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_skip_duplication_removal")
        options = "--input_files {input_files} --output {output} --skip_duplication_removal".format(
            input_files=self.input_files_with_duplicates, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_duplicates,
                            input_vcf=self.test_vcf_with_duplicates, duplicates=0)

    def test_skip_duplication_removal_disabled(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_skip_duplication_removal_disabled")
        options = "--input_files {input_files} --output {output}".format(
            input_files=self.input_files_with_duplicates, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_duplicates,
                            input_vcf=self.test_vcf_with_duplicates, duplicates=7)

    def test_skip_split_mnps(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_skip_split_mnps")
        options = "--input_files {input_files} --output {output} --skip_split_mnps".format(
            input_files=self.input_files_with_mnvs, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_mnvs,
                            input_vcf=self.test_vcf_with_mnvs, decomposed_from_mnvs=0)

    def test_skip_split_mnps_override_decompose_non_blocked_substitutions(self):
        output = pkg_resources.resource_filename(
            __name__, "resources/output_test_skip_split_mnps_override_decompose_non_blocked_substitutions")
        options = "--input_files {input_files} --output {output} --skip_split_mnps " \
              "--decompose_non_blocked_substitutions".format(
            input_files=self.input_files_with_mnvs, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_mnvs,
                            input_vcf=self.test_vcf_with_mnvs, decomposed_from_mnvs=0)

    def test_skip_split_mnps_disabled(self):
        output = pkg_resources.resource_filename(__name__, "resources/output_test_skip_split_mnps_disabled")
        options = "--input_files {input_files} --output {output}".format(
            input_files=self.input_files_with_mnvs, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_mnvs,
                            input_vcf=self.test_vcf_with_mnvs, decomposed_from_mnvs=5)

    def test_skip_split_mnps_disabled_and_decompose_non_blocked_substitutions(self):
        output = pkg_resources.resource_filename(
            __name__, "resources/output_test_skip_split_mnps_disabled_and_decompose_non_blocked_substitutions")
        options = "--input_files {input_files} --output {output} --decompose_non_blocked_substitutions".format(
            input_files=self.input_files_with_mnvs, output=output)
        self._run_nextflow(options=options)
        self._check_results(output_folder=output, sample=self.sample_with_mnvs,
                            input_vcf=self.test_vcf_with_mnvs, decomposed_from_mnvs=6)

    def _check_results(self, output_folder, sample, input_vcf, split_by_type=True, duplicates=0,
                       decomposed_from_mnvs=0, filter=0):
        observed_vcf = os.path.join(output_folder, "{sample}/{sample}.normalized.vcf".format(sample=sample))
        self.assertTrue(os.path.exists(observed_vcf))
        for variant_type in ['snps', 'indels', 'mnps', 'bnd', 'ref', 'other']:
            self.assertEqual(
                os.path.exists(os.path.join(output_folder, "{sample}/{sample}.normalized.{variant_type}.vcf".format(
                    sample=sample, variant_type=variant_type))),
                split_by_type, "Failed with VCF file of variant type {}".format(variant_type))
        self.assertTrue(os.path.exists(
            os.path.join(output_folder, "{sample}/{sample}.normalization.log".format(sample=sample))))
        self.assertTrue(os.path.exists(
            os.path.join(output_folder, "{sample}/{sample}_stats/counts_by_af.snps.dat".format(sample=sample))))
        expected_count_variants = len(list(VCF(input_vcf)))
        observed_count_variants = len(list(VCF(observed_vcf)))
        self.assertEqual(expected_count_variants, observed_count_variants + duplicates - decomposed_from_mnvs + filter,
                         'Unexpected difference in count of variants')

    def _run_nextflow(self, options):
        process = subprocess.Popen(self.module_settings_cmd + self.nextflow_cmd + options,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        print(stdout.decode("utf-8"))
        print(stderr.decode("utf-8"))
        self.assertEqual(process.returncode, 0)

    @staticmethod
    def _write_config_file(filename, entries):
        fd = open(filename, 'w')
        for e in entries:
            fd.write("\t".join(e) + "\n")
        fd.close()
