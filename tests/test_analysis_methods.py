import matplotlib
matplotlib.use('Agg')

import sys
import os

basedir = os.path.dirname(__file__)
sys.path.insert(0, basedir + '/../src/')

import shutil
import unittest

from transit_test import basedir, ctrl_rep1, ctrl_rep2, ctrl_data_txt, mini_wig, combined_wig, samples_metadata, samples_metadata_covariates, samples_metadata_interactions, KO_combined_wig, KO_samples_metadata, exp_rep1, exp_rep2, exp_rep3, exp_data_txt, all_data_list, annotation, small_annotation, output, hist_path, tpp_output_base, tpp_output_paths, reads1, test_multicontig, test_multicontig_reads1, test_multicontig_reads2, h37fna, TransitTestCase, count_hits, significant_pvals_qvals

# fake setup for testing
from pytransit.globals import logging, gui
gui.is_active = False # normally checks sys.argv[] but tests use their own sys.argv

from pytransit.specific_tools import norm_tools
from pytransit.specific_tools import tnseq_tools
from pytransit.specific_tools import console_tools
from pytransit.specific_tools.transit_tools import HAS_R

# Single condition methods
from pytransit.methods.anova  import Method as AnovaMethod
from pytransit.methods.gumbel import Method as GumbelMethod
from pytransit.methods.hmm    import Method as HMMMethod
from pytransit.methods.zinb   import Method as ZinbMethod

# Comparative methods
from pytransit.methods.resampling import Method as ResamplingMethod
from pytransit.methods.utest      import Method as UTestMethod

# Genetic Interactions
from pytransit.methods.gi import Method as GIMethod



class TestMethods(TransitTestCase):
    def test_resampling(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-no-sr"]
        method_object = None
        try:
            ResamplingMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
            import traceback
            traceback.print_exc()
            print(f'''error = {error}''')
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        self.assertLessEqual(
                abs(len(sig_pvals) - 37),
                2,
                "sig_pvals expected in range: %s, actual: %d" % ("[35, 39]", len(sig_qvals)))
        self.assertLessEqual(
                abs(len(sig_qvals) - 35),
                2,
                "sig_qvals expected in range: %s, actual: %d" % ("[33, 37]", len(sig_qvals))) # maybe acceptable range should be expanded to 38

    def test_resampling_combined_wig(self):
        # The conditions in the args should be matched case-insensitively.
        args = [ combined_wig, samples_metadata, small_annotation, "Glycerol", "Cholesterol", output, "-a"]
        try:
            method_object = ResamplingMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
            import traceback
            traceback.print_exc()
            print(f'''error = {error}''')
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        print(len(sig_pvals))
        print(len(sig_qvals))
        self.assertLessEqual(
                abs(len(sig_pvals) - 37),
                2,
                "sig_pvals expected in range: %s, actual: %d" % ("[35, 39]", len(sig_qvals)))
        self.assertLessEqual(
                abs(len(sig_qvals) - 35),
                1,
                "sig_qvals expected in range: %s, actual: %d" % ("[34, 36]", len(sig_qvals)))

    def test_resampling_adaptive(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "-a", "--ctrl_lib", "AA", "-no-sr", "--exp_lib", "AAA"]
        try:
            method_object = ResamplingMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
            import traceback
            traceback.print_exc()
            print(f'''error = {error}''')
        self.assertTrue(os.path.exists(output))
        (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-2, qcol=-1))
        self.assertLessEqual(
                abs(len(sig_pvals) - 37),
                2,
                "sig_pvals expected in range: %s, actual: %d" % ("[35, 39]", len(sig_qvals)))
        self.assertLessEqual(
                abs(len(sig_qvals) - 35),
                2,
                "sig_qvals expected in range: %s, actual: %d" % ("[34, 36]", len(sig_qvals)))

    def test_resampling_histogram(self):
        args = [ctrl_data_txt, exp_data_txt, small_annotation, output, "--s", "1000", "-h"]
        try:
            method_object = ResamplingMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
            import traceback
            traceback.print_exc()
            print(f'''error = {error}''')
        self.assertTrue(os.path.exists(output))
        self.assertTrue(
                os.path.isdir(hist_path),
                "histpath expected: %s" % (hist_path))

    def test_resampling_multistrain(self):
        args = [ctrl_data_txt, exp_data_txt, ','.join([small_annotation, small_annotation]), output, "-h"]
        try:
            method_object = ResamplingMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
            import traceback
            traceback.print_exc()
            print(f'''error = {error}''')
        self.assertTrue(os.path.exists(output))
        self.assertTrue(
                os.path.isdir(hist_path),
                "histpath expected: %s" % (hist_path))
    
    # @unittest.skipUnless(HAS_R, "requires R, rpy2")
    # def test_zinb(self):
    #     args = [combined_wig, samples_metadata, small_annotation, output, "--iN", "5.0", "--iC", "5.0" ]
    #     try:
    #         method_object = ZinbMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #         import traceback
    #         traceback.print_exc()
    #         print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))
    #     (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
    #     sig_qvals.sort()
    #     self.assertEqual(
    #         len(sig_pvals),
    #         31,
    #         "sig_pvals expected: %d, actual: %d" % (31, len(sig_pvals)))
    #     self.assertEqual(
    #         len(sig_qvals),
    #         30,
    #         "sig_qvals expected: %d, actual: %d" % (30, len(sig_qvals)))

    # @unittest.skipUnless(HAS_R, "requires R, rpy2")
    # def test_zinb_covariates(self):
    #     args = [combined_wig, samples_metadata_covariates, small_annotation, output, "--covars", "batch", "--group-by", "NewConditionCol", "--iN", "5.0", "--iC", "5.0", ]
    #     try:
    #         method_object = ZinbMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #         import traceback
    #         traceback.print_exc()
    #         print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))
    #     (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
    #     sig_qvals.sort()
    #     self.assertEqual(len(sig_pvals), 15, "sig_pvals expected: %d, actual: %d" % (15, len(sig_pvals)))
    #     self.assertEqual(len(sig_qvals), 10, "sig_qvals expected: %d, actual: %d" % (10, len(sig_qvals)))

    # @unittest.skipUnless(HAS_R, "requires R, rpy2")
    # def test_zinb_interactions(self):
    #     args = [combined_wig, samples_metadata_interactions, small_annotation,   output, "--covars", "batch", "--interactions", "atm", "--iN", "5.0", "--iC", "5.0", ]
    #     try:
    #         method_object = ZinbMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #         import traceback
    #         traceback.print_exc()
    #         print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))
    #     (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
    #     sig_qvals.sort()
    #     self.assertEqual(
    #         len(sig_pvals),
    #         3,
    #         "sig_pvals expected: %d, actual: %d" % (3, len(sig_pvals)))
    #     self.assertEqual(
    #         len(sig_qvals),
    #         0,
    #         "sig_qvals expected: %d, actual: %d" % (0, len(sig_qvals)))
    
    # def test_anova(self):
    #     args = [combined_wig, samples_metadata, small_annotation,  output]
    #     try:
    #         method_object = AnovaMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #         import traceback
    #         traceback.print_exc()
    #         print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))
    #     (sig_pvals, sig_qvals) = (significant_pvals_qvals(output, pcol=-3, qcol=-2))
    #     sig_qvals.sort()
    #     self.assertEqual(
    #         len(sig_pvals),
    #         30,
    #         "sig_pvals expected: %d, actual: %d" % (30, len(sig_pvals)))
    #     self.assertEqual(
    #         len(sig_qvals),
    #         24,
    #         "sig_qvals expected: %d, actual: %d" % (24, len(sig_qvals)))
    
    def test_gumbel(self):
        args = [ combined_wig, samples_metadata, small_annotation, "Cholesterol", output, "--s", "1000", "--b", "100"]
        try:
            method_object = GumbelMethod.from_args(*console_tools.clean_args(args))
        except Exception as error:
                import traceback
                traceback.print_exc()
                print(f'''error = {error}''')
        self.assertTrue(os.path.exists(output))

    # def test_hmm(self):
    #     args = [mini_wig, small_annotation, output]
    #     try:
    #         method_object = HMMMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #             import traceback
    #             traceback.print_exc()
    #             print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))
    #     genes_path = output.rsplit(".", 1)[0] + "_genes." + output.rsplit(".", 1)[1]
    #     self.assertTrue(os.path.exists(genes_path))


    # def test_utest(self):
    #     args = [ combined_wig, samples_metadata, small_annotation,  "Glycerol", "Cholesterol", output, ]
    #     try:
    #         method_object = UTestMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #             import traceback
    #             traceback.print_exc()
    #             print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))

    # def test_gi(self):
    #     #  usage: {console_tools.subcommand_prefix} gi <combined_wig> <samples_metadata> <conditionA1> <conditionB1> <conditionA2> <conditionB2> <prot_table> <output_file> [optional arguments]
    #     args = [KO_combined_wig, KO_samples_metadata, small_annotation, "H37Rv_day0","H37Rv_day32","Rv2680_day0","Rv2680_day32", output]
    #     try:
    #         method_object = GIMethod.from_args(*console_tools.clean_args(args))
    #     except Exception as error:
    #             import traceback
    #             traceback.print_exc()
    #             print(f'''error = {error}''')
    #     self.assertTrue(os.path.exists(output))

if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestMethods)
    #unittest.TextTestRunner(verbosity=2).run(suite)
