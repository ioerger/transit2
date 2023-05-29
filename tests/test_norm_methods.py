import sys
import os

basedir = os.path.dirname(__file__)
sys.path.insert(0, basedir + '/../src/')

import shutil
import unittest
import os
import numpy

from transit_test import *

# fake setup for testing
from pytransit.globals import logging, gui
gui.is_active = False # normally checks sys.argv[] but tests use their own sys.argv

import pytransit.specific_tools.norm_tools as norm_tools
import pytransit.specific_tools.tnseq_tools as tnseq_tools
import pytransit.specific_tools.console_tools as console_tools

from pytransit.methods.gumbel import Method as GumbelMethod
from pytransit.methods.hmm        import Method as HMMMethod

from pytransit.methods.resampling import Method as ResamplingMethod

# fake setup for testing
from pytransit.globals import logging, gui
gui.is_active = False # normally checks sys.argv[] but tests use their own sys.argv

# RAW STATISTICS:
raw_ctrl_rep1 = (0.41855103545338784, 53.913866362844317, 128.81073464420675, 70.0, 3855.0, 4022244.0, 4.022213863266365, 33.02399374787363)
raw_ctrl_rep2 = (0.51571610481871188, 86.141009315729505, 167.03183885640027, 89.0, 5944.0, 6426550.0, 3.975522680364774, 33.45593924344892)

raw_exp_rep1 = (0.43946116212050129, 52.944722203605657, 120.47645336424084, 56.0, 47546.0, 3949941.0, 54.83729352080048, 4237.703504066973)
raw_exp_rep2 = (0.43891160109912203, 53.013082233094295, 120.78305084745763, 46.0, 217960.0, 3955041.0, 105.78284470780153, 14216.199240166581)
raw_exp_rep3 = (0.35898398230681589, 60.667260907445879, 168.99712493465759, 60.0, 102013.0, 4526081.0, 42.17620960361282, 2327.971405759911)

raw_means = [53.913866362844317, 86.141009315729505, 52.944722203605657, 53.013082233094295, 60.667260907445879]


class TestNormMethods(TransitTestCase):

    def test_nonorm(self):
        data,position = tnseq_tools.CombinedWig.gather_wig_data(all_data_list)
        norm_data,factors = norm_tools.normalize_data(data, "nonorm")
        self.assertTrue((factors == numpy.array([ 1.])).all())
        N = len(all_data_list)
        for k in range(N):
           self.assertEqual(numpy.mean(norm_data[k]), raw_means[k])


    def test_ttr(self):
        N = len(all_data_list)
        data,position = tnseq_tools.CombinedWig.gather_wig_data(all_data_list)
        norm_data,factors = norm_tools.normalize_data(data, "TTR")
        self.assertFalse((factors == numpy.ones(N)).all())
        for k in range(N):
           self.assertNotEqual(numpy.mean(norm_data[k]), raw_means[k])

    def test_resampling_nonorm(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "nonorm", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        pvals, qvals = significant_pvals_qvals(output)
        self.assertLessEqual(len(pvals), 5)
        self.assertLessEqual(len(qvals), 1)


    def test_resampling_ttr(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "TTR", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        pvals, qvals = significant_pvals_qvals(output)
        self.assertLessEqual(len(pvals), 1)
        self.assertLessEqual(len(qvals), 1)


    def test_resampling_nz_mean(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "nzmean", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        pvals, qvals = significant_pvals_qvals(output)
        self.assertLessEqual(len(pvals), 5)
        self.assertLessEqual(len(qvals), 1)

    def test_resampling_tot_reads(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "totreads", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        pvals, qvals = significant_pvals_qvals(output)
        self.assertLessEqual(len(pvals), 5)
        self.assertLessEqual(len(qvals), 1)

    def test_resampling_quantile(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "quantile", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        pvals, qvals = significant_pvals_qvals(output)
        self.assertLessEqual(len(pvals), 5)
        self.assertLessEqual(len(qvals), 1)

    def test_resampling_zinfnb(self):
        args = [ctrl_rep1, ctrl_rep2, small_annotation, output, "--s", "1000", "--n", "zinfnb", "-no-sr"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))

    """
    def test_resampling_bgc(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "--s", "1000", "--n", "betageom"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
        hits = count_hits(output)
        self.assertLessEqual(hits, 20)
    """


    """
    def test_resampling_a_bgc(self):
        args = [ctrl_data_txt, exp_data_txt, annotation, output, "--s", "1000", "--n", "aBGC"]
        G = ResamplingMethod.from_args(*console_tools.clean_args(args))
        self.assertTrue(os.path.exists(output))
    """
    
if __name__ == '__main__':
    unittest.main()