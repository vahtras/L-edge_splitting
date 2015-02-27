import unittest
import mock
import sys
import os
import numpy
from two_p_hole import *

sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
from dalmisc.sirrst import SiriusRestart
from daltools import prop

labels = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

class S_Test(unittest.TestCase):

    def setUp(self):
        self.subdir = 'S'
        self.mol = 'S'
        self.dal_tar_gz = "%s/hf_%s.tar.gz" % (self.subdir, self.mol)
        self.symorb = {2:(1,), 3:(1,), 5:(1,)}
        sys.argv = ['S_Test']


    def tearDown(self):
        sys.argv = ['S_Test']


    def test_all(self):
        er = two_p_eigenvalues(self.dal_tar_gz, self.symorb)
        ref_er = [-0.03389175,  -0.03389175,  0.01691156, 0.01691156, 0.01698019,  0.01698019]
        numpy.testing.assert_allclose(er, ref_er, rtol=1e-6)


    def test_main(self):
        sys.argv.append("{2: (1,), 3: (1,), 5: (1,)}")
        sys.argv.append(self.dal_tar_gz)
        print sys.argv
        main()

class SNosymmetryTest(unittest.TestCase):

    def setUp(self):
        self.subdir = 'S_nosym'
        self.mol = 'S'
        self.dal_tar_gz = "%s/hf_%s.tar.gz" % (self.subdir, self.mol)
        self.symorb = {1: (3, 4, 5)}
        sys.argv = ['SNosymeetrytest']


    def tearDown(self):
        sys.argv = ['SNosymeetrytest']


    def test_all(self):
        er = two_p_eigenvalues(self.dal_tar_gz, self.symorb)
        ref_er = [-0.03389175,  -0.03389175,  0.01691156, 0.01691156, 0.01698019,  0.01698019]
        numpy.testing.assert_allclose(er, ref_er, rtol=1e-6)


    def test_main(self):
        sys.argv.append("{1: (3, 4, 5)}")
        sys.argv.append(self.dal_tar_gz)
        print sys.argv
        main()

if __name__ == "__main__":
    unittest.main()
