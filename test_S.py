import unittest
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
        sirius_rst = "%s/hf_%s.SIRIUS.RST" % (self.subdir, self.mol)
        self.cmo = SiriusRestart(name=sirius_rst).cmo
        aoproper = "%s/hf_%s.AOPROPER" % (self.subdir, self.mol)
        self.ls = prop.read(*labels, filename=aoproper)
        self.symorb = {2:(1,), 3:(1,), 5:(1,)}

    def tearDown(self):
        pass

    def test_setup(self):
        self.assertTupleEqual(tuple(self.cmo.nrow), (9, 6, 6, 2, 6, 2, 2, 1))

    def test_get_2p_orbitals(self):
        p_orbitals = get_orbitals(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        ref = numpy.zeros((34, 3))
        ref[9:15,  0] =  [1.0018, 0.0012, -0.0034, -0.0011, 0.0003, 0.]
        ref[15:21, 1] =  [1.0018, 0.0012, -0.0034, -0.0011, 0., 0.0003]
        ref[23:29, 2] =  [0.9966, -0.0048, 0.0065, 0.0022, 0.0004, 0.]
        
        numpy.testing.assert_allclose(p_orbitals, ref, atol=1e-4)

    def test_get_indices(self):
        orbital_indices = get_orbital_indices(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        self.assertTupleEqual(orbital_indices, (9, 15, 23))

    def test_get_spnorb0(self):
        ls = get_ls(self.cmo, self.symorb )
        np.testing.assert_allclose(
            ls[0], [[0,0,0], [0, 0, 0.03385742 ], [0, -0.03385742, 0]]
        )

    def test_get_spnorb1(self):
        ls = get_ls(self.cmo, self.symorb )
        np.testing.assert_allclose(
            ls[1], [[0, 0, -0.03385742 ],[0,0,0],  [0.03385742, 0, 0]]
        )

    def test_get_spnorb2(self):
        ls = get_ls(self.cmo, self.symorb )
        np.testing.assert_allclose(
            ls[2], [[0, 0.03385742 ], [0.03385742, 0, 0], [0,0,0]]
        )

            

if __name__ == "__main__":
    unittest.main()
