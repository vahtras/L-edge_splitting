import unittest
import numpy

from daltools.sirrst import SiriusRestart
from daltools import prop

from ledges.two_p_hole import *

labels = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')
TOL = (1e-7, 1e-7)

class S_NosymmetryTest(unittest.TestCase):

    def setUp(self):
        self.subdir = 'tests/S_nosym'
        self.mol = 'S'
        self.dal_tar_gz = "%s/hf_%s.tar.gz" % (self.subdir, self.mol)
        sirius_rst = "%s/hf_%s.SIRIUS.RST" % (self.subdir, self.mol)
        self.cmo = SiriusRestart(name=sirius_rst).cmo
        self.aoproper = "%s/hf_%s.AOPROPER" % (self.subdir, self.mol)
        self.lsmo = get_ls1(self.cmo, {1:(3, 4, 5)}, self.aoproper)

    def tearDown(self):
        pass

    def test_setup(self):
        self.assertTupleEqual(tuple(self.cmo.nrow), (34,))

    def test_get_2p_orbitals(self):
        p_orbitals = get_orbitals(self.cmo, {1: (3, 4, 5)})
        ref = numpy.loadtxt('tests/S_nosym/p345.txt')
        
        numpy.testing.assert_allclose(p_orbitals, ref, atol=1e-4)

    def test_get_indices(self):
        orbital_indices = get_orbital_indices(self.cmo, {1:(3, 4, 5,)})
        self.assertTupleEqual(orbital_indices, (2, 3, 4))

    def test_get_x1spnorb(self):
        numpy.testing.assert_allclose(
            self.lsmo[0],
            [[0,0,-0.00000050], [0, 0, 0.03396038], [0.00000050, -0.03396038, 0]],
            *TOL
        )

    def test_get_y1spnorb(self):
        numpy.testing.assert_allclose(
            self.lsmo[1],
            [[0, -0.00000001, -0.03385742 ], [0.00000001, 0, -0.00000050], [0.03385742, 0.00000050, 0]],
            *TOL
        )

    def test_get_z1spnorb(self):
        numpy.testing.assert_allclose(
            self.lsmo[2], [[0, 0.03385742, 0], [-0.03385742,  0, 0], [0, 0, 0]],
            *TOL
        )
