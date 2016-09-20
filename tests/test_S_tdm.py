"""Test general TDM version using reference values from the MO version"""
import unittest
import os
import numpy
from scipy.constants import alpha
import vb.nod
import daltools
import util
from util.full import matrix
from . import ledges
#from ledges import tdm
from ledges.two_p_hole import makeV2

SO_FACTOR = alpha**2/2
TOL = {'rtol': 1e-7, 'atol': 1e-7}

class TDMTest(unittest.TestCase):
    """Repreduce results for S_nosymm calculation with VB method"""

    def setUp(self):
        self.subdir = 'tests/S_nosym'
        fpath = lambda f: os.path.join(self.subdir, f)
        vb.nod.Nod.S = daltools.one.read(
            'OVERLAP',
            fpath('hf_S.AOONEINT')
            ).unpack().unblock()
        rst = daltools.sirrst.SiriusRestart(fpath("hf_S.SIRIUS.RST"))
        cmo = rst.cmo.unblock()
        self.p_x = vb.nod.Nod(
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
            [0, 1, 3, 4, 5, 6, 7, 8],
            C=cmo
            )
        self.p_y = vb.nod.Nod(
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
            [0, 1, 2, 4, 5, 6, 7, 8],
            C=-cmo # The change of sign is from odd permutations to reach the
                   # the hole orbital 3, compared to the mo implementation
            )
        self.p_z = vb.nod.Nod(
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
            [0, 1, 2, 3, 5, 6, 7, 8],
            C=cmo
            )
        self.p = (self.p_x, self.p_y, self.p_z)
        self.so1 = daltools.prop.read(
            'X1SPNORB', 'Y1SPNORB', 'Z1SPNORB',
            filename=fpath("hf_S.AOPROPER")
            )
        for hso in self.so1:
            hso *= SO_FACTOR

        self.spin_density = ledges.tdm.get_transition_spin_densities(*self.p)
        self.lx = numpy.array(
            [[0,0,-0.00000050], [0, 0, 0.03396038], [0.00000050, -0.03396038, 0]]
            )
        self.ly = numpy.array(
            [[0, -0.00000001, -0.03385742 ], [0.00000001, 0, -0.00000050], [0.03385742, 0.00000050, 0]]
            )
        self.lz = numpy.array([[0, 0.03385742, 0], [-0.03385742,  0, 0], [0, 0, 0]])

    def tearDown(self):
        pass

    def test_v_aa(self):
        v_aa = makeV2(self.so1, self.spin_density)[:3, :3]
        numpy.testing.assert_allclose(
            v_aa, 1j/2*self.lz, **TOL
            )

    def test_v_bb(self):
        v_bb = makeV2(self.so1, self.spin_density)[3:, 3:]
        numpy.testing.assert_allclose(
            v_bb, -1j/2*self.lz, **TOL
            )

    def test_v_ab(self):
        v_ab = makeV2(self.so1, self.spin_density)[:3, 3:]
        v_ab_ref = 1j/2*(self.lx + 1j*self.ly)
        numpy.testing.assert_allclose(v_ab, v_ab_ref, **TOL)

    def test_v_ba(self):
        v_ba = makeV2(self.so1, self.spin_density)[3:, :3]
        v_ba_ref = 1j/2*(self.lx - 1j*self.ly)
        numpy.testing.assert_allclose(v_ba, v_ba_ref, **TOL)


    def test_state_overlap(self):
        state_overlap = ledges.tdm.get_state_overlap(*self.p)
        numpy.testing.assert_allclose(state_overlap, numpy.eye(3), **TOL)




if __name__ == "__main__":
    unittest.main()
