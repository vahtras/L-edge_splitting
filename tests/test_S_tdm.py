import unittest
import os
from scipy.constants import alpha
import vb.nod
import daltools
import util

SO_FACTOR = alpha**2/2

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
            C=cmo
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

        nb = self.so1[0].shape[0]
        self.spin_density = util.subblocked.matrix([nb, nb, nb], [nb, nb, nb])
        for i, p_i in enumerate(self.p):
            for j, p_j in enumerate(self.p):
                da_ij, db_ij = vb.nod.ao_transition_matrix(p_i, p_j)
                self.spin_density.subblock[i][j] = (da_ij - db_ij)/2

    def tearDown(self):
        pass

    def test_v_xaxa(self):
        """verify alpha-alpha block"""
        d_xx = vb.nod.Dao(self.p_x, self.p_x)
        v_xaxa = SO_FACTOR*self.so1[2] & self.spin_density.subblock[0][0]
        self.assertAlmostEqual(v_xaxa, 0)

    def test_v_xaya(self):
        """verify alpha-alpha block"""
        d_xy = vb.nod.ao_transition_matrix(self.p_x, self.p_y)
        v_xaya = SO_FACTOR*self.so1[2] & self.spin_density.subblock[0][1]
        self.assertAlmostEqual(v_xaya, 0.01692871)

    def test_v_xaza(self):
        """verify alpha-alpha block"""
        d_xz = vb.nod.ao_transition_matrix(self.p_x, self.p_z)
        v_xaza = SO_FACTOR*self.so1[2] & self.spin_density.subblock[0][2]
        self.assertAlmostEqual(v_xaza, 0.0)

    def test_v_yaxa(self):
        """verify alpha-alpha block"""
        d_yx = vb.nod.ao_transition_matrix(self.p_y, self.p_x)
        v_yaxa = SO_FACTOR*self.so1[2] & self.spin_density.subblock[1][0]
        self.assertAlmostEqual(v_yaxa, -0.01692871)

    def test_v_yaya(self):
        """verify alpha-alpha block"""
        d_yy = vb.nod.ao_transition_matrix(self.p_y, self.p_y)
        v_yaya = SO_FACTOR*self.so1[2] & self.spin_density.subblock[1][1]
        self.assertAlmostEqual(v_yaya, 0.0)

    def test_v_yaza(self):
        """verify alpha-alpha block"""
        d_yz = vb.nod.ao_transition_matrix(self.p_y, self.p_z)
        v_yaza = SO_FACTOR*self.so1[2] & self.spin_density.subblock[1][2]
        self.assertAlmostEqual(v_yaza, 0.0)

    def test_v_zaxa(self):
        """verify alpha-alpha block"""
        d_zx = vb.nod.ao_transition_matrix(self.p_z, self.p_x)
        v_zaxa = SO_FACTOR*self.so1[2] & self.spin_density.subblock[2][0]
        self.assertAlmostEqual(v_zaxa, 0.0)

    def test_v_zaya(self):
        """verify alpha-alpha block"""
        d_zy = vb.nod.ao_transition_matrix(self.p_z, self.p_y)
        v_zaya = SO_FACTOR*self.so1[2] & self.spin_density.subblock[2][1]
        self.assertAlmostEqual(v_zaya, 0.0)

    def test_v_zaza(self):
        """verify alpha-alpha block"""
        d_zz = vb.nod.ao_transition_matrix(self.p_z, self.p_z)
        v_zaza = SO_FACTOR*self.so1[2] & self.spin_density.subblock[2][2]
        self.assertAlmostEqual(v_zaza, 0.0)



if __name__ == "__main__":
    unittest.main()
