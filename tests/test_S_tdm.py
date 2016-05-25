import unittest
import os
import numpy
from scipy.constants import alpha
import vb.nod
import daltools
import util
from . import ledges
#from ledges import tdm

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

        self.spin_density = ledges.tdm.get_transition_spin_densities(*self.p)

    def tearDown(self):
        pass

    def test_v_aa(self):
        v_aa = numpy.ndarray((3, 3), dtype=numpy.complex64)
        for i in range(3):
            for j in range(3):
                v_aa[i, j] = self.so1[2] & self.spin_density.subblock[i][j]
        v_aa *= SO_FACTOR
        numpy.testing.assert_allclose(
            v_aa,
            [[ 0, 0.01692871, 0],[-0.01692871, 0, 0], [0, 0, 0]],
            rtol=1e-7, atol=1e-7
            )




if __name__ == "__main__":
    unittest.main()
