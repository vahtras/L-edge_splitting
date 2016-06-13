import unittest
import mock
import sys
import os
import numpy
import tarfile
import tempfile
from ledges.two_p_hole import *

from daltools.sirrst import SiriusRestart
from daltools import prop

labels = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

class ClTest(unittest.TestCase):

    def setUp(self):
        self.subdir = 'tests/Cl'
        self.mol = 'Cl'
        self.dal_tar_gz = "%s/lr_%s.tar.gz" % (self.subdir, self.mol)
        tgz = tarfile.open(self.dal_tar_gz)
        tmpdir = tempfile.gettempdir()
        tgz.extractall(path=tmpdir) 
        self.sirius_rst = SiriusRestart(
            name=os.path.join(tmpdir, 'SIRIUS.RST')
            )
        self.cmo = self.sirius_rst.cmo
        self.aoproper = os.path.join(tmpdir, 'AOPROPER')
        self.ao2soint = os.path.join(tmpdir, 'AO2SOINT')
        self.ls = prop.read(*labels, filename=self.aoproper)
        self.symorb = {2:(2,), 3:(2,), 5:(2,)}

    def tearDown(self):
        pass

    def test_setup(self):
        self.assertTupleEqual(tuple(self.cmo.nrow), (9, 4, 4, 2, 4, 2, 2, 0))

    def test_get_2p_orbitals(self):
        p_orbitals = get_orbitals(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        ref = numpy.zeros((27, 3))
        ref[9: 13, 0] = [1.00005991, 0.00000495, -0.00000734, -0.00006580 ]
        ref[13: 17, 1] = [1.00005991, 0.00000495, -0.00000734, -0.00006580 ]
        ref[19: 23, 2] = [1.00005991, 0.00000495, -0.00000734, -0.00006580 ]
        
        numpy.testing.assert_allclose(p_orbitals, ref, atol=1e-3)

    def test_get_indices(self):
        orbital_indices = get_orbital_indices(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        self.assertTupleEqual(orbital_indices, (9, 13, 19))

    def test_get_x1spnorb(self):
        ls = get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[0], [[0,0,0], [0, 0, 0.0030434444], [0, -0.0030434444, 0]]
        )

    def test_get_y1spnorb(self):
        ls = get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[1], [[0, 0, -0.0030434444],[0,0,0],  [0.0030434444, 0, 0]]
        )

    def test_get_z1spnorb(self):
        ls = get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[2], [[0, 0.0030434444, 0], [-0.0030434444, 0, 0], [0,0,0]],
            rtol=1e-6
        )

    def test_get_x2spnorb(self):
        ls2 = get_ls2(self.sirius_rst, self.symorb, self.ao2soint,
            reference_occupation=(
                (2, 2, 2), (2, 2), (2, 2), (), (2, 2), (), (), ()
                )
            )
        numpy.testing.assert_allclose(
            ls2[0], [[0,0,0], [0, 0, -0.000562559 ], [0, 0.000562559, 0]],
            rtol=1e-6
        )

    def test_get_y2spnorb(self):
        ls2 = get_ls2(self.sirius_rst, self.symorb, self.ao2soint,
            reference_occupation=(
                (2, 2, 2), (2, 2), (2, 2), (), (2, 2), (), (), ()
                )
            )
        numpy.testing.assert_allclose(
            ls2[1], [[0, 0, 0.000562559], [0,0,0], [-0.000562559, 0, 0]],
            rtol=1e-6
        )

    def test_get_z2spnorb(self):
        ls2 = get_ls2(self.sirius_rst, self.symorb, self.ao2soint,
            reference_occupation=(
                (2, 2, 2), (2, 2), (2, 2), (), (2, 2), (), (), ()
                )
            )
        numpy.testing.assert_allclose(
            ls2[2], [[0, -0.000562559, 0], [0.000562559, 0, 0], [0,0,0]],
            rtol=1e-6
        )
            
    @mock.patch('ledges.two_p_hole.get_eigen')
    @mock.patch('ledges.two_p_hole.makeV')
    @mock.patch('ledges.two_p_hole.get_ls1')
    @mock.patch('ledges.two_p_hole.SiriusRestart')
    @mock.patch('ledges.two_p_hole.tarfile.open')
    def test_untar(self, mock_open, mock_sir, mock_get_ls1, mock_makeV, mock_e):
        two_p_eigenvalues(self.dal_tar_gz, self.symorb)
        mock_open.return_value = mock.Mock()
        mock_open.assert_called_with(self.dal_tar_gz, 'r:gz')

    @mock.patch('ledges.two_p_hole.get_eigen')
    @mock.patch('ledges.two_p_hole.makeV')
    @mock.patch('ledges.two_p_hole.get_ls1')
    @mock.patch('ledges.two_p_hole.SiriusRestart')
    @mock.patch('ledges.two_p_hole.tarfile.open')
    def test_extract(self, mock_open, mock_sir, mock_get_ls1, mock_makeV, mock_e):
        mock_return_object = mock.Mock()
        mock_open.return_value = mock_return_object
        two_p_eigenvalues(self.dal_tar_gz, self.symorb)
        mock_return_object.extractall.assert_called_with(
            path='/tmp',
        )


if __name__ == "__main__":

    unittest.main()
