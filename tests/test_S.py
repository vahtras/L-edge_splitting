import unittest
import mock
import sys
import os
import numpy
import ledges.two_p_hole as l2p

from daltools.sirrst import SiriusRestart
from daltools import prop

labels = ('X1SPNORB', 'Y1SPNORB', 'Z1SPNORB')

class S_Test(unittest.TestCase):

    def setUp(self):
        self.subdir = 'tests/S'
        self.mol = 'S'
        self.dal_tar_gz = "%s/hf_%s.tar.gz" % (self.subdir, self.mol)
        restart_file = "%s/hf_%s.SIRIUS.RST" % (self.subdir, self.mol)
        self.sirius_rst = SiriusRestart(name=restart_file)
        self.cmo = self.sirius_rst.cmo
        self.aoproper = "%s/hf_%s.AOPROPER" % (self.subdir, self.mol)
        self.ao2soint = "%s/hf_%s.AO2SOINT" % (self.subdir, self.mol)
        self.ls = prop.read(*labels, filename=self.aoproper)
        self.symorb = {2:(1,), 3:(1,), 5:(1,)}

    def tearDown(self):
        pass

    def test_setup(self):
        self.assertTupleEqual(tuple(self.cmo.nrow), (9, 6, 6, 2, 6, 2, 2, 1))

    def test_get_2p_orbitals(self):
        p_orbitals = l2p.get_orbitals(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        ref = numpy.zeros((34, 3))
        ref[9:15,  0] =  [1.00179538, 0.00115419, -0.00341405, -0.00112856, 0.00032629, 0.]
        ref[15:21, 1] =  [1.00179538, 0.00115419, -0.00341405, -0.00112856, 0., 0.00032629]
        ref[23:29, 2] =  [0.99656939, -0.00477253, 0.00649723, 0.00218673, 0.00041297, 0.]
        
        numpy.testing.assert_allclose(p_orbitals, ref, atol=1e-4)

    def test_get_indices(self):
        orbital_indices = l2p.get_orbital_indices(self.cmo, {2:(1,), 3:(1,), 5:(1,)})
        self.assertTupleEqual(orbital_indices, (9, 15, 23))

    def test_get_x1spnorb(self):
        ls = l2p.get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[0], [[0,0,0], [0, 0, 0.03385742 ], [0, -0.03385742, 0]]
        )

    def test_get_y1spnorb(self):
        ls = l2p.get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[1], [[0, 0, -0.03385742 ],[0,0,0],  [0.03385742, 0, 0]]
        )

    def test_get_z1spnorb(self):
        ls = l2p.get_ls1(self.cmo, self.symorb, self.aoproper)
        numpy.testing.assert_allclose(
            ls[2], [[0, 0.03396038,  0], [-0.03396038,   0, 0], [0,0,0]],
            rtol=1e-6
        )

    def test_get_x2spnorb(self):
        ls2 = l2p.get_ls2(self.sirius_rst, self.symorb, self.ao2soint)
        numpy.testing.assert_allclose(
            ls2[0], [[0,0,0], [0, 0, -0.00619585 ], [0, 0.00619585, 0]],
            rtol=1e-6
        )

    def test_get_y2spnorb(self):
        ls2 = l2p.get_ls2(self.sirius_rst, self.symorb, self.ao2soint)
        numpy.testing.assert_allclose(
            ls2[1], [[0, 0, 0.00619585], [0,0,0], [-0.00619585, 0, 0]],
            rtol=1e-6
        )

    def test_get_z2spnorb(self):
        ls2 = l2p.get_ls2(self.sirius_rst, self.symorb, self.ao2soint)
        numpy.testing.assert_allclose(
            ls2[2], [[0, -0.00616044, 0], [0.00616044, 0, 0], [0,0,0]],
            rtol=1e-6
        )

    def test_make_V(self):
        units = [i*numpy.identity(3) for i in (2, 4, 6)]
        expect = [[3j, 0, 0, 1j - 2, 0, 0],
                  [0, 3j, 0, 0, 1j - 2, 0],
                  [0, 0, 3j, 0, 0, 1j - 2],
                  [1j + 2, 0, 0, -3j, 0, 0],
                  [0, 1j + 2, 0, 0, -3j, 0],
                  [0, 0, 1j + 2, 0, 0, -3j]]
        numpy.testing.assert_allclose(l2p.makeV(units), expect)   

    def test_eigenvalues(self):
        V = [[0, .5j], [-.5j, 0]]
        numpy.testing.assert_allclose(
            l2p.get_eigen(V),
            [-0.5, 0.5]
        )
            
    @mock.patch('ledges.two_p_hole.get_eigen')
    @mock.patch('ledges.two_p_hole.makeV')
    @mock.patch('ledges.two_p_hole.get_ls1')
    @mock.patch('ledges.two_p_hole.SiriusRestart')
    @mock.patch('ledges.two_p_hole.tarfile.open')
    def test_untar(self, mock_open, mock_sir, mock_get_ls1, mock_makeV, mock_e):
        l2p.two_p_eigenvalues(self.dal_tar_gz, self.symorb)
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
        l2p.two_p_eigenvalues(self.dal_tar_gz, self.symorb)
        mock_return_object.extractall.assert_called_with(
            path='/tmp',
        )


if __name__ == "__main__":

    unittest.main()
