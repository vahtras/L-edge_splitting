import unittest
import mock
import numpy

from  . import ledges
from ledges import peaksplit

class NewTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @mock.patch('numpy.loadtxt')
    def test_24(self, mock_loadtxt):
        mock_loadtxt.return_value = numpy.array(
            [0.9, 1.1, 1.9, 1.95, 2.05, 2.1]
            )
        ms = peaksplit.mean(2, 4, 'dummy')
        self.assertAlmostEqual(ms, 1.0)


if __name__ == "__main__":
        unittest.main()
