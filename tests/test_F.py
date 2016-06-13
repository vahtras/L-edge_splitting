import unittest
import ledges
from ledges.constants import CM



class NewTest(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def full_splitting_test(self):
        """Niclass et al"""
        p_eigenvalues = ledges.two_p_hole.two_p_eigenvalues(
            "tests/F/lr_F.tar.gz",
            {2: (1,), 3: (1,), 5: (1,)},
            all_electron=True,
            reference_occupation=((2, 2), (2,), (2,), (), (2,), (), (), ()),
            )
        splitting = ledges.peaksplit.meansplit(2, 4, p_eigenvalues)
        self.assertAlmostEqual(splitting*CM, 376.53, places=2)


if __name__ == "__main__":
    unittest.main()
