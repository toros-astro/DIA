import unittest
import oisdiff
import numpy as np
import scipy


class TestSubtraction(unittest.TestCase):

    def setUp(self):
        h, w = img_shape = (100, 100)
        n_stars = 10
        self.pos_x = np.random.randint(10, w - 10, n_stars)
        self.pos_y = np.random.randint(10, h - 10, n_stars)
        fluxes = 200.0 + np.random.rand(n_stars) * 300.0
        img = np.random.rand(h, w)
        ref = np.random.rand(h, w)
        for x, y, f in zip(self.pos_x, self.pos_y, fluxes):
            img[y, x] += f
            ref[y, x] += f

        from scipy.ndimage.filters import gaussian_filter
        self.img = gaussian_filter(img, sigma=1.7, mode='constant')
        self.ref = gaussian_filter(ref, sigma=0.8, mode='constant')

    def test_oisdifference(self):
        k_side = 7
        poly_deg = 0
        stamp_side = 23
        subt = oisdiff.oisdifference(
            self.img, self.ref, self.pos_x, self.pos_y, k_side, poly_deg,
            stamp_side,
            )
        subt_norm = np.linalg.norm(subt)
        ref_norm = np.linalg.norm(self.ref)
        self.assertLess(subt_norm / ref_norm, 0.1)


if __name__ == "__main__":
    unittest.main()
