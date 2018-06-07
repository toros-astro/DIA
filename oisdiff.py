import oismodule

__version__ = '1.0a1'


def oisdifference(sci_img, ref_img, xc, yc, k_side=5, poly_deg=1, stamp_side=10):
    return oismodule.subtract(
        sci_img.astype('f8'), ref_img.astype('f8'),
        xc.astype('i4'), yc.astype('i4'),
        k_side, poly_deg, stamp_side)
