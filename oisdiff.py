import oismodule

__version__ = '1.0a0'


def oisdifference(sci_img, ref_img, xc, yc, k_side, poly_deg, stamp_side):
    return oismodule.subtract(
        sci_img.astype('f8'), ref_img.astype('f8'),
        xc.astype('i4'), yc.astype('i4'),
        k_side, poly_deg, stamp_side)
