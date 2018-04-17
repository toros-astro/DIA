__version__ = '1.0a0'

import oismodule
import numpy as np

def oisdifference():
    zero = 2 * np.ones((10, 10))
    one = 3 * np.ones((10, 10))
    result = 5 * np.ones((10, 10))
    suma, greeting = oismodule.subtract(zero, one, 1, 1)
    is_it_zero = np.sum(suma - result) == 0.0
    print(greeting)
    print("Is it zero: {}".format(is_it_zero))
