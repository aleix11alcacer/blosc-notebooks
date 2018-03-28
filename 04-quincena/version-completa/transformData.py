'''
Módulo que llama a las funciones de la libreria tData.
'''

import numpy as np
from tData import ffi, lib


def tData(src, sub_shp, inverse=False):

    dest = np.zeros(src.size, dtype=src.dtype).reshape(src.shape)

    typesize = src.dtype.itemsize
    shape = src.shape
    dimension = len(shape)

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    lib.tData(src2, dest2, typesize, sub_shp, shape, dimension, inverse)

    return dest
