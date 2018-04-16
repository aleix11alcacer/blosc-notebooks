'''
Implements functions to use tData library.
'''

import numpy as np
from tData import ffi, lib


def tData(src, sub_shp, inverse=False):
    """
    Execute tData function of tData library.

    Parameters
    ----------
    src : np.array
        Data to transform.
    sub_shp: int[] or tuple
        Data partition shape.
    inverse: bool, optional

    Returns
    -------
    dest: np.array
        Data transformed.
    """

    typesize = src.dtype.itemsize
    shape = src.shape
    dimension = len(shape)

    padShape = []

    for i in range(len(shape)):
        d = 0
        if shape[i] % sub_shp[i] != 0:
            d = sub_shp[i] - shape[i] % sub_shp[i]
        padShape.append((0, d))

    src = np.pad(src, padShape, 'constant')

    dest = np.zeros(src.size, dtype=src.dtype)

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    shape = src.shape

    lib.tData(src2, dest2, typesize, sub_shp, shape, dimension, inverse)

    return dest.reshape(shape)


def createIndexation(s, sb):
    dic = {}
    cont = 0
    for i in range(0, s[0], sb[0]):
        for j in range(0, s[1], sb[1]):
            for k in range(0, s[2], sb[2]):

                K = k + j*s[2] + i*s[1]*s[2]
                dic[K] = cont
                cont += 1
    return dic


def obtainIndex(x, y, z, dic, s, sb):

    ind = []

    for i in range(x[0]//32*32, x[1], sb[0]):
        for j in range(y[0]//32*32, y[1], sb[1]):
            for k in range(z[0]//32*32, z[1], sb[2]):

                K = k + j*s[2] + i*s[1]*s[2]

                ind.append(((i, j, k), dic[K]))
    return ind
