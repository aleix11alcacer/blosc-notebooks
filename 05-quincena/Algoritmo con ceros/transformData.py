'''
Implements functions to use tData library.
'''

import numpy as np
from tData import ffi, lib
import pycblosc2 as cb2


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

    return dest, shape


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


# Compression/decompression functions

def compress(cparams, dparams, src):

    bsize = src.size * src.dtype.itemsize

    schunk = cb2.blosc2_new_schunk(cparams, dparams)

    nchunks = cb2.blosc2_append_buffer(schunk, bsize, src)

    return schunk


def decompress(schunk, item_size, s, x=slice(0, None), y=slice(0, None), z=slice(0, None)):

    size = np.prod(s)
    bsize = size * item_size

    dest = np.zeros(size, dtype=np.int32).reshape(s)

    cb2.blosc2_decompress_chunk(schunk, 0, dest, bsize)
    cb2.blosc2_free_schunk(schunk)

    return dest[x, y, z]


def compress_trans(cparams, dparams, srct, ts, sb):

    a_size = np.prod(sb)
    a_bsize = a_size * srct.dtype.itemsize

    schunk = cb2.blosc2_new_schunk(cparams, dparams)

    for i in range(ts[0]//sb[0] * ts[1]//sb[1] * ts[2]//sb[2]):
        aux = srct[i * a_size:(i+1) * a_size]
        nchunks = cb2.blosc2_append_buffer(schunk, a_bsize, aux)

    return schunk


def decompress_trans(schunk, indexation, item_size, s, ts, sb, x=slice(0, None), y=slice(0, None),
                     z=slice(0, None)):
    index_aux = [1, 1, 1]
    xl, yl, zl = ts
    xi, yi, zi = (0, ts[0]), (0, ts[1]), (0, ts[2])

    if x != slice(0, None):
        xl = sb[0]
        xi = (x, x+1)
        index_aux[0] = 0
        x = x % xl
    if y != slice(0, None):
        yl = sb[1]
        yi = (y, y+1)
        index_aux[1] = 0
        y = y % yl
    if z != slice(0, None):
        zl = sb[2]
        zi = (z, z+1)
        index_aux[2] = 0
        z = z % zl
    SUBPL = [xl, yl, zl]

    ind = obtainIndex(xi, yi, zi, indexation, ts, sb)

    dest = np.zeros(np.prod(SUBPL), dtype=np.int32).reshape(SUBPL)

    AUX_SIZE = np.prod(sb)
    AUX_bsize = AUX_SIZE * item_size

    aux = np.zeros(AUX_SIZE, dtype=np.int32).reshape(sb)

    for index, n in ind:
        i, j, k = [index[q]*index_aux[q] for q in range(3)]
        cb2.blosc2_decompress_chunk(schunk, n, aux, AUX_bsize)
        dest[i:i+sb[0], j:j+sb[1], k:k+sb[2]] = aux

    return dest[x, y, z]
