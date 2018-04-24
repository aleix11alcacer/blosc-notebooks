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

    src = padData(src, sub_shp)

    dest = np.zeros(src.size, dtype=src.dtype)

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    shape = src.shape

    lib.tData(src2, dest2, typesize, sub_shp, shape, dimension, inverse)

    return dest, shape


def padData(src, sub_shp):
    """
    Pad matrix with 0. Similar to numpy.pad()

    Parameters
    ----------
    src : np.array
        Data to transform.
    sub_shp: int[] or tuple
        Data partition shape.

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
        d = shape[i]
        if shape[i] % sub_shp[i] != 0:
            d = shape[i] + sub_shp[i] - shape[i] % sub_shp[i]
        padShape.append(d)

    dest = np.zeros(np.prod(padShape), dtype=src.dtype).reshape(padShape)

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    lib.padData(src2, dest2, typesize, shape, padShape, dimension)

    return dest


def reorder_dim(dim, l):
    aux = []
    for i in range(len(dim)):
        aux.append(dim[(i+l) % len(dim)])
    return aux


def create_shape(s):
    s_aux = [1, 1, 1, 1, 1, 1, 1, 1]

    for i in range(len(s)):
        s_aux[-len(s) + i] = s[i]
    return s_aux


def createIndexation(s, sb):

    s = create_shape(s)
    sb = create_shape(sb)

    dic = {}
    cont = 0

    for a in range(0, s[0], sb[0]):
        for b in range(0, s[1], sb[1]):
            for c in range(0, s[2], sb[2]):
                for d in range(0, s[3], sb[3]):
                    for e in range(0, s[4], sb[4]):
                        for f in range(0, s[5], sb[5]):
                            for g in range(0, s[6], sb[6]):
                                for h in range(0, s[7], sb[7]):

                                    K = (h
                                         + g*s[7]
                                         + f*s[7]*s[6]
                                         + e*s[7]*s[6]*s[5]
                                         + d*s[7]*s[6]*s[5]*s[4]
                                         + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                         + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                         + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1])

                                    dic[K] = cont
                                    cont += 1
    return dic


def obtainIndex(dim, dic, s, sb):

    s = create_shape(s)
    sb = create_shape(sb)

    ind = []

    for a in range(dim[0][0]//sb[0]*sb[0], dim[0][1], sb[0]):
        for b in range(dim[1][0]//sb[1]*sb[1], dim[1][1], sb[1]):
            for c in range(dim[2][0]//sb[2]*sb[2], dim[2][1], sb[2]):
                for d in range(dim[3][0]//sb[3]*sb[3], dim[3][1], sb[3]):
                    for e in range(dim[4][0]//sb[4]*sb[4], dim[4][1], sb[4]):
                        for f in range(dim[5][0]//sb[5]*sb[5], dim[5][1], sb[5]):
                            for g in range(dim[6][0]//sb[6]*sb[6], dim[6][1], sb[6]):
                                for h in range(dim[7][0]//sb[7]*sb[7], dim[7][1], sb[7]):

                                    K = (h
                                         + g*s[7]
                                         + f*s[7]*s[6]
                                         + e*s[7]*s[6]*s[5]
                                         + d*s[7]*s[6]*s[5]*s[4]
                                         + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                         + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                         + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1])

                                    ind.append(((a, b, c, d, e, f, g, h), dic[K]))
    return ind


# Compression/decompression functions

def compress(cparams, dparams, src):

    bsize = src.size * src.dtype.itemsize

    schunk = cb2.blosc2_new_schunk(cparams, dparams)

    nchunks = cb2.blosc2_append_buffer(schunk, bsize, src)

    return schunk


def decompress(schunk, item_size, s, a=slice(0, None), b=slice(0, None), c=slice(0, None),
               d=slice(0, None), e=slice(0, None), f=slice(0, None), g=slice(0, None),
               h=slice(0, None)):

    dim = reorder_dim([a, b, c, d, e, f, g, h], len(s))

    s = create_shape(s)

    size = np.prod(s)
    bsize = size * item_size

    dest = np.zeros(size, dtype=np.int32).reshape(s)

    cb2.blosc2_decompress_chunk(schunk, 0, dest, bsize)
    cb2.blosc2_free_schunk(schunk)

    return dest[dim]


def compress_trans(cparams, dparams, srct, ts, sb):

    a_size = np.prod(sb)
    a_bsize = a_size * srct.dtype.itemsize

    schunk = cb2.blosc2_new_schunk(cparams, dparams)

    for i in range(np.prod([ts[j]//sb[j] for j in range(len(ts))])):
        aux = srct[i * a_size:(i+1) * a_size]
        nchunks = cb2.blosc2_append_buffer(schunk, a_bsize, aux)

    return schunk


def decompress_trans(schunk, indexation, item_size, s, ts, sb, a=None, b=None, c=None, d=None,
                     e=None, f=None, g=None, h=None):

    dim = reorder_dim([a, b, c, d, e, f, g, h], len(s))

    s = create_shape(s)
    ts = create_shape(ts)
    sb = create_shape(sb)

    index_aux = [1, 1, 1, 1, 1, 1, 1, 1]
    ul = np.copy(ts)
    ui = [(0, ts[0]), (0, ts[1]), (0, ts[2]), (0, ts[3]), (0, ts[4]), (0, ts[5]), (0, ts[6]),
          (0, ts[7])]

    for i in range(len(index_aux)):
        if dim[i] is not None:
            ul[i] = sb[i]
            ui[i] = (dim[i], dim[i]+1)
            index_aux[i] = 0
            dim[i] = dim[i] % ul[i]
        else:
            dim[i] = slice(0, s[i])

    subpl = ul

    ind = obtainIndex(ui, indexation, ts, sb)

    dest = np.zeros(np.prod(subpl), dtype=np.int32).reshape(subpl)

    aux_size = np.prod(sb)
    AUX_bsize = aux_size * item_size

    aux = np.zeros(aux_size, dtype=np.int32).reshape(sb)

    for index, n in ind:

        index = [index[q]*index_aux[q] for q in range(len(dim))]
        cb2.blosc2_decompress_chunk(schunk, n, aux, AUX_bsize)

        dest[index[0]:index[0]+sb[0],
             index[1]:index[1]+sb[1],
             index[2]:index[2]+sb[2],
             index[3]:index[3]+sb[3],
             index[4]:index[4]+sb[4],
             index[5]:index[5]+sb[5],
             index[6]:index[6]+sb[6],
             index[7]:index[7]+sb[7]] = aux

    return dest[dim]
