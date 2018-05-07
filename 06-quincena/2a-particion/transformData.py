"""
Implements functions to use tData library.
"""

import numpy as np
from tData import ffi, lib
import pycblosc2 as cb2


def tData(src, ps, inverse=False):
    """
    Apply a data transformation based in reorganize data partitions.

    Parameters
    ----------
    src : np.array
        Data to transform.
    ps: int[] or tuple
        Data partition shape.
    inverse: bool, optional

    Returns
    -------
    dest: np.array
        Data transformed.
    d: dict
        Dictionary with partitions indexation.
    """

    # Obtain src parameters
    print("empieza aquí")
    typesize = src.dtype.itemsize
    shape = src.shape
    dimension = len(shape)

    # Calculate the extended dataset parameters

    ts = []

    for i in range(len(shape)):
        d = shape[i]
        if shape[i] % ps[i] != 0:
            d = shape[i] + ps[i] - shape[i] % ps[i]
        ts.append(d)

    size = np.prod(ts)

    # Create destination dataset

    dest = np.empty(size, dtype=src.dtype).reshape(ts)

    # Transform datasets to buffers (for use in cffi)

    src_b = ffi.from_buffer(src)
    dest_b = ffi.from_buffer(dest)

    # Execute the transformation

    inv = 1 if inverse else 0

    lib.tData(src_b, dest_b, typesize, shape, ts, ps, size, dimension, inv)

    d = createIndexation(ts, ps)

    return dest, d


def createIndexation(shape, ps):
    """
    Create an indexation of data partitions.

    Parameters
    ----------
    shape : int[] or tuple
        Data shape.
    ps: int[] or tuple
        Data partition shape.

    Returns
    -------
    indexation: dict
        Dictianary containing the indexation.
        key: partition number
        value: schunk number to decompress
    """

    dimension = len(shape)

    # Calculate the partitions the size

    size = int(np.prod(shape)/np.prod(ps))

    keys = np.zeros(size, dtype=np.int64)

    k_b = ffi.from_buffer(keys)

    # Create the indexation

    lib.createIndexation(k_b, shape, ps, dimension)

    # Create a dicttionary with keys and values

    d = dict([(k, v) for v, k in enumerate(keys)])

    return d


def obtainIndex(dim, dic, s, ps):

    dimension = len(s)

    s_aux = [1, 1, 1, 1, 1, 1, 1, 1]
    ps_aux = [1, 1, 1, 1, 1, 1, 1, 1]

    for i in range(dimension):
        s_aux[-dimension + i] = s[i]
        ps_aux[-dimension + i] = ps[i]

    ps = ps_aux
    s = s_aux

    ind = []

    for a in range(dim[0][0]//ps[0]*ps[0], dim[0][1], ps[0]):
        for b in range(dim[1][0]//ps[1]*ps[1], dim[1][1], ps[1]):
            for c in range(dim[2][0]//ps[2]*ps[2], dim[2][1], ps[2]):
                for d in range(dim[3][0]//ps[3]*ps[3], dim[3][1], ps[3]):
                    for e in range(dim[4][0]//ps[4]*ps[4], dim[4][1], ps[4]):
                        for f in range(dim[5][0]//ps[5]*ps[5], dim[5][1], ps[5]):
                            for g in range(dim[6][0]//ps[6]*ps[6], dim[6][1], ps[6]):
                                for h in range(dim[7][0]//ps[7]*ps[7], dim[7][1], ps[7]):

                                    k = (h
                                         + g*s[7]
                                         + f*s[7]*s[6]
                                         + e*s[7]*s[6]*s[5]
                                         + d*s[7]*s[6]*s[5]*s[4]
                                         + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                         + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                         + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1])

                                    ind.append((k, dic[k]))
    return ind


# Compression/decompression functions


def compress(src, ps):
    """
    Compress data.

    Parameters
    ----------
    src : np.array
        Data to compress.
    ps: int[] or tuple
        Data partition shape.

    Returns
    -------
    dest : chunk
        Data compressed.
    """

    size = src.size
    itemsize = src.dtype.itemsize
    bsize = size * itemsize

    dest = np.empty(size, dtype=src.dtype)

    cb2.blosc_set_blocksize(np.prod(ps) * itemsize)

    size_c = cb2.blosc_compress(5, 1, itemsize, bsize, src, dest, bsize)

    return dest


def decompress(comp, s, itemsize, dtype):
    """
    Decompress data.

    Parameters
    ----------
    comp : chunk
        Data compressed.
    s: int[] or tuple
        Original data shape.
    itemsize: int
        Data item size.
    dtype: np.type
        Data type.

    Returns
    -------
    dest : np.array
        Data decompressed.
    """

    size = np.prod(s)
    bsize = size * itemsize

    dest = np.empty(size, dtype=dtype).reshape(s)

    size_d = cb2.blosc_decompress(comp, dest, bsize)

    return dest


def get_block(comp, index, ps, dtype):
    """
    Extract a block of compressed data.

    Parameters
    ----------
    comp : chunk
        Data compressed.
    index: int
        Number of block to decompress.
    ps: int[] or tuple
        Data partition shape.
    dtype: np.type
        Data type.

    Returns
    -------
    dest : np.array
        Data compressed.
    """

    size = np.prod(ps)

    dest = np.empty(size, dtype=dtype)

    size_d = cb2.blosc_getitem(comp, index * size, size, dest)

    return dest


def decompress_trans(comp, indexation, dtype, s, ts, ps, a=-1, b=-1, c=-1, d=-1,
                     e=-1, f=-1, g=-1, h=-1):

    dimension = len(s)
    dim = [a, b, c, d, e, f, g, h][:-dimension]

    s_aux = [1, 1, 1, 1, 1, 1, 1, 1]
    ts_aux = [1, 1, 1, 1, 1, 1, 1, 1]
    ps_aux = [1, 1, 1, 1, 1, 1, 1, 1]
    dim_aux = [-1, -1, -1, -1, -1, -1, -1, -1]
    index_aux = [1, 1, 1, 1, 1, 1, 1, 1]

    for i in range(dimension):
        s_aux[-dimension + i] = s[i]
        ts_aux[-dimension + i] = ts[i]
        ps_aux[-dimension + i] = ps[i]
        dim_aux[-dimension - i] = dim[i]
    ts = ts_aux
    ps = ps_aux
    s = s_aux
    dim = dim_aux

    ul = np.copy(ts)
    ui = [(0, ts[0]), (0, ts[1]), (0, ts[2]), (0, ts[3]), (0, ts[4]), (0, ts[5]), (0, ts[6]),
          (0, ts[7])]

    for i in range(len(index_aux)):
        if dim[i] != -1:
            ul[i] = ps[i]
            ui[i] = (dim[i], dim[i]+1)
            index_aux[i] = 0

    subpl = ul

    ind = obtainIndex(ui, indexation, ts, ps)

    print(np.prod(subpl))

    dest = np.zeros(np.prod(subpl), dtype=dtype)

    for i, (k, n) in enumerate(ind):

        aux = get_block(comp, n, ps, dtype)

        h = k % subpl[7]
        g = k // (subpl[7]) % subpl[6]
        f = k // (subpl[7]*subpl[6]) % subpl[5]
        e = k // (subpl[7]*subpl[6]*subpl[5]) % subpl[4]
        d = k // (subpl[7]*subpl[6]*subpl[5]*subpl[4]) % subpl[3]
        c = k // (subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]) % subpl[2]
        b = k // (subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]) % subpl[1]
        a = k // (subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]*subpl[1]) % subpl[0]

        cont = 0

        for ra in range(a, a + ps[0]):
            for rb in range(b, b + ps[1]):
                for rc in range(c, c + ps[2]):
                    for rd in range(d, d + ps[3]):
                        for re in range(e, e + ps[4]):
                            for rf in range(f, f + ps[5]):
                                for rg in range(g, g + ps[6]):

                                    n = (h
                                         + rg*subpl[7]
                                         + rf*subpl[7]*subpl[6]
                                         + re*subpl[7]*subpl[6]*subpl[5]
                                         + rd*subpl[7]*subpl[6]*subpl[5]*subpl[4]
                                         + rc*subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]
                                         + rb*subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]
                                         + ra*subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]*subpl[1]) // ps[7]

                                    dest[n % np.prod(subpl) * ps[7]: (n % np.prod(subpl)+1) * ps[7]] = aux[cont * ps[7]: (cont+1) * ps[7]]
                                    cont += 1

    return dest.reshape(subpl[-dimension:])
