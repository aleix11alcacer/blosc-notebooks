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

    return dest


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


def decompress_trans(comp, ts, ps, a=-1, b=-1, c=-1, d=-1, e=-1, f=-1, g=-1, h=-1):
    """
    Decompress partitioned data.

    Parameters
    ----------
    comp : chunk
        Data compressed.
    ts: int[] or tuple
        Partitioned data shape.
    ps: int[] or tuple
        Data partition shape.
    a, b, c, d, e, f, g, h: int, optional
        Defines the subset of data desired

    Returns
    -------
    dest : np.array
     Data decompressed.
    """

    dimension = len(ps)
    dim = [a, b, c, d, e, f, g, h][:dimension]
    b_size = np.prod(ps)

    # Calculate desired data shape

    subpl = [1] * dimension

    for i in range(dimension):
        if i < dimension:
            if dim[i] != -1:
                subpl[i] = ps[i]
            else:
                subpl[i] = ts[i]

    dest = np.empty(np.prod(subpl), dtype=np.int32).reshape(subpl)

    dest_b = ffi.from_buffer(dest)
    comp_b = ffi.from_buffer(comp)

    lib.decompress_trans(comp_b, dest_b, ts, ps,  dim, subpl, dimension, b_size, dest.dtype.itemsize)

    return dest
