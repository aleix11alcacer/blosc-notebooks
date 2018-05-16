"""
Test para comprobar que el algoritmo implementado de decompresión funciona
"""

import transformData as td
import numpy as np
import pytest


@pytest.fixture(scope="module",
                params=[
                        ([16, 16, 16], [4, 4, 4], [5, -1, 12, -1, -1, -1, -1, -1]),
                        ([128, 128, 128], [8, 8, 8], [-1, 127, -1, -1, -1, -1, -1, -1]),
                        ([64, 64, 64, 64], [8, 8, 8, 8], [23, -1, -1, 56, -1, -1, -1, -1]),
                        ([123, 124, 158], [32, 34, 31], [13, -1, 119, -1, -1, -1, -1, -1]),
                        ([4, 4, 4, 4], [2, 2, 3, 3], [-1, 3, 3, -1, -1, -1, -1, -1]),
                        ([8, 8, 8, 8, 8, 8, 8, 8], [2, 5, 7, 2, 3, 4, 6, 3],
                         [-1, 4, 6, -1, 5, -1, -1, 0])
                       ])
def shapes(request):
    shape, part_shape, index = request.param
    yield (shape, part_shape, index)


# --------------------------------------------------------------------------

def test_algoritm(shapes):
    shape, part_shape, index = shapes

    # Create slices

    slices = []
    for i in range(len(shape)):
        if index[i] != -1:
            slices.append(slice(index[i], index[i] + 1))
        else:
            slices.append(slice(0, shape[i]))
    slices = tuple(slices)

    # Create indexes

    a = index[0]
    b = index[1]
    c = index[2]
    d = index[3]
    e = index[4]
    f = index[5]
    g = index[6]
    h = index[7]

    # Test

    size = np.prod(shape)

    src = np.arange(size, dtype=np.int32).reshape(shape)

    src_trans = td.tData(src, part_shape)

    trans_shape = src_trans.shape

    comp = td.compress(src_trans, part_shape)

    dest = td.decompress_trans(comp, shape, trans_shape, part_shape, a, b, c, d, e, f, g, h)

    np.testing.assert_array_equal(src[slices], dest)
