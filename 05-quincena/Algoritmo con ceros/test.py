import transformData3 as td
import pycblosc2 as cb2
import numpy as np
import pytest


@pytest.fixture(scope="module",
                params=[
                        ([16, 16, 16], [4, 4, 4]),
                        ([128, 128, 128], [8, 8, 8]),
                        ([64, 64, 64, 64], [8, 8, 8, 8]),
                        ([123, 124, 158], [32, 34, 31]),
                        ([42, 50, 74, 22], [32, 34, 31, 11])
                       ])
def shapes(request):
    SHAPE, PART_SHAPE = request.param
    yield (SHAPE, PART_SHAPE)


# --------------------------------------------------------------------------

def test_algoritm(shapes):
    SHAPE, PART_SHAPE = shapes
    size = np.prod(np.array(SHAPE))

    src = np.arange(size, dtype=np.int32).reshape(SHAPE)
    ITEMSIZE = src.dtype.itemsize

    src_part, indexation = td.tData(src, PART_SHAPE, inverse=False)
    TSHAPE = src_part.shape

    KB = 1024

    BLOSC_BLOCKSIZE = 16 * KB
    BLOSC_TYPESIZE = 4
    BLOSC_CODE = 5

    cparams = cb2.blosc2_create_cparams(compcode=BLOSC_CODE, clevel=8, use_dict=0,
                                        typesize=BLOSC_TYPESIZE, nthreads=4,
                                        blocksize=BLOSC_BLOCKSIZE, schunk=None,
                                        filters=[0, 0, 0, 0, 1], filters_meta=[0, 0, 0, 0, 0])

    dparams = cb2.blosc2_create_dparams(nthreads=4, schunk=None)

    schunk = td.compress(cparams, dparams, src)
    res = td.decompress(schunk, ITEMSIZE, SHAPE, b=15)

    schunk_t = td.compress_trans(cparams, dparams, src_part, TSHAPE, PART_SHAPE)
    res2 = td.decompress_trans(schunk_t, indexation, ITEMSIZE, SHAPE, TSHAPE, PART_SHAPE,
                               b=15)

    np.testing.assert_array_equal(res, res2)
