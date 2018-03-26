import transformData as td
import numpy as np
import pytest


@pytest.fixture(scope="module",
                params=[
                        # ([6, 4, 6, 8], [2, 2, 2, 2]),
                        # ([6, 4, 6, 8], [3, 2, 2, 4]),
                        # ([5, 17, 13, 8], [2, 9, 3, 3]),  # 4 dimensiones
                        #
                        # ([12, 18, 6, 9, 6], [3, 3, 3, 3, 3]),
                        # ([8, 12, 6, 9, 6], [3, 4, 3, 3, 2]),
                        # ([5, 17, 5, 3, 7], [2, 13, 2, 2, 3]),  # 5 dimensiones
                        #
                        # ([6, 8, 4, 6, 8, 6], [2, 2, 2, 2, 2, 2]),
                        # ([6, 8, 4, 6, 8, 6], [3, 2, 2, 3, 4, 2]),
                        # ([5, 4, 4, 6, 8], [2, 3, 2, 2, 3]),
                        # ([4, 4, 5, 3, 6, 6], [2, 3, 2, 2, 3, 3]),  # 6 dimensiones
                        #
                        # ([6, 8, 6, 8, 6, 8, 6, 4], [2, 2, 2, 2, 2, 2, 2, 2]),
                        ([10, 13, 11, 17, 5, 8, 2, 3], [6, 3, 4, 6, 3, 2, 2, 2])  # 8 dimensiones
                       ])
def shapes(request):
    shape, subshape = request.param
    yield (shape, subshape)


# --------------------------------------------------------------------------

def test_algoritm(shapes):
    shape, subshape = shapes
    size = np.prod(np.array(shape))
    src = np.arange(size, dtype=np.int32).reshape(shape)

    destN = td.t_dataN(src, subshape)
    destP = td.t_dataP(src, subshape)
    destC = td.t_dataP(src, subshape)

    np.testing.assert_array_equal(destC, destP)
    np.testing.assert_array_equal(destN, destC)
    np.testing.assert_array_equal(destN, destP)


def test_inverse(shapes):
    shape, subshape = shapes
    size = np.prod(np.array(shape))
    src = np.arange(size, dtype=np.int32).reshape(shape)
    dest = np.zeros(size, dtype=src.dtype).reshape(src.shape)
    src2 = np.zeros(size, dtype=src.dtype).reshape(src.shape)

    dest = td.t_dataC(src, subshape, inverse=False)

    src2 = td.t_dataC(dest, subshape, inverse=True)

    np.testing.assert_array_equal(src, src2)
