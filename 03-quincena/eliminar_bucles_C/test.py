from transform_data import ffi, lib
import numpy as np


def transform_data(src, dest, sub_shape, inverse=0):

    typesize = src.dtype.itemsize

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    shape = src.shape

    dimension = len(src.shape)

    lib.transform_data(src2, dest2, typesize, sub_shape, shape, dimension, inverse)


src = np.arange(128*512*64, dtype=np.int64).reshape(128, 512, 64)
dest = np.empty(128*512*64, dtype=np.int64).reshape(128, 512, 64)
src2 = np.empty(128*512*64, dtype=np.int64).reshape(128, 512, 64)

sub_shape = [8, 4, 16]


transform_data(src, dest, sub_shape)
transform_data(dest, src2, sub_shape, inverse=1)

np.testing.assert_array_equal(src, src2)
