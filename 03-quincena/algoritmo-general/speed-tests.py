import transformData as td
import numpy as np
import time as t

shapes = ([6, 6, 16, 8, 8, 6, 8, 10], [2, 2, 4, 2, 4, 3, 2, 5])

shape, subshape = shapes

datasize = np.prod(np.array(shape))
datatype = np.int32

src = np.arange(datasize, dtype=datatype).reshape(shape)

bsize = datasize * src.dtype.itemsize

print('Dataset size is {:d} bytes'.format(bsize))

start = t.perf_counter()
dest = td.t_dataC(src, subshape, inverse=False)
end = t.perf_counter()
time = end - start
print('General algorithm time: {:.4f}s'.format(time))

start = t.perf_counter()
destS = td.t_dataC_simple(src, subshape, inverse=False)
end = t.perf_counter()
timeS = end - start
print('Simple algorithm time: {:.4f}s'.format(timeS))

ratio = time / timeS
print('The speed up is {:.2f}x'.format(ratio))

np.testing.assert_array_equal(dest, destS)
