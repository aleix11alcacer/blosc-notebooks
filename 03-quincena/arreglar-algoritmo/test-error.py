import transformData as td
import numpy as np

DIM = 6
shapesN = ([4, 4, 5, 3, 6, 6], [2, 3, 2, 2, 3, 3])

shapes = shapesN

shape, subshape = shapes
shapeN, subshapeN = shapesN

size = np.prod(np.array(shape))

srcN = np.arange(size, dtype=np.int32).reshape(shapeN)
src = np.arange(size, dtype=np.int32).reshape(shape)

destN = td.t_dataN(srcN, subshapeN).reshape(shape[:DIM])
destP = td.t_dataP(src, subshape).reshape(shape[:DIM])


for i in range(destN.size):
    a = destN.flatten()[i]
    b = destP.flatten()[i]
    if a != b:
        print(i, a, b)
np.testing.assert_array_equal(destN, destP)
