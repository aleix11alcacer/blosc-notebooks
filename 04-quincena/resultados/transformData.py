'''
Módulo que contiene los tres algoritmos (Numpy, Python, C) para ser usados desde Python.
'''

import numpy as np
from tdg import ffi as ffi
from tdg import lib as lib

from tdgSimple import ffi as ffiS
from tdgSimple import lib as libS


def t_dataP(data, sub_shp, inverse=False):

    dim = [0, 0, 0, 0, 0, 0, 0, 0]
    sho = [1, 1, 1, 1, 1, 1, 1, 1]
    shp = [1, 1, 1, 1, 1, 1, 1, 1]
    sub = [1, 1, 1, 1, 1, 1, 1, 1]

    dest = np.zeros(data.shape, dtype=data.dtype).flatten()

    org_shp = list(data.shape)

    dta_shp = [org_shp[i] - org_shp[i] % sub_shp[i] for i in range(len(org_shp))]

    for i in range(len(sub_shp)):
        dim[len(sub_shp) - 1 - i] = 1
        shp[len(sub_shp) - 1 - i] = dta_shp[i]
        sub[len(sub_shp) - 1 - i] = sub_shp[i]
        sho[len(sub_shp) - 1 - i] = org_shp[i]

    print(dim, shp, sub, sho)
    cont = 0

    K = 0

    data_aux = np.copy(data.flatten())

    if not inverse:

        while K < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]:

            if (
                (K+cont) % sho[0] >= sho[0] - sho[0] % sub[0]
                or (K+cont) % (sho[0]*sho[1])//(sho[0]) >= sho[1] - sho[1] % sub[1]
                or (K+cont) % (sho[0]*sho[1]*sho[2])//(sho[0]*sho[1]) >= sho[2] - sho[2] % sub[2]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3])//(sho[0]*sho[1]*sho[2]) >= sho[3] - sho[3] % sub[3]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4])//(sho[0]*sho[1]*sho[2]*sho[3]) >= sho[4] - sho[4] % sub[4]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]) >= sho[5] - sho[5] % sub[5]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]) >= sho[6] - sho[6] % sub[6]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]) >= sho[7] - sho[7] % sub[7]):

                dest[K + cont] = data_aux[K + cont]
                cont += 1

            else:

                L = K + cont

                J = (dim[0] * ((K) % sub[0] + K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[0]//sub[0])*sub[0])
                     +
                     dim[1] * (K//(sub[0]) % sub[1]*shp[0] + K//(shp[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[1]//sub[1])*shp[0]*sub[1])
                     +
                     dim[2] * (K//(sub[0]*sub[1]) % sub[2]*shp[1]*shp[0] + K//(shp[0]*shp[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[2]//sub[2])*shp[0]*shp[1]*sub[2])
                     +
                     dim[3] * (K//(sub[0]*sub[1]*sub[2]) % sub[3]*shp[0]*shp[1]*shp[2] + K//(shp[0]*shp[1]*shp[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[3]//sub[3])*shp[0]*shp[1]*shp[2]*sub[3])
                     +
                     dim[4] * (K//(sub[0]*sub[1]*sub[2]*sub[3]) % sub[4]*shp[0]*shp[1]*shp[2]*shp[3] + K//(shp[0]*shp[1]*shp[2]*shp[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[4]//sub[4])*shp[0]*shp[1]*shp[2]*shp[3]*sub[4])
                     +
                     dim[5] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]) % sub[5]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5]*sub[6]*sub[7]) % (shp[5]//sub[5])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5])
                     +
                     dim[6] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]) % sub[6]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6]*sub[7]) % (shp[6]//sub[6])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6])
                     +
                     dim[7] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]) % sub[7]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7]) % (shp[7]//sub[7])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7])
                     )

                inc = (dim[0] * (J)//shp[0] * (sho[0] % sub[0])
                       + dim[1] * (J)//(shp[1]*shp[0]) * (sho[1] % sub[1])*sho[0]
                       + dim[2] * (J)//(shp[2]*shp[1] * shp[0]) * (sho[2] % sub[2]) * sho[0] * sho[1]
                       + dim[3] * (J)//(shp[3]*shp[2]*shp[1]*shp[0]) * (sho[3] % sub[3]) * sho[0]*sho[1]*sho[2]
                       + dim[4] * (J)//(shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[4] % sub[4]) * sho[0]*sho[1]*sho[2]*sho[3]
                       + dim[5] * (J)//(shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[5] % sub[5]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]
                       + dim[6] * (J)//(shp[6]*shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[6] % sub[6]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]
                       + dim[7] * (J)//(shp[7]*shp[6]*shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[7] % sub[7]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6])

                dest[L] = data_aux[J + inc]

                K += 1

        for L in range(K + cont, sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7]):
            dest[L] = data_aux[L]

    else:

        while K < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]:

            if ((K+cont) % sho[0] >= sho[0] - sho[0] % sub[0]
                or (K+cont) % (sho[0]*sho[1])//(sho[0]) >= sho[1] - sho[1] % sub[1]
                or (K+cont) % (sho[0]*sho[1]*sho[2])//(sho[0]*sho[1]) >= sho[2] - sho[2] % sub[2]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3])//(sho[0]*sho[1]*sho[2]) >= sho[3] - sho[3] % sub[3]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4])//(sho[0]*sho[1]*sho[2]*sho[3]) >= sho[4] - sho[4] % sub[4]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]) >= sho[5] - sho[5] % sub[5]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]) >= sho[6] - sho[6] % sub[6]
                or (K+cont) % (sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7])//(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]) >= sho[7] - sho[7] % sub[7]):

                dest[K + cont] = data_aux[K + cont]
                cont += 1

            else:

                L = K + cont

                J = (dim[0] * ((K) % sub[0] + K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[0]//sub[0])*sub[0])
                     +
                     dim[1] * (K//(sub[0]) % sub[1]*shp[0] + K//(shp[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[1]//sub[1])*shp[0]*sub[1])
                     +
                     dim[2] * (K//(sub[0]*sub[1]) % sub[2]*shp[1]*shp[0] + K//(shp[0]*shp[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[2]//sub[2])*shp[0]*shp[1]*sub[2])
                     +
                     dim[3] * (K//(sub[0]*sub[1]*sub[2]) % sub[3]*shp[0]*shp[1]*shp[2] + K//(shp[0]*shp[1]*shp[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[3]//sub[3])*shp[0]*shp[1]*shp[2]*sub[3])
                     +
                     dim[4] * (K//(sub[0]*sub[1]*sub[2]*sub[3]) % sub[4]*shp[0]*shp[1]*shp[2]*shp[3] + K//(shp[0]*shp[1]*shp[2]*shp[3]*sub[4]*sub[5]*sub[6]*sub[7]) % (shp[4]//sub[4])*shp[0]*shp[1]*shp[2]*shp[3]*sub[4])
                     +
                     dim[5] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]) % sub[5]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5]*sub[6]*sub[7]) % (shp[5]//sub[5])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5])
                     +
                     dim[6] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]) % sub[6]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6]*sub[7]) % (shp[6]//sub[6])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6])
                     +
                     dim[7] * (K//(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]) % sub[7]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6] + K//(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7]) % (shp[7]//sub[7])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7])
                     )

                inc = (dim[0] * (J)//shp[0] * (sho[0] % sub[0])
                       + dim[1] * (J)//(shp[1]*shp[0]) * (sho[1] % sub[1])*sho[0]
                       + dim[2] * (J)//(shp[2]*shp[1] * shp[0]) * (sho[2] % sub[2]) * sho[0] * sho[1]
                       + dim[3] * (J)//(shp[3]*shp[2]*shp[1]*shp[0]) * (sho[3] % sub[3]) * sho[0]*sho[1]*sho[2]
                       + dim[4] * (J)//(shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[4] % sub[4]) * sho[0]*sho[1]*sho[2]*sho[3]
                       + dim[5] * (J)//(shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[5] % sub[5]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]
                       + dim[6] * (J)//(shp[6]*shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[6] % sub[6]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]
                       + dim[7] * (J)//(shp[7]*shp[6]*shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[7] % sub[7]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6])

                dest[J + inc] = data_aux[L]

                K += 1

        for L in range(K + cont, sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7]):
            dest[L] = data_aux[L]

    return dest.reshape(data.shape)


def t_dataC(src, sub_shp, inverse=False):

    dest = np.zeros(src.size, dtype=src.dtype).reshape(src.shape)

    typesize = src.dtype.itemsize
    shape = src.shape
    dimension = len(shape)

    src2 = ffi.from_buffer(src)
    dest2 = ffi.from_buffer(dest)

    lib.transform_data(src2, dest2, typesize, sub_shp, shape, dimension, inverse)

    return dest


def t_dataC_simple(src, sub_shp, inverse=False):

    dest = np.zeros(src.size, dtype=src.dtype).reshape(src.shape)

    typesize = src.dtype.itemsize
    shape = src.shape
    dimension = len(shape)

    src2 = ffiS.from_buffer(src)
    dest2 = ffiS.from_buffer(dest)

    libS.transform_data(src2, dest2, typesize, sub_shp, shape, dimension, inverse)

    return dest
