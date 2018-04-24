'''
Creates CFFI library to implement a data transformation algorithm.


calculate_j(k, dim, shp, sub)

    Calculate the value of j and returns it.

    Parameters
    ----------
        k: size_t
        dim: size_t[]
        shp: size_t[]
        sub: size_t[]

    Returns
    -------
        j: size_t


calculate_i(j, shp, sho, sub)

    Calculate the value to increment and returns it.

    Parameters
    ----------
        k: size_t
        shp: size_t[]
        sho: size_t[]
        sub: size_t[]

    Returns
    -------
        i: size_t


calculate_cond(k, c, sho, sub)
    Calculate the condition and returns it.

    Parameters
    ----------
        k: size_t
        c: size_t
        sho: size_t[]
        sub: size_t[]

    Returns
    -------
        cond: int

tData(src, dest, typesize, sub_shape, shape, dimension, inverse)
    Choose if the algorithm needed is general or simple.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: size_t
            Data element size.
        sub_shape: size_t[]
            Data partition shape.
        shape: size_t[]
            Data shape.
        dimension: size_t
            Data dimension.
        inverse: size_t


tData_general(src, dest, typesize, sub_shape, shape, dimension, inverse)
    Calculate the data transformation in case that the algorithm needed is general.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: size_t
            Data element size.
        sub_shape: size_t[]
            Data partition shape.
        shape: size_t[]
            Data shape.
        dimension: size_t
            Data dimension.
        inverse: size_t


tData_general(src, dest, typesize, sub_shape, shape, dimension, inverse)
    Calculate the data transformation in case that the algorithm needed is simple.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: size_t
            Data element size.
        sub_shape: size_t[]
            Data partition shape.
        shape: size_t[]
            Data shape.
        dimension: size_t
            Data dimension.
        inverse: size_t
'''


from cffi import FFI
ffibuilder = FFI()

ffibuilder.set_source("tData",
'''
#include <stdio.h>
#include <string.h>

size_t calculate_j(size_t k, size_t dim[], size_t shp[], size_t sub[]) {

    size_t j = dim[0]*((k)%sub[0] + k/(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[0]/sub[0])*sub[0])
        +
        dim[1]*(k/(sub[0])%sub[1]*shp[0] + k/(shp[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[1]/sub[1])*shp[0]*sub[1])
        +
        dim[2]*(k/(sub[0]*sub[1])%sub[2]*shp[1]*shp[0] + k/(shp[0]*shp[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[2]/sub[2])*shp[0]*shp[1]*sub[2])
        +
        dim[3]*(k/(sub[0]*sub[1]*sub[2])%sub[3]*shp[0]*shp[1]*shp[2] + k/(shp[0]*shp[1]*shp[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[3]/sub[3])*shp[0]*shp[1]*shp[2]*sub[3])
        +
        dim[4]*(k/(sub[0]*sub[1]*sub[2]*sub[3])%sub[4]*shp[0]*shp[1]*shp[2]*shp[3] + k/(shp[0]*shp[1]*shp[2]*shp[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[4]/sub[4])*shp[0]*shp[1]*shp[2]*shp[3]*sub[4])
        +
        dim[5]*(k/(sub[0]*sub[1]*sub[2]*sub[3]*sub[4])%sub[5]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4] + k/(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5]*sub[6]*sub[7])%(shp[5]/sub[5])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*sub[5])
        +
        dim[6]*(k/(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5])%sub[6]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5] + k/(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6]*sub[7])%(shp[6]/sub[6])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*sub[6])
        +
        dim[7]*(k/(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6])%sub[7]*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6] + k/(shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7])%(shp[7]/sub[7])*shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*sub[7]);

        return j;
}

void tData_simple(char* src, char* dest, size_t typesize, size_t sub_shape[], size_t shape[], size_t dimension, size_t inverse) {

    size_t MAX_DIM = 8;
    size_t DIM = dimension;

    size_t dim[MAX_DIM], shp[MAX_DIM], sub[MAX_DIM];

    for (size_t i = 0; i < MAX_DIM; i++) {
        if (i < DIM) {
            dim[i] = 1;
            shp[DIM - i - 1] = shape[i];
            sub[DIM - i - 1] = sub_shape[i];
        } else {
            dim[i] = 0;
            shp[i] = 1;
            sub[i] = 1;
        }
    }

    size_t j;

    if (inverse == 0) {
        for (size_t k = 0; k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]; k++) {

            j = calculate_j(k, dim, shp, sub);


            memcpy(&dest[k * typesize], &src[j * typesize], typesize);
        }
    } else {
        for (size_t k = 0; k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]; k++) {

            j = calculate_j(k, dim, shp, sub);

            memcpy(&dest[j * typesize], &src[k * typesize], typesize);
        }
    }
}

int tData(char* src, char* dest, size_t typesize, size_t sub_shape[], size_t shape[], size_t dimension, size_t inverse) {

    tData_simple(src, dest, typesize, sub_shape, shape, dimension, inverse);
    return 1;
}

int k, k2;
void padData(char* src, char* dest, int typesize, int shape[], int pad_shape[], int dimension) {

    int MAX_DIM = 8;
    int DIM = dimension;

    int s[MAX_DIM], ps[MAX_DIM];

    for (int i = 0; i < MAX_DIM; i++) {
        if (i < DIM) {
            s[MAX_DIM + i - DIM] = shape[i];
            ps[MAX_DIM + i - DIM] = pad_shape[i];
        } else {
            s[MAX_DIM - i - 1] = 1;
            ps[MAX_DIM - i - 1] = 1;
        }
    }


    for (int a = 0; a < s[0]; a++) {
        for (int b = 0; b < s[1]; b++) {
            for (int c = 0; c < s[2]; c++) {
                for (int d = 0; d < s[3]; d++) {
                    for (int e = 0; e < s[4]; e++) {
                        for (int f = 0; f < s[5]; f++) {
                            for (int g = 0; g < s[6]; g++) {

                                k =  g*s[7]
                                     + f*s[7]*s[6]
                                     + e*s[7]*s[6]*s[5]
                                     + d*s[7]*s[6]*s[5]*s[4]
                                     + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                     + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                     + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1];

                                k2 = g*ps[7]
                                     + f*ps[7]*ps[6]
                                     + e*ps[7]*ps[6]*ps[5]
                                     + d*ps[7]*ps[6]*ps[5]*ps[4]
                                     + c*ps[7]*ps[6]*ps[5]*ps[4]*ps[3]
                                     + b*ps[7]*ps[6]*ps[5]*ps[4]*ps[3]*ps[2]
                                     + a*ps[7]*ps[6]*ps[5]*ps[4]*ps[3]*ps[2]*ps[1];

                                memcpy(&dest[k2 * typesize], &src[k *  typesize], s[7] * typesize);
                            }
                        }
                    }
                }
            }
        }
    }
}
''',
libraries=[])

ffibuilder.cdef(
                '''
                void tData(char* src, char* dest, size_t typesize, size_t sub_shape[],
                           size_t shape[], size_t dimension, size_t inverse);

                void tData_simple(char* src, char* dest, size_t typesize, size_t sub_shape[],
                                  size_t shape[], size_t dimension, size_t inverse);

                size_t calculate_j(size_t k, size_t dim[], size_t shp[], size_t sub[]);

                void padData(char* src, char* dest, int typesize, int shape[], int pad_shape[],
                             int dimension);

                '''
                )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
