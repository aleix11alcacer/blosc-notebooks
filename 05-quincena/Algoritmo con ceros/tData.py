'''
Creates CFFI library to implement a data transformation algorithm.


calculate_j(k, dim, shp, sub)

    Calculate the value of j and returns it.

    Parameters
    ----------
        k: int
        dim: int[]
        shp: int[]
        sub: int[]

    Returns
    -------
        j: int


tData(src, dest, typesize, sub_shape, shape, dimension, inverse)
    Choose if the algorithm needed is general or simple.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: int
            Data element size.
        sub_shape: int[]
            Data partition shape.
        shape: int[]
            Data shape.
        dimension: int
            Data dimension.
        inverse: int


tData_simple(src, dest, typesize, sub_shape, shape, dimension, inverse)
    Calculate the data transformation in case that the algorithm needed is simple.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: int
            Data element size.
        sub_shape: int[]
            Data partition shape.
        shape: int[]
            Data shape.
        dimension: int
            Data dimension.
        inverse: int
'''


from cffi import FFI
ffibuilder = FFI()

ffibuilder.set_source("tData",
'''
#include <stdio.h>
#include <string.h>
#include <stdint.h>

int calculate_j(int k, int dim[], int shp[], int sub[]) {

    int j = dim[0]*((k)%sub[0] + k/(sub[0]*sub[1]*sub[2]*sub[3]*sub[4]*sub[5]*sub[6]*sub[7])%(shp[0]/sub[0])*sub[0])
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

void createIndexation(char* keys, int shape[], int sub_shape[], int dimension){

    int MAX_DIM = 8;
    int DIM = dimension;

    int s[MAX_DIM], sb[MAX_DIM];

    for (int i = 0; i < MAX_DIM; i++) {
        if (i < DIM) {
            s[MAX_DIM + i - DIM] = shape[i];
            sb[MAX_DIM + i - DIM] = sub_shape[i];
        } else {
            s[MAX_DIM - i - 1] = 1;
            sb[MAX_DIM - i - 1] = 1;
        }
    }

    int cont = 0;
    uint64_t k;

    for (int a=0; a < s[0]; a += sb[0]){
        for (int b=0; b < s[1]; b += sb[1]){
            for (int c=0; c < s[2]; c += sb[2]){
                for (int d=0; d < s[3]; d += sb[3]){
                    for (int e=0; e < s[4]; e += sb[4]){
                        for (int f=0; f < s[5]; f += sb[5]){
                            for (int g=0; g < s[6]; g += sb[6]){
                                for (int h=0; h < s[7]; h += sb[7]){

                                     k = h
                                          + g*s[7]
                                          + f*s[7]*s[6]
                                          + e*s[7]*s[6]*s[5]
                                          + d*s[7]*s[6]*s[5]*s[4]
                                          + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                          + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                          + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1];

                                     memcpy(&keys[cont * 8], &k, 8);
                                     cont++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

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

    int k, k2;

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

void tData_simple(char* src, char* dest, int typesize, int sub_shape[], int shape[], int dimension,
                  int inverse) {

    int MAX_DIM = 8;
    int DIM = dimension;

    int dim[MAX_DIM], shp[MAX_DIM], sub[MAX_DIM];

    for (int i = 0; i < MAX_DIM; i++) {
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

    int j;

    if (inverse == 0) {
        for (int k = 0; k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]; k++) {

            j = calculate_j(k, dim, shp, sub);


            memcpy(&dest[k * typesize], &src[j * typesize], typesize);
        }
    } else {
        for (int k = 0; k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]; k++) {

            j = calculate_j(k, dim, shp, sub);

            memcpy(&dest[j * typesize], &src[k * typesize], typesize);
        }
    }
}

int tData(char* src, char* dest, int typesize, int shape[], int pad_shape[],
          int sub_shape[], int size, int dimension, int inverse) {


    char* src_aux = (char *) calloc (size, typesize);

    padData(src, src_aux, typesize, shape, pad_shape, dimension);

    tData_simple(src_aux, dest, typesize, sub_shape, pad_shape, dimension, inverse);

    return 1;
}
''',
libraries=[])

ffibuilder.cdef(
                '''
                void tData(char* src, char* dest, int typesize, int shape[],
                            int pad_shape[], int sub_shape[], int size, int dimension,
                           int inverse);

                void tData_simple(char* src, char* dest, int typesize, int sub_shape[],
                                  int shape[], int dimension, int inverse);

                void padData(char* src, char* dest, int typesize, int shape[], int pad_shape[],
                             int dimension);


                void createIndexation(char* keys, int shape[], int sub_shape[],  int dimension);

                int calculate_j(int k, int dim[], int shp[], int sub[]);

                '''
                )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
