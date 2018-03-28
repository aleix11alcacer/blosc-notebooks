'''
Crea una libreria que ciene en una función los dos algoritmos juntos: el simple y el general.
'''

from cffi import FFI
ffibuilder = FFI()

ffibuilder.set_source("tData0",
'''
#include <stdio.h>
#include <string.h>

int calculate_j(size_t k, size_t dim[], size_t shp[], size_t sub[]) {

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

int calculate_i(size_t j, size_t shp[], size_t sho[], size_t sub[]) {
    size_t i = j/shp[0] * (sho[0]%sub[0])
                 + j/(shp[1] * shp[0]) * (sho[1]%sub[1])*sho[0]
                 + j/(shp[2] * shp[1] * shp[0]) * (sho[2]%sub[2]) * sho[0] * sho[1]
                 + j/(shp[3] * shp[2]*shp[1]*shp[0]) * (sho[3]%sub[3]) * sho[0]*sho[1]*sho[2]
                 + j/(shp[4] * shp[3]*shp[2]*shp[1]*shp[0]) * (sho[4]%sub[4]) * sho[0]*sho[1]*sho[2]*sho[3]
                 + j/(shp[5] * shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[5]%sub[5]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]
                 + j/(shp[6] * shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[6]%sub[6]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]
                 + j/(shp[7] * shp[6]*shp[5]*shp[4]*shp[3]*shp[2]*shp[1]*shp[0]) * (sho[7]%sub[7]) * sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6];
     return i;

}

int calculate_cond(size_t k, size_t c, size_t sho[], size_t sub[]) {
    int cond = (((k + c)%sho[0] >= sho[0] - sho[0]%sub[0])
            || ((k + c)%(sho[0]*sho[1])/(sho[0]) >= sho[1] - sho[1]%sub[1])
            || ((k + c)%(sho[0]*sho[1]*sho[2])/(sho[0]*sho[1]) >= sho[2] - sho[2]%sub[2])
            || ((k + c)%(sho[0]*sho[1]*sho[2]*sho[3])/(sho[0]*sho[1]*sho[2]) >= sho[3] - sho[3]%sub[3])
            || ((k + c)%(sho[0]*sho[1]*sho[2]*sho[3]*sho[4])/(sho[0]*sho[1]*sho[2]*sho[3]) >= sho[4] - sho[4]%sub[4])
            || ((k + c)%(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5])/(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]) >= sho[5] - sho[5]%sub[5])
            || ((k + c)%(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6])/(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]) >= sho[6] - sho[6]%sub[6])
            || ((k + c)%(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7])/(sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]) >= sho[7] - sho[7]%sub[7]));
    return cond;
}

void tData_general(void* src, void* dest, size_t typesize, size_t sub_shape[], size_t shape[], size_t dimension, size_t inverse) {

    size_t MAX_DIM = 8;
    size_t DIM = dimension;


    size_t dim[MAX_DIM], sho[MAX_DIM], shp[MAX_DIM], sub[MAX_DIM];
    size_t i;

    for (size_t i = 0; i < MAX_DIM; i++) {

        if (i < DIM) {
            dim[i] = 1;
            sho[DIM - i - 1] = shape[i];
            sub[DIM - i - 1] = sub_shape[i];
            shp[DIM - i - 1] = shape[i]/sub_shape[i]*sub_shape[i];

        } else {

            dim[i] = 0;
            shp[i] = 1;
            sub[i] = 1;
            sho[i] = 1;

        }
    }

    size_t c = 0;
    size_t j = 0;
    size_t k = 0;
    size_t l = 0;

    if (inverse == 0) {
        while (k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]){

            if (calculate_cond(k, c, sho, sub)){

                    memcpy(&dest[(k + c) * typesize], &src[(k + c) * typesize], typesize);
                    c = c + 1;

            } else {

                l = k + c;

                j = calculate_j(k, dim, shp, sub);

                i = calculate_i(j, shp, sho, sub);

                memcpy(&dest[l * typesize], &src[(j + i) * typesize], typesize);
                k += 1;
            }
        }
        for (size_t l = k + c; l < sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7]; l++) {
            memcpy(&dest[l * typesize], &src[l * typesize], typesize);
        }
    }
    else {
        while (k < shp[0]*shp[1]*shp[2]*shp[3]*shp[4]*shp[5]*shp[6]*shp[7]){

            if (calculate_cond(k, c, sho, sub)){

                    memcpy(&dest[(k + c) * typesize], &src[(k + c) * typesize], typesize);
                    c = c + 1;

            } else {

                l = k + c;

                j = calculate_j(k, dim, shp, sub);

                i = calculate_i(j, shp, sho, sub);

                memcpy(&dest[(j + i) * typesize], &src[l * typesize], typesize);
                k += 1;
            }
        }
        for (size_t l = k + c; l < sho[0]*sho[1]*sho[2]*sho[3]*sho[4]*sho[5]*sho[6]*sho[7]; l++) {
            memcpy(&dest[l * typesize], &src[l * typesize], typesize);
        }
    }
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

    size_t DIM = dimension;

    for (size_t i = 0; i < DIM; i++) {
        if (shape[i]%sub_shape[i] != 0){
            tData_general(src, dest, typesize, sub_shape, shape, dimension, inverse);
            return 0;
        }
    }

    tData_simple(src, dest, typesize, sub_shape, shape, dimension, inverse);
    return 1;
}

''',
libraries=[])

ffibuilder.cdef(
                '''
                void tData(char* src, char* dest, size_t typesize, size_t sub_shape[],
                           size_t shape[], size_t dimension, size_t inverse);

                void tData_general(char* src, char* dest, size_t typesize, size_t sub_shape[],
                                   size_t shape[], size_t dimension, size_t inverse);

                void tData_simple(char* src, char* dest, size_t typesize, size_t sub_shape[],
                                  size_t shape[], size_t dimension, size_t inverse);

                '''
                )

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
