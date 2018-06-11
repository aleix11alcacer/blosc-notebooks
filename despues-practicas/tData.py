'''
Creates CFFI library to implement data transformation and compression algorithms.


calculate_j(k, dim, s, sb)

    Calculate the value of j and returns it.

    Parameters
    ----------
        k: int
        dim: int[]
        s: int[]
        sb: int[]

    Returns
    -------
        j: int


tData(src, dest, typesize, shape, pad_shape, sub_shape, size, dimension, inverse)

    Resize data of needed dimensions and transform it.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: int
            Data element size.
        shape: int[]
            Data shape.
        pad_shape: int[]
            Data transformed shape.
        sub_shape: int[]
            Data partition shape.
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


padData(src, dest, typesize, shape, pad_shape, dimension)

    Resize orignial data padding it with 0 to obtain expandend data  wich dimension is multiple
        of partition dimension.

    Parameters
    ----------
        src: char*
            Location of data to transform pointer.
        dest: char*
            Location of data transformed pointer.
        typesize: int
            Data element size.
        shape: int[]
            Original data shape.
        shape: int[]
            Expanded data shape.
        dimension: int
            Data dimension.


decompress_trans(comp, dest, trans_shape, part_shape, dimensions, sub_trans, dimension, b_size,
               typesize)

    Decompress partitioned data decompressing only the desired data.

    Parameters
    ----------
        comp: char*
            Location of compressed data pointer.
        dest: char*
            Location of data transformed pointer.
        trans_shape: int[]
            Partitioned data shape.
        part_shape: int[]
            A data partition shape.
        dimensions: int[]
            Defines data desired subset.
        sub_trans: int[]
            Desired data shape.
        dimension: int
            Data dimension.
        b_size: int
            Data partition size.
        typesize: int
            Data element size.
'''

from cffi import FFI

ffibuilder = FFI()

ffibuilder.set_source("tData",
                      '''
                      #include <stdio.h>
                      #include <stdint.h>
                      #include <blosc.h>
                      #include <time.h>
                      
                      int calculate_j(int k, int dim[], int s[], int sb[]) {
                      

                            int j = dim[0]*((k)%sb[0] + k/(sb[0]*sb[1]*sb[2]*sb[3]*sb[4]*sb[5]*sb[6]*sb[7])%(s[0]/sb[0])*sb[0])
                                +
                                dim[1]*(k/(sb[0])%sb[1]*s[0] + k/(s[0]*sb[1]*sb[2]*sb[3]*sb[4]*sb[5]*sb[6]*sb[7])%(s[1]/sb[1])*s[0]*sb[1])
                                +
                                dim[2]*(k/(sb[0]*sb[1])%sb[2]*s[1]*s[0] + k/(s[0]*s[1]*sb[2]*sb[3]*sb[4]*sb[5]*sb[6]*sb[7])%(s[2]/sb[2])*s[0]*s[1]*sb[2])
                                +
                                dim[3]*(k/(sb[0]*sb[1]*sb[2])%sb[3]*s[0]*s[1]*s[2] + k/(s[0]*s[1]*s[2]*sb[3]*sb[4]*sb[5]*sb[6]*sb[7])%(s[3]/sb[3])*s[0]*s[1]*s[2]*sb[3])
                                +
                                dim[4]*(k/(sb[0]*sb[1]*sb[2]*sb[3])%sb[4]*s[0]*s[1]*s[2]*s[3] + k/(s[0]*s[1]*s[2]*s[3]*sb[4]*sb[5]*sb[6]*sb[7])%(s[4]/sb[4])*s[0]*s[1]*s[2]*s[3]*sb[4])
                                +
                                dim[5]*(k/(sb[0]*sb[1]*sb[2]*sb[3]*sb[4])%sb[5]*s[0]*s[1]*s[2]*s[3]*s[4] + k/(s[0]*s[1]*s[2]*s[3]*s[4]*sb[5]*sb[6]*sb[7])%(s[5]/sb[5])*s[0]*s[1]*s[2]*s[3]*s[4]*sb[5])
                                +
                                dim[6]*(k/(sb[0]*sb[1]*sb[2]*sb[3]*sb[4]*sb[5])%sb[6]*s[0]*s[1]*s[2]*s[3]*s[4]*s[5] + k/(s[0]*s[1]*s[2]*s[3]*s[4]*s[5]*sb[6]*sb[7])%(s[6]/sb[6])*s[0]*s[1]*s[2]*s[3]*s[4]*s[5]*sb[6])
                                +
                                dim[7]*(k/(sb[0]*sb[1]*sb[2]*sb[3]*sb[4]*sb[5]*sb[6])%sb[7]*s[0]*s[1]*s[2]*s[3]*s[4]*s[5]*s[6] + k/(s[0]*s[1]*s[2]*s[3]*s[4]*s[5]*s[6]*sb[7])%(s[7]/sb[7])*s[0]*s[1]*s[2]*s[3]*s[4]*s[5]*s[6]*sb[7]);
                        
                                return j;
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
                                       
                          int k, j;
                                                   
                          if (inverse == 0) {
                              for (int a = 0; a < s[0]; a++) {
                                  for (int b = 0; b < s[1]; b++) {
                                      for (int c = 0; c < s[2]; c++) {
                                          for (int d = 0; d < s[3]; d++) {
                                              for (int e = 0; e < s[4]; e++) {
                                                  for (int f = 0; f < s[5]; f++) {
                                                      for (int g = 0; g < s[6]; g++) {
                                                          for (int h = 0; h < s[7]; h+=sb[7]) {
                      
                                                              k =  h
                                                                   + g*s[7]
                                                                   + f*s[7]*s[6]
                                                                   + e*s[7]*s[6]*s[5]
                                                                   + d*s[7]*s[6]*s[5]*s[4]
                                                                   + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                                                   + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                                                   + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1];
                                                                                                                                 
                                                              j = calculate_j(k, dim, shp, sub);
                      
                                                              memcpy(&dest[k * typesize], &src[j * typesize], typesize * sb[7]);
                      
                      
                                                          }
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                  }
                              }
                      
                          } else {
                              for (int a = 0; a < s[0]; a++) {
                                  for (int b = 0; b < s[1]; b++) {
                                      for (int c = 0; c < s[2]; c++) {
                                          for (int d = 0; d < s[3]; d++) {
                                              for (int e = 0; e < s[4]; e++) {
                                                  for (int f = 0; f < s[5]; f++) {
                                                      for (int g = 0; g < s[6]; g++) {
                                                          for (int h = 0; h < s[7]; h+=sb[7]) {
                      
                                                              k =  h
                                                                   + g*s[7]
                                                                   + f*s[7]*s[6]
                                                                   + e*s[7]*s[6]*s[5]
                                                                   + d*s[7]*s[6]*s[5]*s[4]
                                                                   + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                                                   + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                                                   + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1];
                                                                   
                                                              printf("%d - ", h);
                                                              
                                                              j = calculate_j(k, dim, s, sb);
                      
                                                              memcpy(&dest[j * typesize], &src[k * typesize], typesize * sb[7]);
                      
                      
                                                          }
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                  }
                              }
                          }
                      }
                      
                      void tData(char* src, char* dest, int typesize, int shape[], int pad_shape[],
                                int sub_shape[], int size, int dimension, int inverse) {
                      
                      
                          char* src_aux = (char *) calloc (size, typesize);
                      
                          padData(src, src_aux, typesize, shape, pad_shape, dimension);
                      
                          tData_simple(src_aux, dest, typesize, sub_shape, pad_shape, dimension, inverse);
                      
                      }
                      
                      void decompress_trans(char* comp, char* dest2, int shape[], int trans_shape[], int part_shape[],
                                            int final_s[], int dimensions[], int dimension, int b_size, int typesize){
                      
                          int MAX_DIM = 8;
                          int DIM = dimension;
                      
                          // Calculate dimensions data
                      
                          int subpl[MAX_DIM], ts[MAX_DIM], ps[MAX_DIM], dim[MAX_DIM], sd[MAX_DIM], fs[MAX_DIM], s[MAX_DIM];
                      
                          for (int i = 0; i < MAX_DIM; i++) {
                              if (i < DIM) {
                                  s[MAX_DIM + i - DIM] = shape[i];
                                  ts[MAX_DIM + i - DIM] = trans_shape[i];
                                  ps[MAX_DIM + i - DIM] = part_shape[i];
                                  dim[MAX_DIM + i - DIM] = dimensions[i];
                                  fs[MAX_DIM + i - DIM] = final_s[i];
                                  sd[MAX_DIM + i - DIM] = trans_shape[i]/part_shape[i];
                      
                              } else {
                                  s[MAX_DIM - i - 1] = 1;
                                  ts[MAX_DIM - i - 1] = 1;
                                  ps[MAX_DIM - i - 1] = 1;
                                  dim[MAX_DIM - i - 1] = -1;
                                  sd[MAX_DIM - i - 1] = 1;
                                  fs[MAX_DIM - i - 1] = 1;
                      
                              }
                      
                          }
                      
                          int oi_s[MAX_DIM], oi_f[MAX_DIM];
                      
                      
                          for (int i = 0; i < MAX_DIM; i++) {
                              if (dim[i] != -1) {
                                  subpl[i] = ps[i];
                                  oi_s[i] = dim[i];
                                  oi_f[i] = dim[i] + 1;
                              } else {
                                  subpl[i] = ts[i];
                                  oi_s[i] = 0;
                                  oi_f[i] = ts[i];
                              }
                          }
                      
                          int subpl_size = 1;
                      
                          for (int i = 0; i < MAX_DIM; i++) {
                              subpl_size *= subpl[i];
                          }
                      
                          // Malloc buffers
                      
                          char *aux = malloc(b_size * typesize);
                          char *dest = malloc(subpl_size * typesize);
                      
                          // Calculate index to decompress
                      
                          int k, n;
                          int h2, g2, f2, e2, d2, c2, b2, a2;
                          int cont;
                      
                          for (int a = oi_s[0]/ps[0]*ps[0]; a < oi_f[0]; a += ps[0]) {
                              for (int b = oi_s[1]/ps[1]*ps[1]; b < oi_f[1]; b += ps[1]) {
                                  for (int c = oi_s[2]/ps[2]*ps[2]; c < oi_f[2]; c += ps[2]) {
                                      for (int d = oi_s[3]/ps[3]*ps[3]; d < oi_f[3]; d += ps[3]) {
                                          for (int e = oi_s[4]/ps[4]*ps[4]; e < oi_f[4]; e += ps[4]) {
                                              for (int f = oi_s[5]/ps[5]*ps[5]; f < oi_f[5]; f += ps[5]) {
                                                  for (int g = oi_s[6]/ps[6]*ps[6]; g < oi_f[6]; g += ps[6]) {
                                                      for (int h = oi_s[7]/ps[7]*ps[7]; h < oi_f[7]; h += ps[7]) {
                      
                                                          k = h
                                                               + g*ts[7]
                                                               + f*ts[7]*ts[6]
                                                               + e*ts[7]*ts[6]*ts[5]
                                                               + d*ts[7]*ts[6]*ts[5]*ts[4]
                                                               + c*ts[7]*ts[6]*ts[5]*ts[4]*ts[3]
                                                               + b*ts[7]*ts[6]*ts[5]*ts[4]*ts[3]*ts[2]
                                                               + a*ts[7]*ts[6]*ts[5]*ts[4]*ts[3]*ts[2]*ts[1];
                      
                                                           n = a/ps[0]*sd[1]*sd[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                                              + b/ps[1]*sd[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                                              + c/ps[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                                              + d/ps[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                                              + e/ps[4]*sd[5]*sd[6]*sd[7]
                                                              + f/ps[5]*sd[6]*sd[7]
                                                              + g/ps[6]*sd[7]
                                                              + h/ps[7];
                      
                      
                                                          blosc_getitem(comp, n * b_size, b_size, aux);
                      
                                                          h2 = k % ts[7] % subpl[7];
                                                          g2 = k / (ts[7]) % subpl[6];
                                                          f2 = k / (ts[7]*ts[6]) % subpl[5];
                                                          e2 = k / (ts[7]*ts[6]*ts[5]) % subpl[4];
                                                          d2 = k / (ts[7]*ts[6]*ts[5]*ts[4]) % subpl[3];
                                                          c2 = k / (ts[7]*ts[6]*ts[5]*ts[4]*ts[3]) % subpl[2];
                                                          b2 = k / (ts[7]*ts[6]*ts[5]*ts[4]*ts[3]*ts[2]) % subpl[1];
                                                          a2 = k / (ts[7]*ts[6]*ts[5]*ts[4]*ts[3]*ts[2]*ts[1]) % subpl[0];
                      
                                                          // Copy block to final data
                      
                                                          cont = 0;
                      
                                                          for (int ra = a2; ra < a2 + ps[0]; ra++) {
                                                              for (int rb = b2; rb < b2 + ps[1]; rb++) {
                                                                  for (int rc = c2; rc < c2 + ps[2]; rc++) {
                                                                      for (int rd = d2; rd < d2 + ps[3]; rd++) {
                                                                          for (int re = e2; re < e2 + ps[4]; re++) {
                                                                              for (int rf = f2; rf < f2 + ps[5]; rf++) {
                                                                                  for (int rg = g2; rg < g2 + ps[6]; rg++) {
                      
                                                                                      k = h2
                                                                                           + rg * subpl[7]
                                                                                           + rf * subpl[7]*subpl[6]
                                                                                           + re * subpl[7]*subpl[6]*subpl[5]
                                                                                           + rd * subpl[7]*subpl[6]*subpl[5]*subpl[4]
                                                                                           + rc * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]
                                                                                           + rb * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]
                                                                                           + ra * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]*subpl[1];
                      
                                                                                      memcpy(&dest[k * typesize], &aux[cont * typesize], ps[7] * typesize);
                                                                                      cont += ps[7];
                                                                                  }
                                                                              }
                                                                          }
                                                                      }
                                                                  }
                                                              }
                                                          }
                      
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                  }
                             }
                      
                      
                          }
                      
                          // Reduce final data
                      
                          int ini[MAX_DIM], fin[MAX_DIM];
                      
                          for (int i = 0; i < MAX_DIM - 1; i += 1) {
                              if (dim[i] != -1) {
                                  ini[i] = dim[i] % ps[i];
                                  fin[i] = dim[i] % ps[i] + 1;
                              }
                              else {
                                  ini[i] = 0;
                                  fin[i] = s[i];
                              }
                          }
                      
                          if (dim[MAX_DIM - 1] != -1) {
                              ini[MAX_DIM - 1] = dim[MAX_DIM - 1] % ps[MAX_DIM - 1];
                              fin[MAX_DIM - 1] = dim[MAX_DIM - 1] % ps[MAX_DIM - 1] + 1 ;
                          } else {
                              ini[MAX_DIM - 1] = 0;
                              fin[MAX_DIM - 1] = 1;
                          }
                      
                          cont = 0;
                      
                          for (int a = ini[0]; a < fin[0]; a += 1) {
                              for (int b = ini[1]; b < fin[1]; b += 1) {
                                  for (int c = ini[2]; c < fin[2]; c += 1) {
                                      for (int d = ini[3]; d < fin[3]; d += 1) {
                                          for (int e = ini[4]; e < fin[4]; e += 1) {
                                              for (int f = ini[5]; f < fin[5]; f += 1) {
                                                  for (int g = ini[6]; g < fin[6]; g += 1) {
                                                      for (int h = ini[7]; h < fin[7]; h += 1) {
                                                          k = h
                                                               + g * subpl[7]
                                                               + f * subpl[7]*subpl[6]
                                                               + e * subpl[7]*subpl[6]*subpl[5]
                                                               + d * subpl[7]*subpl[6]*subpl[5]*subpl[4]
                                                               + c * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]
                                                               + b * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]
                                                               + a * subpl[7]*subpl[6]*subpl[5]*subpl[4]*subpl[3]*subpl[2]*subpl[1];
                      
                                                          memcpy(&dest2[cont * typesize], &dest[k * typesize], fs[7] * typesize);
                      
                                                          cont += fs[7];
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                  }
                              }
                          }
                          free(aux);
                          free(dest);
                      }
                      
                      ''', libraries=['blosc'])

ffibuilder.cdef(
    '''
    void tData(char* src, char* dest, int typesize, int shape[],
                int pad_shape[], int sub_shape[], int size, int dimension,
               int inverse);

    void tData_simple(char* src, char* dest, int typesize, int sub_shape[],
                      int shape[], int dimension, int inverse);

    void padData(char* src, char* dest, int typesize, int shape[], int pad_shape[],
                 int dimension);


    int calculate_j(int k, int dim[], int s[], int sb[]);

    void decompress_trans(char* comp, char* dest2, int shape[], int trans_shape[],
                          int part_shape[], int final_s[], int dimensions[],
                          int dimension, int b_size, int typesize);

    '''
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
