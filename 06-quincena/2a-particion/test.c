#include <stdlib.h>
#include <stdio.h>


int main( int argc, char *argv[ ]) {
    int s[8], sb[8], sd[8];
    for (int i = 0; i < 5; i++) {
        s[i] = 1;
        sb[i] = 1;
        sd[i] = 1;
    }
    for (int i = 5; i < 8; i++) {
        s[i] = 16;
        sb[i] = 4;
        sd[i] = s[i]/sb[i];
    }
    sb[6] = 5;
    s[6] = 20;


    int co = 0;
    for (int a=0; a < s[0]; a += sb[0]){
        for (int b=0; b < s[1]; b += sb[1]){
            for (int c=0; c < s[2]; c += sb[2]){
                for (int d=0; d < s[3]; d += sb[3]){
                    for (int e=0; e < s[4]; e += sb[4]){
                        for (int f=0; f < s[5]; f += sb[5]){
                            for (int g=0; g < s[6]; g += sb[6]){
                                for (int h=0; h < s[7]; h += sb[7]){

                                     int k = h
                                          + g*s[7]
                                          + f*s[7]*s[6]
                                          + e*s[7]*s[6]*s[5]
                                          + d*s[7]*s[6]*s[5]*s[4]
                                          + c*s[7]*s[6]*s[5]*s[4]*s[3]
                                          + b*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]
                                          + a*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1];

                                     int cont = a/sb[0]*sd[1]*sd[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                            + b/sb[1]*sd[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                            + c/sb[2]*sd[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                            + d/sb[3]*sd[4]*sd[5]*sd[6]*sd[7]
                                            + e/sb[4]*sd[5]*sd[6]*sd[7]
                                            + f/sb[5]*sd[6]*sd[7]
                                            + g/sb[6]*sd[7]
                                            + h/sb[7];

                                    printf("%d - %d\n", cont, co);
                                    co++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
