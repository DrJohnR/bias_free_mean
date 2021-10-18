#include <stdio.h>
#include <math.h>

#define SampleSize 90000 // = (N - ND) from Gaussian.c

int main(void) {

  FILE *fp = fopen("xnData.dat", "r"); // opens the data file created by Gaussian.c
  FILE *fp1 = fopen("HalfSample1.dat", "w"); // two files opened, each will contain half of the data set.
  FILE *fp2 = fopen("HalfSample2.dat", "w"); // if total number of data elements is odd, this file will contain one extra element.

  int i = 0;
  double x[SampleSize];

 for(i = 0; i < SampleSize; i++) {
   fscanf(fp, "%lf", &x[i]);  // import the complete data set from Gaussian.c
 }

 for(i = 0; i < (SampleSize / 2); i++) {
    fprintf(fp1, "%f\n", x[i]);   // first half of data written to HalfSample1.dat
  }
  for(i = (SampleSize / 2); i < SampleSize; i++) {
    fprintf(fp2, "%f\n", x[i]);   // second half of data written to HalfSample2.dat
  }

  fclose(fp);  // files are closed.
  fclose(fp1);
  fclose(fp2);

  return 0;
}
