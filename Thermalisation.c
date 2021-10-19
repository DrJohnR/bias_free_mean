#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SampleSize 90000 // = (N - ND) from Gaussian.c

int main(void) {

  FILE *fp1 = fopen("HalfSample1.dat", "r");
  FILE *fp2 = fopen("HalfSample2.dat", "r");

  int i = 0, halfSample = SampleSize/2;
  double mean1 = 0, mean2 = 0, variance1 = 0, variance2 = 0, stand_dev1 = 0, stand_dev2 = 0;
  double x1[halfSample], x2[halfSample];

  for(i = 0; i < halfSample; i++) {
    fscanf(fp1, "%lf", &x1[i]);
    fscanf(fp2, "%lf", &x2[i]);

    mean1 += x1[i];
    mean2 += x2[i];
 }
  mean1 /= halfSample;
  mean2 /= halfSample;

  for(i = 0; i < halfSample; i++) {
    variance1 += (mean1 - x1[i])*(mean1 - x1[i]);
    variance2 += (mean2 - x2[i])*(mean2 - x2[i]);
  }
  variance1 /= halfSample;
  variance2 /= halfSample;
  stand_dev1 = sqrt(variance1);
  stand_dev2 = sqrt(variance2);

 fclose(fp1);
 fclose(fp2);

 printf("----- Thermalisation statistics for the two data subsets -----\n");
 printf("Mean_1: %lf  Mean_2: %lf\n", mean1, mean2);
 printf("Variance_1: %lf  Variance_2: %lf\n", variance1, variance2);
 printf("Stand_Dev_1: %lf  Stand_Dev_2: %lf\n\n", stand_dev1, stand_dev2);

 if( (mean1 - stand_dev1) <= mean2 && mean2 <= (mean1 + stand_dev1) &&
      (mean2 - stand_dev2) <= mean1 && mean1 <= (mean2 + stand_dev2)  ){
   printf("Data is adequately thermalised.");
 }
 else{
   printf("Data is not adequately thermalised.");
 }

  return 0;
}
