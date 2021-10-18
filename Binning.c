#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SampleSize 90000 // = (N - ND) from Gaussian.c

int main(void) {

  FILE *fp = fopen("xnData.dat", "r");

  //srand(0);
  long time();
  srand(time(0));

  int dat = 0, b = 1, bins = 0, bin_no = 0, it = 0, max_it = 0;
  float bin_sum = 0, boot_sum = 0, boot_mean = 0, bf_mean = 0; // bf --> bias-free, boot --> bootstrap
  double arr[SampleSize];

  printf("Input desired number of data bins:\n");
  scanf("%d", &bins);
  if(SampleSize % bins != 0) {
    printf("Ensure that number of bins exactly divides the number of data points. Terminating.");
    return 1;
  }
  printf("Input desired number of resamples (iterations):\n");
  scanf("%d", &max_it);


  for(it = 1; it < max_it; it++) {
    boot_sum = 0, boot_mean = 0; // reinitialise before next iteration
    for(b = 1; b <= bins; b++) {
      bin_sum = 0; // empties the bin before adding to it again
      bin_no = 1 + (rand()%bins);  // randomly generates an integer in [1, 2, ..., bins-1, bins]

      for(dat = ((bin_no - 1) * SampleSize / bins); dat < (bin_no * SampleSize / bins); dat++) {  // Loop runs over data within bin range
        fscanf(fp, "%lf", &arr[dat]);  // data set elements within bin range are read from .dat file
        bin_sum += arr[dat];  // data set elements within the bin range are summed and added to bin_sum
      } // end of bin range

    boot_sum += bin_sum; // total sum of the bin contributes to the boot_sum
  } // loop runs until requested number of bin sums have been added to boot_sum

  boot_mean = boot_sum / SampleSize;  // computes the mean of each resample
  bf_mean += boot_mean; // sums each of the resample means
  }

  bf_mean /= max_it;
  printf("Bias-free mean of data sample after %d resamples using %d bins is %f\n", max_it, bins, bf_mean);
  fclose(fp);

  return 0;
}
