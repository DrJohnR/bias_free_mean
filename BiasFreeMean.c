#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void){

long time();
srand(time(0));


/*//////////////// (1) Generate the data set ////////////////*/
FILE *fp; // create a file pointer
fp = fopen("xnData.dat", "w"); // the data set will be exported to this file.

char therm;
float i = 0, N = 0, ND = 0, j = 0, acceptance = 0;
double x0 = 0, x1 = 0, x = 0, y = 0, x_squared = 0;
double alpha = 0, beta = 0, r = rand()/RAND_MAX, P = 0;

alpha = 0.3; // 0.3
beta = 0.6; // 0.1 -> 0.6
x0 = 1; // 1

printf("Input desired number of sweeps:\n");  // 100 000
scanf("%f", &N);
printf("Input desired number of discards:\n");  // 10 000
scanf("%f", &ND);
printf("Impose data thermalisation? (y/n)\n");
scanf(" %c", &therm);
if(therm != 'y' && therm != 'n'){
  printf("Unrecognised command, type y for yes, or n for no. Terminating.");
  return 1;
}

int SampleSize = (N - ND);
if( (int)(N - ND)%2 != 0 && therm == 'y'){
  printf("Warning:\n data set contains an odd number of elements, thermalisation using two data subsets may be less accurate.\n\n");
}
  
for(i = 1; i <= SampleSize; i++) {
  y = x0 - (2 * alpha * r) + alpha; // candidate sample value given the previous sample x0
  P = exp( (-1 * beta) * ((y*y) - (x*x)) ); // proportional to a Gaussian distribution

  if(P > 1) { // if acceptance = P(y)/P(x0) = P(y) is >1 (more probable), then candidate value is immediately accepted
    x1 = y;
    j += 1; // j is an acceptance counter, increments when a new value is accepted by the algorithm.
  }
  else { // generate r in [0,1], if acceptance = P(y) is >r then accept, else reject and return to inital sample value x0
    r = (double) (rand()) / (RAND_MAX);
    if(P >= r) {
      x1 = y;
      j += 1;
    }
    else {
      x1 = x0;
    }
  }

  x_squared = x1*x1;
  fprintf(fp, "%f\n", x_squared);
  x0 = x1; // assigns the variable x0 to its new value before the loop repeats

} // end for loop on line 40
  
  acceptance = j/(N - ND); // acceptance is calculated and then printed to the user, typically aim for between 0.5 and 0.8.
  printf("Acceptance = %f\n\n", acceptance);
  fclose(fp);
  /*///////////////////////////////////////////////////////*/



/*//////////////// (2) Subroutine runs if thermalisation is requested ////////////////*/
if(therm == 'y'){

/*//////// (2a) Dividing the data set ////////*/
FILE *fp = fopen("xnData.dat", "r"); // opens the data file created in (1)
FILE *fp1 = fopen("HalfSample1.dat", "w"); // two files opened, each will contain half of the data set.
FILE *fp2 = fopen("HalfSample2.dat", "w"); // if total number of data elements is odd, this file will contain one extra element.

int SampleSize = (N - ND), i = 0;
double x[SampleSize];

for(i = 0; i < SampleSize; i++) {
  fscanf(fp, "%lf", &x[i]);  // import the complete data set
}
for(i = 0; i < (SampleSize / 2); i++) {
   fprintf(fp1, "%f\n", x[i]);  // first half of data written to HalfSample1.dat
 }
 for(i = (SampleSize / 2); i < SampleSize; i++) {
   fprintf(fp2, "%f\n", x[i]);  // second half of data written to Halfsample2.dat
 }

 fclose(fp);  // files are closed.
 fclose(fp1);
 fclose(fp2);
  

/*//////// (2b) Computing thermalisation statistics ////////*/
 fp1 = fopen("HalfSample1.dat", "r");
 fp2 = fopen("HalfSample2.dat", "r");

 double mean1 = 0, mean2 = 0, variance1 = 0, variance2 = 0, stand_dev1 = 0, stand_dev2 = 0;
 int halfSample = SampleSize/2;
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

if( (mean1 - stand_dev1) <= mean2 <= (mean1 + stand_dev1) && (mean2 - stand_dev2) <= mean1 <= (mean2 + stand_dev2) ){
  printf("Data is adequately thermalised, proceeding...\n\n");
}
else{
  printf("Data is inadequately thermalised, try altering the parameters alpha and beta. Terminating.");
  return 1;
}

} // end of optional thermalisation subroutine
/*///////////////////////////////////////////////////////*/



/*//////////////// (3) Binning and resampling data set to compute the bias-free mean ////////////////*/
fp = fopen("xnData.dat", "r");

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
} // end for loop on line 162

bf_mean /= max_it;
printf("Bias-free mean of data sample after %d resamples using %d bins is %f\n", max_it, bins, bf_mean);
fclose(fp);


return 0;
} // end of main()
