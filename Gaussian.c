#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void){

  FILE *fp;
  fp = fopen("xnData.dat", "w"); // Opens a file into which data will be written.

  //srand(0);
  long time();
  srand(time(0));

  float i = 0, N = 0, ND = 0, j = 0, acceptance = 0;
  double x0 = 0, x1 = 0, x = 0, y = 0, x_squared = 0;
  double alpha = 0, beta = 0, r = rand()/RAND_MAX, P = 0;

  alpha = 0.3; // 0.3
  beta = 0.4; // 0.1 -> 0.6
  x0 = 1; // 1

  printf("Input desired number of sweeps:\n");  // 100 000
  scanf("%f", &N);
  printf("Input desired number of discards:\n");  // 10 000
  scanf("%f", &ND);

  int SampleSize = (N - ND);
  for(i = 1; i <= SampleSize; i++) {

    y = x0 - (2 * alpha * r) + alpha; // candidate sample value given the previous sample x0
    P = exp( (-1 * beta) * ((y*y) - (x*x)) ); // proportional to a Gaussian distribution

    if(P > 1) { // if acceptance = P(y)/P(x0) = P(y) is >1 (more probable), then candidate value is immediately accepted
      x1 = y;
      j += 1; // j is an acceptance counter, increments when a new value is accepted by the algorithm.
    }
    else {  // generate r in [0,1], if acceptance = P(y) is >r then accept, else reject and return to inital sample value x0
      r = (double) (rand()) / (RAND_MAX);
      if(P >= r) {
        x1 = y;
        j += 1;
      }
      else {
        x1 = x0;
      }
    }

    x_squared = (x1)*(x1);
    fprintf(fp, "%f\n", x_squared);
    x0 = x1; // assigns the variable x0 to its new value before the loop repeats

  }

  acceptance = j/(N - ND); // acceptance is calculated and then printed to the user, typically aim for between 0.5 and 0.8.
  printf("Acceptance = %f", acceptance);
  fclose(fp);

return 0;
}
