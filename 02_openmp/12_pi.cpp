#include <cstdio>
#include <omp.h>

int main() {
  int n = 10;
  double dx = 1. / n;
  double pi = 0;
  double pilist[n];
  double piresult[n];

#pragma omp parallel for reduction(+:pi) 
  for(int i=0; i<n; i++){
    double x = (i + 0.5) * dx;
    pi += 4.0 / (1.0 + x * x) * dx;
  }
  printf("%17.15f\n", pi);


/* Using prefix sum:
#pragma omp parallel for
  for (int i=0; i<n; i++) {
    double x = (i + 0.5) * dx;
    pilist[i] = 4.0 / (1.0 + x * x) * dx;
  }
  printf("\n");

#pragma omp parallel
  for(int j=1; j<n; j <<=1){

#pragma omp for
    for(int i =0; i<n; i++)
	piresult[i] = pilist[i];
#pragma omp for
    for(int i=j; i<n; i++){
	pilist[i] += piresult[i-j];
    }
  }
  printf("%17.15f\n",pilist[n-1]);
*/
}
